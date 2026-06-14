/*
orcwa_datalib.c
OpenRCWA データ作成ライブラリ (C) (OpenFDTD の ofd_datalib に準拠)
python/datalib/orcwa_datalib.py と同じ書式の入力データファイルを出力する
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "orcwa_datalib.h"

#define VERSION "OpenRCWA 4 2"

#define MAXSECTION   100
#define MAXMATERIAL  256
#define MAXGEOMETRY 1000
#define MAXFEED      100
#define MAXLOAD      100
#define MAXPOINT     100
#define MAXFAR1D     100
#define MAXNEAR1D    100
#define MAXNEAR2D    100

#define MAXTITLE     256
#define MAXPROP       16

/* === 計算 === */

static char   Title[MAXTITLE];

static int    NXsection, NYsection, NZsection;
static double Xsection[MAXSECTION], Ysection[MAXSECTION], Zsection[MAXSECTION];
static int    NXdivision, NYdivision, NZdivision;
static int    Xdivision[MAXSECTION], Ydivision[MAXSECTION], Zdivision[MAXSECTION];

typedef struct {
	int    type;                       /* 1: material, 2: material_dispersion */
	double d[4];
} material_t;
static int        NMaterial;
static material_t MaterialData[MAXMATERIAL];

typedef struct {
	int    material;
	int    shape;
	int    ng;                         /* 6 or 8 */
	double g[8];
} geometry_t;
static int        NGeometry;
static geometry_t GeometryData[MAXGEOMETRY];

typedef struct {
	char   dir;
	double x, y, z, amp, delay, z0;
} feed_t;
static int    NFeed;
static feed_t FeedData[MAXFEED];

static int    IPlanewave;
static double PlanewaveTheta, PlanewavePhi;
static int    PlanewavePol;

typedef struct {
	char   dir;
	double x, y, z;
	char   type;
	double rcl;
} load_t;
static int    NLoad;
static load_t LoadData[MAXLOAD];

typedef struct {
	char   dir;
	double x, y, z;
	char   prop[MAXPROP];
} point_t;
static int     NPoint;
static point_t PointData[MAXPOINT];

static int    IRfeed,      IPulsewidth,    ITimestep;
static double RfeedValue,  PulsewidthValue, TimestepValue;

static int    IPml;
static int    PmlL;
static double PmlM, PmlR0;

static int    IPbc;
static int    PbcX, PbcY, PbcZ;

static int    IFreq1, IFreq2;
static double Freq1[2], Freq2[2];
static int    Freq1Div, Freq2Div;

static int    ISolver;
static int    SolverMaxiter, SolverNout;
static double SolverConverg;

/* === ポスト処理 === */

static int    IMatchingloss, Matchingloss;
static int    IPlotiter,     Plotiter;
static int    IPlotfeed,     Plotfeed;
static int    IPlotpoint,    Plotpoint;
static int    IPlotsmith,    Plotsmith;

typedef struct {
	int    scale;
	double smin, smax;
	int    sdiv;
} plotscale_t;
static int         IPlotzin,      IPlotyin,      IPlotref,      IPlotspara,      IPlotcoupling;
static plotscale_t Plotzin,       Plotyin,       Plotref,       Plotspara,       Plotcoupling;

static int         IPlotfar0d;
static double      Far0dTheta, Far0dPhi;
static plotscale_t Plotfar0d;

static int    IFreqdiv, Freqdiv;

typedef struct {
	char   dir;
	int    div;
	double angle;
} far1d_t;
static int     NFar1d;
static far1d_t Far1dData[MAXFAR1D];

static int    IFar1dstyle,     Far1dstyle;
static int    IFar1dcomponent, Far1dcomponent[3];
static int    IFar1ddb,        Far1ddb;
static int    IFar1dnorm,      Far1dnorm;
static int    IFar1dscale;
static double Far1dscale[2];
static int    Far1dscaleDiv;

static int    IPlotfar2d, Far2dDivtheta, Far2dDivphi;
static int    IFar2dcomponent, Far2dcomponent[7];
static int    IFar2ddb,        Far2ddb;
static int    IFar2dscale;
static double Far2dscale[2];
static int    Far2dscaleDiv;
static int    IFar2dobj;
static double Far2dobj;

typedef struct {
	char   component[MAXPROP];
	char   dir;
	double p1, p2;
} near1d_t;
static int      NNear1d;
static near1d_t Near1dData[MAXNEAR1D];

static int    INear1ddb,    Near1ddb;
static int    INear1dnoinc, Near1dnoinc;
static int    INear1dscale;
static double Near1dscale[2];
static int    Near1dscaleDiv;

typedef struct {
	char   component[MAXPROP];
	char   dir;
	double p;
} near2d_t;
static int      NNear2d;
static near2d_t Near2dData[MAXNEAR2D];

static int    INear2ddim,     Near2ddim[2];
static int    INear2dframe,   Near2dframe;
static int    INear2ddb,      Near2ddb;
static int    INear2dscale;
static double Near2dscale[2];
static int    INear2dcontour, Near2dcontour;
static int    INear2dobj,     Near2dobj;
static int    INear2dnoinc,   Near2dnoinc;
static int    INear2dzoom;
static double Near2dzoom[4];

/* === 計算 === */

void orcwa_init(void)
{
	memset(Title, 0, sizeof(Title));
	NXsection  = NYsection  = NZsection  = 0;
	NXdivision = NYdivision = NZdivision = 0;
	NMaterial  = 0;
	NGeometry  = 0;
	NFeed      = 0;
	IPlanewave = 0;
	NLoad      = 0;
	NPoint     = 0;
	IRfeed     = IPulsewidth = ITimestep = 0;
	IPml       = 0;
	IPbc       = 0;
	IFreq1     = IFreq2 = 0;
	ISolver    = 0;

	IMatchingloss = IPlotiter = IPlotfeed = IPlotpoint = IPlotsmith = 0;
	IPlotzin = IPlotyin = IPlotref = IPlotspara = IPlotcoupling = 0;
	IPlotfar0d = 0;
	IFreqdiv   = 0;
	NFar1d     = 0;
	IFar1dstyle = IFar1dcomponent = IFar1ddb = IFar1dnorm = IFar1dscale = 0;
	IPlotfar2d  = IFar2dcomponent = IFar2ddb = IFar2dscale = IFar2dobj  = 0;
	NNear1d    = 0;
	INear1ddb  = INear1dnoinc = INear1dscale = 0;
	NNear2d    = 0;
	INear2ddim = INear2dframe = INear2ddb = INear2dscale = 0;
	INear2dcontour = INear2dobj = INear2dnoinc = INear2dzoom = 0;
}

void orcwa_title(const char title[])
{
	strncpy(Title, title, MAXTITLE - 1);
}

static void setsection(const char fname[], int *num, double sec[], int n, va_list ap)
{
	if ((n < 2) || (n > MAXSECTION)) {
		printf("*** %sの引数が不正です\n", fname);
		return;
	}
	*num = n;
	for (int i = 0; i < n; i++) {
		sec[i] = va_arg(ap, double);
	}
}

static void setdivision(const char fname[], int *num, int div[], int n, va_list ap)
{
	if ((n < 1) || (n > MAXSECTION - 1)) {
		printf("*** %sの引数が不正です\n", fname);
		return;
	}
	*num = n;
	for (int i = 0; i < n; i++) {
		div[i] = va_arg(ap, int);
	}
}

void orcwa_xsection(int n, ...)
{
	va_list ap;
	va_start(ap, n);
	setsection("xsection", &NXsection, Xsection, n, ap);
	va_end(ap);
}

void orcwa_ysection(int n, ...)
{
	va_list ap;
	va_start(ap, n);
	setsection("ysection", &NYsection, Ysection, n, ap);
	va_end(ap);
}

void orcwa_zsection(int n, ...)
{
	va_list ap;
	va_start(ap, n);
	setsection("zsection", &NZsection, Zsection, n, ap);
	va_end(ap);
}

void orcwa_xdivision(int n, ...)
{
	va_list ap;
	va_start(ap, n);
	setdivision("xdivision", &NXdivision, Xdivision, n, ap);
	va_end(ap);
}

void orcwa_ydivision(int n, ...)
{
	va_list ap;
	va_start(ap, n);
	setdivision("ydivision", &NYdivision, Ydivision, n, ap);
	va_end(ap);
}

void orcwa_zdivision(int n, ...)
{
	va_list ap;
	va_start(ap, n);
	setdivision("zdivision", &NZdivision, Zdivision, n, ap);
	va_end(ap);
}

void orcwa_material(double epsr, double esgm, double amur, double msgm)
{
	if (NMaterial >= MAXMATERIAL) return;
	material_t *m = &MaterialData[NMaterial++];
	m->type = 1;
	m->d[0] = epsr;
	m->d[1] = esgm;
	m->d[2] = amur;
	m->d[3] = msgm;
}

void orcwa_material_dispersion(double einf, double ae, double be, double ce)
{
	if (NMaterial >= MAXMATERIAL) return;
	material_t *m = &MaterialData[NMaterial++];
	m->type = 2;
	m->d[0] = einf;
	m->d[1] = ae;
	m->d[2] = be;
	m->d[3] = ce;
}

void orcwa_geometry(int material, int shape, ...)
{
	if (NGeometry >= MAXGEOMETRY) return;
	const int s = shape;
	const int ng = ((s == 31) || (s == 32) || (s == 33)
	             || (s == 41) || (s == 42) || (s == 43)
	             || (s == 51) || (s == 52) || (s == 53)) ? 8 : 6;
	geometry_t *g = &GeometryData[NGeometry++];
	g->material = material;
	g->shape    = shape;
	g->ng       = ng;
	va_list ap;
	va_start(ap, shape);
	for (int i = 0; i < ng; i++) {
		g->g[i] = va_arg(ap, double);
	}
	va_end(ap);
}

void orcwa_feed(char dir, double x, double y, double z, double amp, double delay, double z0)
{
	if (NFeed >= MAXFEED) return;
	feed_t *f = &FeedData[NFeed++];
	f->dir   = dir;
	f->x     = x;
	f->y     = y;
	f->z     = z;
	f->amp   = amp;
	f->delay = delay;
	f->z0    = z0;
}

void orcwa_planewave(double theta, double phi, int pol)
{
	IPlanewave     = 1;
	PlanewaveTheta = theta;
	PlanewavePhi   = phi;
	PlanewavePol   = pol;
}

void orcwa_load(char dir, double x, double y, double z, char type, double rcl)
{
	if (NLoad >= MAXLOAD) return;
	load_t *l = &LoadData[NLoad++];
	l->dir  = dir;
	l->x    = x;
	l->y    = y;
	l->z    = z;
	l->type = type;
	l->rcl  = rcl;
}

void orcwa_point(char dir, double x, double y, double z, const char prop[])
{
	if (NPoint >= MAXPOINT) return;
	point_t *p = &PointData[NPoint++];
	p->dir = dir;
	p->x   = x;
	p->y   = y;
	p->z   = z;
	memset(p->prop, 0, sizeof(p->prop));
	if (prop != NULL) {
		strncpy(p->prop, prop, MAXPROP - 1);
	}
}

void orcwa_rfeed(double rfeed)
{
	IRfeed     = 1;
	RfeedValue = rfeed;
}

void orcwa_pulsewidth(double pulsewidth)
{
	IPulsewidth     = 1;
	PulsewidthValue = pulsewidth;
}

void orcwa_timestep(double timestep)
{
	ITimestep     = 1;
	TimestepValue = timestep;
}

void orcwa_pml(int l, double m, double r0)
{
	IPml  = 1;
	PmlL  = l;
	PmlM  = m;
	PmlR0 = r0;
}

void orcwa_pbc(int pbcx, int pbcy, int pbcz)
{
	IPbc = 1;
	PbcX = pbcx;
	PbcY = pbcy;
	PbcZ = pbcz;
}

void orcwa_frequency1(double fstart, double fend, int div)
{
	IFreq1   = 1;
	Freq1[0] = fstart;
	Freq1[1] = fend;
	Freq1Div = div;
}

void orcwa_frequency2(double fstart, double fend, int div)
{
	IFreq2   = 1;
	Freq2[0] = fstart;
	Freq2[1] = fend;
	Freq2Div = div;
}

void orcwa_solver(int maxiter, int nout, double converg)
{
	ISolver       = 1;
	SolverMaxiter = maxiter;
	SolverNout    = nout;
	SolverConverg = converg;
}

/* === ポスト処理 === */

void orcwa_matchingloss(int i0)
{
	IMatchingloss = 1;
	Matchingloss  = i0;
}

void orcwa_plotiter(int i0)
{
	IPlotiter = 1;
	Plotiter  = i0;
}

void orcwa_plotfeed(int i0)
{
	IPlotfeed = 1;
	Plotfeed  = i0;
}

void orcwa_plotpoint(int i0)
{
	IPlotpoint = 1;
	Plotpoint  = i0;
}

void orcwa_plotsmith(int i0)
{
	IPlotsmith = 1;
	Plotsmith  = i0;
}

static void setplotscale(int *flag, plotscale_t *p, int scale, double smin, double smax, int sdiv)
{
	*flag   = 1;
	p->scale = scale;
	p->smin  = smin;
	p->smax  = smax;
	p->sdiv  = sdiv;
}

void orcwa_plotzin(int scale, double smin, double smax, int sdiv)
{
	setplotscale(&IPlotzin, &Plotzin, scale, smin, smax, sdiv);
}

void orcwa_plotyin(int scale, double smin, double smax, int sdiv)
{
	setplotscale(&IPlotyin, &Plotyin, scale, smin, smax, sdiv);
}

void orcwa_plotref(int scale, double smin, double smax, int sdiv)
{
	setplotscale(&IPlotref, &Plotref, scale, smin, smax, sdiv);
}

void orcwa_plotspara(int scale, double smin, double smax, int sdiv)
{
	setplotscale(&IPlotspara, &Plotspara, scale, smin, smax, sdiv);
}

void orcwa_plotcoupling(int scale, double smin, double smax, int sdiv)
{
	setplotscale(&IPlotcoupling, &Plotcoupling, scale, smin, smax, sdiv);
}

void orcwa_plotfar0d(double theta, double phi, int scale, double smin, double smax, int sdiv)
{
	Far0dTheta = theta;
	Far0dPhi   = phi;
	setplotscale(&IPlotfar0d, &Plotfar0d, scale, smin, smax, sdiv);
}

void orcwa_freqdiv(int freqdiv)
{
	IFreqdiv = 1;
	Freqdiv  = freqdiv;
}

/* 遠方界1D */

void orcwa_plotfar1d(char dir, int div, double angle)
{
	if (NFar1d >= MAXFAR1D) return;
	far1d_t *f = &Far1dData[NFar1d++];
	f->dir   = dir;
	f->div   = div;
	f->angle = angle;
}

void orcwa_far1dstyle(int i0)
{
	IFar1dstyle = 1;
	Far1dstyle  = i0;
}

void orcwa_far1dcomponent(int i0, int i1, int i2)
{
	IFar1dcomponent   = 1;
	Far1dcomponent[0] = i0;
	Far1dcomponent[1] = i1;
	Far1dcomponent[2] = i2;
}

void orcwa_far1ddb(int i0)
{
	IFar1ddb = 1;
	Far1ddb  = i0;
}

void orcwa_far1dnorm(int i0)
{
	IFar1dnorm = 1;
	Far1dnorm  = i0;
}

void orcwa_far1dscale(double smin, double smax, int sdiv)
{
	IFar1dscale   = 1;
	Far1dscale[0] = smin;
	Far1dscale[1] = smax;
	Far1dscaleDiv = sdiv;
}

/* 遠方界2D */

void orcwa_plotfar2d(int divtheta, int divphi)
{
	IPlotfar2d    = 1;
	Far2dDivtheta = divtheta;
	Far2dDivphi   = divphi;
}

void orcwa_far2dcomponent(int i0, int i1, int i2, int i3, int i4, int i5, int i6)
{
	IFar2dcomponent   = 1;
	Far2dcomponent[0] = i0;
	Far2dcomponent[1] = i1;
	Far2dcomponent[2] = i2;
	Far2dcomponent[3] = i3;
	Far2dcomponent[4] = i4;
	Far2dcomponent[5] = i5;
	Far2dcomponent[6] = i6;
}

void orcwa_far2ddb(int i0)
{
	IFar2ddb = 1;
	Far2ddb  = i0;
}

void orcwa_far2dscale(double smin, double smax, int sdiv)
{
	IFar2dscale   = 1;
	Far2dscale[0] = smin;
	Far2dscale[1] = smax;
	Far2dscaleDiv = sdiv;
}

void orcwa_far2dobj(double obj)
{
	IFar2dobj = 1;
	Far2dobj  = obj;
}

/* 近傍界1D */

void orcwa_plotnear1d(const char component[], char dir, double p1, double p2)
{
	if (NNear1d >= MAXNEAR1D) return;
	near1d_t *n = &Near1dData[NNear1d++];
	memset(n->component, 0, sizeof(n->component));
	strncpy(n->component, component, MAXPROP - 1);
	n->dir = dir;
	n->p1  = p1;
	n->p2  = p2;
}

void orcwa_near1ddb(int i0)
{
	INear1ddb = 1;
	Near1ddb  = i0;
}

void orcwa_near1dnoinc(int i0)
{
	INear1dnoinc = 1;
	Near1dnoinc  = i0;
}

void orcwa_near1dscale(double smin, double smax, int sdiv)
{
	INear1dscale   = 1;
	Near1dscale[0] = smin;
	Near1dscale[1] = smax;
	Near1dscaleDiv = sdiv;
}

/* 近傍界2D */

void orcwa_plotnear2d(const char component[], char dir, double p)
{
	if (NNear2d >= MAXNEAR2D) return;
	near2d_t *n = &Near2dData[NNear2d++];
	memset(n->component, 0, sizeof(n->component));
	strncpy(n->component, component, MAXPROP - 1);
	n->dir = dir;
	n->p   = p;
}

void orcwa_near2ddim(int i0, int i1)
{
	INear2ddim   = 1;
	Near2ddim[0] = i0;
	Near2ddim[1] = i1;
}

void orcwa_near2dframe(int i0)
{
	INear2dframe = 1;
	Near2dframe  = i0;
}

void orcwa_near2ddb(int i0)
{
	INear2ddb = 1;
	Near2ddb  = i0;
}

void orcwa_near2dscale(double smin, double smax)
{
	INear2dscale   = 1;
	Near2dscale[0] = smin;
	Near2dscale[1] = smax;
}

void orcwa_near2dcontour(int i0)
{
	INear2dcontour = 1;
	Near2dcontour  = i0;
}

void orcwa_near2dobj(int i0)
{
	INear2dobj = 1;
	Near2dobj  = i0;
}

void orcwa_near2dnoinc(int i0)
{
	INear2dnoinc = 1;
	Near2dnoinc  = i0;
}

void orcwa_near2dzoom(double p0, double p1, double p2, double p3)
{
	INear2dzoom   = 1;
	Near2dzoom[0] = p0;
	Near2dzoom[1] = p1;
	Near2dzoom[2] = p2;
	Near2dzoom[3] = p3;
}

/* ファイル出力 */

static void outmesh(FILE *fp, const char key[], int nsec, const double sec[], int ndiv, const int div[])
{
	if (nsec < 2) return;
	if (ndiv != nsec - 1) {
		printf("*** %s : section(%d)とdivision(%d)の個数が不整合です\n", key, nsec, ndiv);
		return;
	}
	fprintf(fp, "%s =", key);
	for (int i = 0; i < nsec - 1; i++) {
		fprintf(fp, " %g %d", sec[i], div[i]);
	}
	fprintf(fp, " %g\n", sec[nsec - 1]);
}

static void outplotscale(FILE *fp, const char key[], int flag, const plotscale_t *p)
{
	if (!flag) return;
	fprintf(fp, "%s = %d", key, p->scale);
	if (p->scale == 2) {
		fprintf(fp, " %g %g %d", p->smin, p->smax, p->sdiv);
	}
	fprintf(fp, "\n");
}

void orcwa_outdata(const char filename[])
{
	FILE *fp = fopen(filename, "wt");
	if (fp == NULL) {
		printf("*** ファイルを開けません : %s\n", filename);
		return;
	}

	/* === 計算 === */

	/* ヘッダー */
	fprintf(fp, "%s\n", VERSION);

	/* タイトル */
	if (Title[0] != '\0') {
		fprintf(fp, "title = %s\n", Title);
	}

	/* メッシュ */
	outmesh(fp, "xmesh", NXsection, Xsection, NXdivision, Xdivision);
	outmesh(fp, "ymesh", NYsection, Ysection, NYdivision, Ydivision);
	outmesh(fp, "zmesh", NZsection, Zsection, NZdivision, Zdivision);

	/* 物性値 */
	for (int n = 0; n < NMaterial; n++) {
		const material_t *m = &MaterialData[n];
		fprintf(fp, "material = %d %g %g %g %g\n", m->type, m->d[0], m->d[1], m->d[2], m->d[3]);
	}

	/* 物体形状 */
	for (int n = 0; n < NGeometry; n++) {
		const geometry_t *g = &GeometryData[n];
		fprintf(fp, "geometry = %d %d", g->material, g->shape);
		for (int i = 0; i < g->ng; i++) {
			fprintf(fp, " %g", g->g[i]);
		}
		fprintf(fp, "\n");
	}

	/* 給電点 */
	for (int n = 0; n < NFeed; n++) {
		const feed_t *f = &FeedData[n];
		fprintf(fp, "feed = %c %g %g %g %g %g %g\n", f->dir, f->x, f->y, f->z, f->amp, f->delay, f->z0);
	}

	/* 平面波入射 */
	if (IPlanewave) {
		fprintf(fp, "planewave = %g %g %d\n", PlanewaveTheta, PlanewavePhi, PlanewavePol);
	}

	/* 負荷 */
	for (int n = 0; n < NLoad; n++) {
		const load_t *l = &LoadData[n];
		fprintf(fp, "load = %c %g %g %g %c %g\n", l->dir, l->x, l->y, l->z, l->type, l->rcl);
	}

	/* 観測点 */
	for (int n = 0; n < NPoint; n++) {
		const point_t *p = &PointData[n];
		fprintf(fp, "point = %c %g %g %g %s\n", p->dir, p->x, p->y, p->z, p->prop);
	}

	/* その他 */
	if (IRfeed) {
		fprintf(fp, "rfeed = %g\n", RfeedValue);
	}
	if (IPulsewidth) {
		fprintf(fp, "pulsewidth = %g\n", PulsewidthValue);
	}
	if (ITimestep) {
		fprintf(fp, "timestep = %g\n", TimestepValue);
	}

	/* PML/PBC */
	if (IPml) {
		fprintf(fp, "abc = 1 %d %g %g\n", PmlL, PmlM, PmlR0);
	}
	if (IPbc) {
		fprintf(fp, "pbc = %d %d %d\n", PbcX, PbcY, PbcZ);
	}

	/* 周波数 */
	if (IFreq1) {
		fprintf(fp, "frequency1 = %g %g %d\n", Freq1[0], Freq1[1], Freq1Div);
	}
	if (IFreq2) {
		fprintf(fp, "frequency2 = %g %g %d\n", Freq2[0], Freq2[1], Freq2Div);
	}

	/* 計算条件 */
	if (ISolver) {
		fprintf(fp, "solver = %d %d %g\n", SolverMaxiter, SolverNout, SolverConverg);
	}

	/* === ポスト処理 === */

	/* 周波数特性 */
	if (IMatchingloss) {
		fprintf(fp, "matchingloss = %d\n", Matchingloss);
	}
	if (IPlotiter) {
		fprintf(fp, "plotiter = %d\n", Plotiter);
	}
	if (IPlotfeed) {
		fprintf(fp, "plotfeed = %d\n", Plotfeed);
	}
	if (IPlotpoint) {
		fprintf(fp, "plotpoint = %d\n", Plotpoint);
	}
	if (IPlotsmith) {
		fprintf(fp, "plotsmith = %d\n", Plotsmith);
	}
	outplotscale(fp, "plotzin",      IPlotzin,      &Plotzin);
	outplotscale(fp, "plotyin",      IPlotyin,      &Plotyin);
	outplotscale(fp, "plotref",      IPlotref,      &Plotref);
	outplotscale(fp, "plotspara",    IPlotspara,    &Plotspara);
	outplotscale(fp, "plotcoupling", IPlotcoupling, &Plotcoupling);
	if (IPlotfar0d) {
		fprintf(fp, "plotfar0d = %g %g %d", Far0dTheta, Far0dPhi, Plotfar0d.scale);
		if (Plotfar0d.scale == 2) {
			fprintf(fp, " %g %g %d", Plotfar0d.smin, Plotfar0d.smax, Plotfar0d.sdiv);
		}
		fprintf(fp, "\n");
	}
	if (IFreqdiv) {
		fprintf(fp, "freqdiv = %d\n", Freqdiv);
	}

	/* 遠方界1D */
	for (int n = 0; n < NFar1d; n++) {
		const far1d_t *f = &Far1dData[n];
		fprintf(fp, "plotfar1d = %c %d", f->dir, f->div);
		if (f->angle != 0) {
			fprintf(fp, " %g", f->angle);
		}
		fprintf(fp, "\n");
	}
	if (IFar1dstyle) {
		fprintf(fp, "far1dstyle = %d\n", Far1dstyle);
	}
	if (IFar1dcomponent) {
		fprintf(fp, "far1dcomponent = %d %d %d\n", Far1dcomponent[0], Far1dcomponent[1], Far1dcomponent[2]);
	}
	if (IFar1ddb) {
		fprintf(fp, "far1ddb = %d\n", Far1ddb);
	}
	if (IFar1dnorm) {
		fprintf(fp, "far1dnorm = %d\n", Far1dnorm);
	}
	if (IFar1dscale) {
		fprintf(fp, "far1dscale = %g %g %d\n", Far1dscale[0], Far1dscale[1], Far1dscaleDiv);
	}

	/* 遠方界2D */
	if (IPlotfar2d) {
		fprintf(fp, "plotfar2d = %d %d\n", Far2dDivtheta, Far2dDivphi);
	}
	if (IFar2dcomponent) {
		fprintf(fp, "far2dcomponent = %d %d %d %d %d %d %d\n",
			Far2dcomponent[0], Far2dcomponent[1], Far2dcomponent[2], Far2dcomponent[3],
			Far2dcomponent[4], Far2dcomponent[5], Far2dcomponent[6]);
	}
	if (IFar2ddb) {
		fprintf(fp, "far2ddb = %d\n", Far2ddb);
	}
	if (IFar2dscale) {
		fprintf(fp, "far2dscale = %g %g %d\n", Far2dscale[0], Far2dscale[1], Far2dscaleDiv);
	}
	if (IFar2dobj) {
		fprintf(fp, "far2dobj = %g\n", Far2dobj);
	}

	/* 近傍界1D */
	for (int n = 0; n < NNear1d; n++) {
		const near1d_t *p = &Near1dData[n];
		fprintf(fp, "plotnear1d = %s %c %g %g\n", p->component, p->dir, p->p1, p->p2);
	}
	if (INear1ddb) {
		fprintf(fp, "near1ddb = %d\n", Near1ddb);
	}
	if (INear1dnoinc) {
		fprintf(fp, "near1dnoinc = %d\n", Near1dnoinc);
	}
	if (INear1dscale) {
		fprintf(fp, "near1dscale = %g %g %d\n", Near1dscale[0], Near1dscale[1], Near1dscaleDiv);
	}

	/* 近傍界2D */
	for (int n = 0; n < NNear2d; n++) {
		const near2d_t *p = &Near2dData[n];
		fprintf(fp, "plotnear2d = %s %c %g\n", p->component, p->dir, p->p);
	}
	if (INear2ddim) {
		fprintf(fp, "near2ddim = %d %d\n", Near2ddim[0], Near2ddim[1]);
	}
	if (INear2dframe) {
		fprintf(fp, "near2dframe = %d\n", Near2dframe);
	}
	if (INear2ddb) {
		fprintf(fp, "near2ddb = %d\n", Near2ddb);
	}
	if (INear2dscale) {
		fprintf(fp, "near2dscale = %g %g\n", Near2dscale[0], Near2dscale[1]);
	}
	if (INear2dcontour) {
		fprintf(fp, "near2dcontour = %d\n", Near2dcontour);
	}
	if (INear2dobj) {
		fprintf(fp, "near2dobj = %d\n", Near2dobj);
	}
	if (INear2dnoinc) {
		fprintf(fp, "near2dnoinc = %d\n", Near2dnoinc);
	}
	if (INear2dzoom) {
		fprintf(fp, "near2dzoom = %g %g %g %g\n", Near2dzoom[0], Near2dzoom[1], Near2dzoom[2], Near2dzoom[3]);
	}

	fprintf(fp, "end\n");

	fclose(fp);

	printf("output -> %s\n", filename);
}
