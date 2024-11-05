/*
input2.c
*/

#include "orcwa.h"
#include "orcwa_prototype.h"


// index and length
static void getindex(char dir, double x, double y, double z,
    int *i, int *j, int *k, double *dx, double *dy, double *dz)
{
    if      (dir == 'X') {
        *i = nearest(x, 0, Nx - 1, Xc);
        *j = nearest(y, 0, Ny,     Yn);
        *k = nearest(z, 0, Nz,     Zn);
        *dx = Xn[*i + 1] - Xn[*i];
        *dy = (Ny > 1) ? (Yc[*j] - Yc[*j - 1]) : (Yn[1] - Yn[0]);
        *dz = (Nz > 1) ? (Zc[*k] - Zc[*k - 1]) : (Zn[1] - Zn[0]);
    }
    else if (dir == 'Y') {
        *j = nearest(y, 0, Ny - 1, Yc);
        *k = nearest(z, 0, Nz,     Zn);
        *i = nearest(x, 0, Nx,     Xn);
        *dy = Yn[*j + 1] - Yn[*j];
        *dz = (Nz > 1) ? (Zc[*k] - Zc[*k - 1]) : (Zn[1] - Zn[0]);
        *dx = (Nx > 1) ? (Xc[*i] - Xc[*i - 1]) : (Xn[1] - Xn[0]);
    }
    else if (dir == 'Z') {
        *k = nearest(z, 0, Nz - 1, Zc);
        *i = nearest(x, 0, Nx,     Xn);
        *j = nearest(y, 0, Ny,     Yn);
        *dz = Zn[*k + 1] - Zn[*k];
        *dx = (Nx > 1) ? (Xc[*i] - Xc[*i - 1]) : (Xn[1] - Xn[0]);
        *dy = (Ny > 1) ? (Yc[*j] - Yc[*j - 1]) : (Yn[1] - Yn[0]);
    }
}


// number of cells
void setup_cells(int nxr, int nyr, int nzr, int *dxr, int *dyr, int *dzr)
{
    int xsum = 0;
    for (int i = 0; i < nxr; i++) {
        xsum += dxr[i];
    }
    Nx = xsum;

    int ysum = 0;
    for (int j = 0; j < nyr; j++) {
        ysum += dyr[j];
    }
    Ny = ysum;

    int zsum = 0;
    for (int k = 0; k < nzr; k++) {
        zsum += dzr[k];
    }
    Nz = zsum;
}


// node
void setup_node(int nxr, int nyr, int nzr, double *xr, double *yr, double *zr, int *dxr, int *dyr, int *dzr)
{
    if ((nxr <= 0) || (nyr <= 0) || (nzr <= 0)) return;

    Xn = (double *)malloc((Nx + 1) * sizeof(double));
    Yn = (double *)malloc((Ny + 1) * sizeof(double));
    Zn = (double *)malloc((Nz + 1) * sizeof(double));

    int xid = 0;
    for (int mx = 0; mx < nxr; mx++) {
        double dx = (xr[mx + 1] - xr[mx]) / dxr[mx];
        for (int i = 0; i < dxr[mx]; i++) {
            Xn[xid++] = xr[mx] + (i * dx);
        }
    }
    Xn[xid] = xr[nxr];

    int yid = 0;
    for (int my = 0; my < nyr; my++) {
        double dy = (yr[my + 1] - yr[my]) / dyr[my];
        for (int j = 0; j < dyr[my]; j++) {
            Yn[yid++] = yr[my] + (j * dy);
        }
    }
    Yn[yid] = yr[nyr];

    int zid = 0;
    for (int mz = 0; mz < nzr; mz++) {
        double dz = (zr[mz + 1] - zr[mz]) / dzr[mz];
        for (int k = 0; k < dzr[mz]; k++) {
            Zn[zid++] = zr[mz] + (k * dz);
        }
    }
    Zn[zid] = zr[nzr];
/*
    // debug
    for (int i = 0; i <= Nx; i++) {
        printf("Xn[%d]=%.5f\n", i, Xn[i] * 1e3);
    }
    for (int j = 0; j <= Ny; j++) {
        printf("Yn[%d]=%.5f\n", j, Yn[j] * 1e3);
    }
    for (int k = 0; k <= Nz; k++) {
        printf("Zn[%d]=%.5f\n", k, Zn[k] * 1e3);
    }
*/
}


// cell center
void setup_center(void)
{
    Xc = (double *)malloc(Nx * sizeof(double));
    Yc = (double *)malloc(Ny * sizeof(double));
    Zc = (double *)malloc(Nz * sizeof(double));

    for (int i = 0; i < Nx; i++) {
        Xc[i] = (Xn[i] + Xn[i + 1]) / 2;
    }
    for (int j = 0; j < Ny; j++) {
        Yc[j] = (Yn[j] + Yn[j + 1]) / 2;
    }
    for (int k = 0; k < Nz; k++) {
        Zc[k] = (Zn[k] + Zn[k + 1]) / 2;
    }
/*
    // debug
    for (int i = 0; i < Nx; i++) {
        printf("Xc[%d]=%.5f\n", i, Xc[i] * 1e3);
    }
    for (int j = 0; j < Ny; j++) {
        printf("Yc[%d]=%.5f\n", j, Yc[j] * 1e3);
    }
    for (int k = 0; k < Nz; k++) {
        printf("Zc[%d]=%.5f\n", k, Zc[k] * 1e3);
    }
*/
}


// fit geometry without thickness to the nearest node
void fitgeometry(void)
{
    if ((Nx <= 0) || (Ny <= 0) || (Nz <= 0)) return;

    const double d0 = EPS * (
        fabs(Xn[Nx] - Xn[0]) +
        fabs(Yn[Ny] - Yn[0]) +
        fabs(Zn[Nz] - Zn[0]));

    for (int n = 0; n < NGeometry; n++) {
        int    shape = Geometry[n].shape;
        double *g = Geometry[n].g;

        if      ((shape == 1) || (shape == 2)) {
            if (fabs(g[0] - g[1]) < d0) {
                int i = nearest(g[0], 0, Nx, Xn);
                g[0] = g[1] = Xn[i];
            }
            if (fabs(g[2] - g[3]) < d0) {
                int j = nearest(g[2], 0, Ny, Yn);
                g[2] = g[3] = Yn[j];
            }
            if (fabs(g[4] - g[5]) < d0) {
                int k = nearest(g[4], 0, Nz, Zn);
                g[4] = g[5] = Zn[k];
            }
        }
        //else if ((shape == 11) || (shape == 64) || (shape == 65)) {
        else if (shape == 11) {
            if (fabs(g[0] - g[1]) < d0) {
                int i = nearest(g[0], 0, Nx, Xn);
                g[0] = g[1] = Xn[i];
            }
        }
        //else if ((shape == 12) || (shape == 61) || (shape == 66)) {
        else if (shape == 12) {
            if (fabs(g[2] - g[3]) < d0) {
                int j = nearest(g[2], 0, Ny, Yn);
                g[2] = g[3] = Yn[j];
            }
        }
        //else if ((shape == 13) || (shape == 62) || (shape == 63)) {
        else if (shape == 13) {
            if (fabs(g[4] - g[5]) < d0) {
                int k = nearest(g[4], 0, Nz, Zn);
                g[4] = g[5] = Zn[k];
            }
        }
    }
}

double waveform(double t)
{
    const double f0 = (Freq2[0] + Freq2[NFreq2 - 1]) / 2;

    // ガウスパルスの例
    double tau = 1e-15; // パルス幅
    double t0 = 5 * tau; // パルスの中心時間
    return exp(-pow((t - t0) / tau, 2)) * cos(2 * PI * f0 * (t - t0));
}

void setup_planewave2(void)
{
    // 光の波長と周波数の設定
    //const double wavelength = 500e-9; // 500 nm
    //const double f0 = C / wavelength; // 周波数
    const double f0 = (Freq2[0] + Freq2[NFreq2 - 1]) / 2;

    // 入射角の設定（必要に応じて変更）
    //Planewave.theta = 0.0; // 入射角 θ
    //Planewave.phi = 0.0;   // 入射角 φ

    // ラジアンへの変換
    const double cost = cos(Planewave.theta * DTOR);
    const double sint = sin(Planewave.theta * DTOR);
    const double cosp = cos(Planewave.phi   * DTOR);
    const double sinp = sin(Planewave.phi   * DTOR);

    // 単位ベクトルの計算
    double r1[3], t1[3], p1[3];
    r1[0] = + sint * cosp;
    r1[1] = + sint * sinp;
    r1[2] = + cost;
    t1[0] = + cost * cosp;
    t1[1] = + cost * sinp;
    t1[2] = - sint;
    p1[0] = - sinp;
    p1[1] = + cosp;
    p1[2] = 0;

    // 伝搬ベクトルの設定
    Planewave.ri[0] = - r1[0];
    Planewave.ri[1] = - r1[1];
    Planewave.ri[2] = - r1[2];

    // 電界の設定
    if (Planewave.pol == 1) {
        // 垂直偏波
        Planewave.ei[0] = - t1[0];
        Planewave.ei[1] = - t1[1];
        Planewave.ei[2] = - t1[2];
    }
    else if (Planewave.pol == 2) {
        // 水平偏波
        Planewave.ei[0] = + p1[0];
        Planewave.ei[1] = + p1[1];
        Planewave.ei[2] = + p1[2];
    }

    // 磁界の計算（H = E × r）
    Planewave.hi[0] = (Planewave.ei[1] * r1[2]) - (Planewave.ei[2] * r1[1]);
    Planewave.hi[1] = (Planewave.ei[2] * r1[0]) - (Planewave.ei[0] * r1[2]);
    Planewave.hi[2] = (Planewave.ei[0] * r1[1]) - (Planewave.ei[1] * r1[0]);

    // 初期位置の設定
    const double r = sqrt((Xn[0] - Xn[Nx]) * (Xn[0] - Xn[Nx]) +
                          (Yn[0] - Yn[Ny]) * (Yn[0] - Yn[Ny]) +
                          (Zn[0] - Zn[Nz]) * (Zn[0] - Zn[Nz])) / 2 + (0.5 * C / f0);
    printf("%f\n", r);
    Planewave.r0[0] = ((Xn[0] + Xn[Nx]) / 2 - (r * Planewave.ri[0]));
    Planewave.r0[1] = ((Yn[0] + Yn[Ny]) / 2 - (r * Planewave.ri[1]));
    Planewave.r0[2] = ((Zn[0] + Zn[Nz]) / 2 - (r * Planewave.ri[2]));
    printf("%f %f %f\n", Planewave.r0[0], Planewave.r0[1], Planewave.r0[2]);

    // 波形パラメータの設定
    Planewave.ai = 1.0; // 振幅を適切に設定

    // 周波数配列の設定（単一周波数の場合）
    //Freq2[0] = f0;
    //NFreq2 = 1;
}

// plane wave
void setup_planewave(void)
{
    printf("setup_planewave.\n");

    const double cost = cos(Planewave.theta * DTOR);
    const double sint = sin(Planewave.theta * DTOR);
    const double cosp = cos(Planewave.phi   * DTOR);
    const double sinp = sin(Planewave.phi   * DTOR);

    // unit vector in (r, theta, phi)
    double r1[3], t1[3], p1[3];
    r1[0] = + sint * cosp;
    r1[1] = + sint * sinp;
    r1[2] = + cost;
    t1[0] = + cost * cosp;
    t1[1] = + cost * sinp;
    t1[2] = - sint;
    p1[0] = - sinp;
    p1[1] = + cosp;
    p1[2] = 0;

    // propagation vector = - r1
    Planewave.ri[0] = - r1[0];
    Planewave.ri[1] = - r1[1];
    Planewave.ri[2] = - r1[2];
    printf("propagation vector.\n");
    printf("%.10f %.10f %.10f\n", Planewave.r0[0], Planewave.r0[1], Planewave.r0[2]);

    // E
    if      (Planewave.pol == 1) {
        // V-pol
        Planewave.ei[0] = - t1[0];
        Planewave.ei[1] = - t1[1];
        Planewave.ei[2] = - t1[2];
    }
    else if (Planewave.pol == 2) {
        // H-pol
        Planewave.ei[0] = + p1[0];
        Planewave.ei[1] = + p1[1];
        Planewave.ei[2] = + p1[2];
    }
    printf("E.\n");
    printf("%.10f %.10f %.10f\n", Planewave.ei[0], Planewave.ei[1], Planewave.ei[2]);

    // H = E X r
    Planewave.hi[0] = (Planewave.ei[1] * r1[2]) - (Planewave.ei[2] * r1[1]);
    Planewave.hi[1] = (Planewave.ei[2] * r1[0]) - (Planewave.ei[0] * r1[2]);
    Planewave.hi[2] = (Planewave.ei[0] * r1[1]) - (Planewave.ei[1] * r1[0]);
    printf("H.\n");
    printf("%.10f %.10f %.10f\n", Planewave.hi[0], Planewave.hi[1], Planewave.hi[2]);

    // initial position
    const double f0 = (Freq2[0] + Freq2[NFreq2 - 1]) / 2;
    const double r = sqrt((Xn[0] - Xn[Nx]) * (Xn[0] - Xn[Nx]) +
                          (Yn[0] - Yn[Ny]) * (Yn[0] - Yn[Ny]) +
                          (Zn[0] - Zn[Nz]) * (Zn[0] - Zn[Nz])) / 2 + (0.5 * C / f0);
    //const double r = 10;
    printf("%.10f\n", r); // 0.000003

    Planewave.r0[0] = ((Xn[0] + Xn[Nx]) / 2 - (r * Planewave.ri[0]));
    Planewave.r0[1] = ((Yn[0] + Yn[Ny]) / 2 - (r * Planewave.ri[1]));
    Planewave.r0[2] = ((Zn[0] + Zn[Nz]) / 2 - (r * Planewave.ri[2]));
    Planewave.r0[2] = 0.0000005;
    printf("initial position.\n");
    printf("%.10f %.10f %.10f\n", Planewave.r0[0], Planewave.r0[1], Planewave.r0[2]);

    // waveform parameter(ガウス微分パルス)
    Planewave.ai = 4 / (1.27 / f0);
}


// feed
void setup_feed(const double *x, const double *y, const double *z)
{
    for (int n = 0; n < NFeed; n++) {
        getindex(Feed[n].dir, x[n], y[n], z[n],
            &Feed[n].i, &Feed[n].j, &Feed[n].k, &Feed[n].dx, &Feed[n].dy, &Feed[n].dz);
    }
}


// point
void setup_point(const double *x, const double *y, const double *z, const char strprop[])
{
    for (int n = 0; n < NPoint; n++) {
        getindex(Point[n].dir, x[n], y[n], z[n],
            &Point[n].i, &Point[n].j, &Point[n].k, &Point[n].dx, &Point[n].dy, &Point[n].dz);
    }

    // add 1+, 1- points
    Point = (point_t *)realloc(Point, (NPoint + 2) * sizeof(point_t));
    Point[NPoint] = Point[NPoint + 1] = Point[0];
    if      (!strcmp(strprop, "+X") || !strcmp(strprop, "-X")) {
        Point[NPoint].i     = Point[0].i + (!strcmp(strprop, "+X") ? +1 : -1);
        Point[NPoint + 1].i = Point[0].i + (!strcmp(strprop, "+X") ? -1 : +1);
        Point[NPoint].j = Point[NPoint + 1].j = Point[0].j;
        Point[NPoint].k = Point[NPoint + 1].k = Point[0].k;
    }
    else if (!strcmp(strprop, "+Y") || !strcmp(strprop, "-Y")) {
        Point[NPoint].j     = Point[0].j + (!strcmp(strprop, "+Y") ? +1 : -1);
        Point[NPoint + 1].j = Point[0].j + (!strcmp(strprop, "+Y") ? -1 : +1);
        Point[NPoint].k = Point[NPoint + 1].k = Point[0].k;
        Point[NPoint].i = Point[NPoint + 1].i = Point[0].i;
    }
    else if (!strcmp(strprop, "+Z") || !strcmp(strprop, "-Z")) {
        Point[NPoint].k     = Point[0].k + (!strcmp(strprop, "+Z") ? +1 : -1);
        Point[NPoint + 1].k = Point[0].k + (!strcmp(strprop, "+Z") ? -1 : +1);
        Point[NPoint].i = Point[NPoint + 1].i = Point[0].i;
        Point[NPoint].j = Point[NPoint + 1].j = Point[0].j;
    }
/*
    for (int n = 0; n < NPoint + 2; n++) {
        printf("%d %c %d %d %d %f %f %f\n", n, Point[n].dir, Point[n].i, Point[n].j, Point[n].k, Point[n].dx, Point[n].dy, Point[n].dz);
    }
*/
}


// load
void setup_load(int nload, char *dload, double *xload, double *yload, double *zload, char *cload, double *pload, int array_inc)
{
    for (int n = 0; n < nload; n++) {
        const char dir = dload[n];
        const char rcl = cload[n];
        int    i = 0, j = 0, k = 0;
        double dx = 0, dy = 0, dz = 0;
        getindex(dload[n], xload[n], yload[n], zload[n], &i, &j, &k, &dx, &dy, &dz);
        const double dlds = (dir == 'X') ? (dx / (dy * dz))
                          : (dir == 'Y') ? (dy / (dz * dx))
                          : (dir == 'Z') ? (dz / (dx * dy)) : 0;

        if      ((rcl == 'R') || (rcl == 'C')) {
            // material
            if (NMaterial % array_inc == 2) {
                Material = (material_t *)realloc(Material, (NMaterial + array_inc) * sizeof(material_t));
            }
            if (NMaterial >= MAXMATERIAL) {
                fprintf(stderr, "*** too many load(R,C)\n");
                exit(1);
            }
            Material[NMaterial].type = 1;
            Material[NMaterial].epsr = (rcl == 'R') ? 1 : (pload[n] / EPS0) * dlds;
            Material[NMaterial].esgm = (rcl == 'R') ? (1 / pload[n]) * dlds : 0;
            Material[NMaterial].amur = 1;
            Material[NMaterial].msgm = 0;
            //printf("%zd %d %f %f %f %f\n", NMaterial, Material[NMaterial].type, Material[NMaterial].epsr, Material[NMaterial].esgm, Material[NMaterial].amur, Material[NMaterial].msgm);
            // geometry
            if (NGeometry % array_inc == 0) {
                Geometry = (geometry_t *)realloc(Geometry, (NGeometry + array_inc) * sizeof(geometry_t));
            }
            geometry_t *ptr = &Geometry[NGeometry];
            ptr->m = (id_t)NMaterial;
            ptr->shape = 1;
            if      (dir == 'X') {
                ptr->g[0] = Xn[i    ];
                ptr->g[1] = Xn[i + 1];
                ptr->g[2] =
                ptr->g[3] = Yn[j];
                ptr->g[4] =
                ptr->g[5] = Zn[k];
            }
            else if (dir == 'Y') {
                ptr->g[2] = Yn[j    ];
                ptr->g[3] = Yn[j + 1];
                ptr->g[4] =
                ptr->g[5] = Zn[k];
                ptr->g[0] =
                ptr->g[1] = Xn[i];
            }
            else if (dir == 'Z') {
                ptr->g[4] = Zn[k    ];
                ptr->g[5] = Zn[k + 1];
                ptr->g[0] =
                ptr->g[1] = Xn[i];
                ptr->g[2] =
                ptr->g[3] = Yn[j];
            }
            NMaterial++;
            NGeometry++;
        }
        else if (rcl == 'L') {
            Inductor = (inductor_t *)realloc(Inductor, (NInductor + 1) * sizeof(inductor_t));
            Inductor[NInductor].dir = dir;
            Inductor[NInductor].i = i;
            Inductor[NInductor].j = j;
            Inductor[NInductor].k = k;
            Inductor[NInductor].dx = dx;
            Inductor[NInductor].dy = dy;
            Inductor[NInductor].dz = dz;
            Inductor[NInductor].fctr = MU0 * dlds / pload[n];
            //printf("%d %c %d %d %d %e\n", NInductor, Inductor[NInductor].dir, Inductor[NInductor].i, Inductor[NInductor].j, Inductor[NInductor].k, Inductor[NInductor].fctr);
            Inductor[NInductor].e =
            Inductor[NInductor].esum = 0;
            NInductor++;
        }
    }
}


// setup geometry lines
void setup_geomline3d(void)
{
    if (!NGeometry) return;

    // alloc
    int *shape     =         (int *)malloc(NGeometry * sizeof(int));
    int *mid       =         (int *)malloc(NGeometry * sizeof(int));
    double (*g)[8] = (double (*)[8])malloc(NGeometry * 8 * sizeof(double));
    for (int n = 0; n < NGeometry; n++) {
        shape[n] = Geometry[n].shape;
        mid[n] = Geometry[n].m;
        memcpy(g[n], Geometry[n].g, 8 * sizeof(double));
    }

    // number of lines
    NGline = geomlines(0, (int)NGeometry, shape, NULL, NULL, NULL, NULL, 0);

    // alloc
    Gline  = (double (*)[2][3])malloc(NGline * 2 * 3 * sizeof(double));
    MGline =           (id_t *)malloc(NGline * sizeof(id_t));
    int *mgline =       (int *)malloc(NGline * sizeof(int));
    memset(Gline, 0, NGline * 2 * 3 * sizeof(double));
    memset(mgline, 0, NGline * sizeof(int));

    // set geometry
    const double eps = EPS * (fabs(Xn[Nx] - Xn[0]) + fabs(Yn[Ny] - Yn[0]) + fabs(Zn[Nz] - Zn[0]));
    NGline = geomlines(1, (int)NGeometry, shape, mid, g, Gline, mgline, eps);

    // set index
    for (int n = 0; n < NGline; n++) {
        MGline[n] = (id_t)mgline[n];
    }

    // free
    free(shape);
    free(mid);
    free(g);
    free(mgline);
}

