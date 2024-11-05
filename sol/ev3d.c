/*
ev3d.c

output: ev3 or HTML(3D)

ev3 format: (ascii or binary)
-1                                            new page
2 I X1 Y1 Z1 X2 Y2 Z2 R G B                   draw a line
3 I X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 R G B          fill a triangle
4 I X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 X4 Y4 Z4 R G B fill a quadrangle
-4                                            draw a title (support multiple data)
characters                                    title

id, I   : int
X, Y, Z : float
R, G, B : unsigned char
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

typedef struct {
	int idx;  // 2:line, 3:triangle(fill), 4:quadrangle(fill), -4:title
	double xvt[4], yvt[4], zvt[4];
	unsigned char rgb[3];
	char *str;
	int num;
} ev3d_t;

int ev3d_type = 1;
int ev3d_binary = 0;
char ev3d_fn[BUFSIZ] = "ev.ev3";
int ev3d_page = -1;
int ev3d_num = 0;
int *ev3d_ndata;
ev3d_t **ev3d_data;
unsigned char ev3d_rgb[3] = {0, 0, 0};
int ev3d_html_width = 500;
int ev3d_html_height = 500;
int ev3d_html_fontsize = 13;
double ev3d_html_theta = 60;
double ev3d_html_phi = 30;
char ev3d_html_fontname[BUFSIZ] = "monospace";
const int ev3d_increment = 10000;


// alloc increment
static void alloc_inc(void)
{
	if ((ev3d_ndata[ev3d_page]) > 0 && (ev3d_ndata[ev3d_page] % ev3d_increment == 0)) {
		ev3d_data[ev3d_page] = (ev3d_t *)realloc(ev3d_data[ev3d_page], (ev3d_ndata[ev3d_page] + ev3d_increment) * sizeof(ev3d_t));
	}
}


// initialize
void ev3d_init(void)
{
	// page
	ev3d_page = -1;

	// alloc
	ev3d_ndata = (int *)    malloc(1 * sizeof(int));
	ev3d_data  = (ev3d_t **)malloc(1 * sizeof(ev3d_t *));
}


// new page
void ev3d_newPage(void)
{
	// skip no data page
	if ((ev3d_page >= 0) && !ev3d_ndata[ev3d_page]) return;

	// page++
	ev3d_page++;

	// alloc
	if (ev3d_page > 0) {
		ev3d_ndata = (int *)    realloc(ev3d_ndata, (ev3d_page + 1) * sizeof(int));
		ev3d_data  = (ev3d_t **)realloc(ev3d_data,  (ev3d_page + 1) * sizeof(ev3d_t *));
	}
	ev3d_ndata[ev3d_page] = 0;
	ev3d_data[ev3d_page] = (ev3d_t *)malloc(ev3d_increment * sizeof(ev3d_t));

	// clear
	ev3d_rgb[0] = ev3d_rgb[1] = ev3d_rgb[2] = 0;
	ev3d_num = 0;
}


// end (ev.ev3, ascii)
static void ev3d_end_ev3_ascii(FILE *fp)
{
	for (int page = 0; page <= ev3d_page; page++) {
		if (!ev3d_ndata[page]) continue;
		fprintf(fp, "-1\n");
		for (int n = 0; n < ev3d_ndata[page]; n++) {
			ev3d_t *ptr = &ev3d_data[page][n];
			fprintf(fp, "%d", ptr->idx);
			if (ptr->idx > 0) {
				fprintf(fp, " %d", ptr->num);
				for (int m = 0; m < ptr->idx; m++) {
					fprintf(fp, " %g %g %g", ptr->xvt[m], ptr->yvt[m], ptr->zvt[m]);
				}
				fprintf(fp, " %u %u %u\n", ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
			}
			else if (ptr->idx == -4) {
				// title
				fprintf(fp, "\n%s\n", ptr->str);
			}
		}
	}
}


// end (ev.ev3, binary)
static void ev3d_end_ev3_binary(FILE *fp)
{
	float pos[3];

	for (int page = 0; page <= ev3d_page; page++) {
		if (!ev3d_ndata[page]) continue;
		const int i1 = -1;
		fwrite(&i1, sizeof(int), 1, fp);
		for (int n = 0; n < ev3d_ndata[page]; n++) {
			ev3d_t *ptr = &ev3d_data[page][n];
			const int i2 = ptr->idx;
			const int i3 = ptr->num;
			fwrite(&i2, sizeof(int), 1, fp);
			fwrite(&i3, sizeof(int), 1, fp);
			if (ptr->idx > 0) {
				for (int m = 0; m < ptr->idx; m++) {
					pos[0] = (float)ptr->xvt[m];
					pos[1] = (float)ptr->yvt[m];
					pos[2] = (float)ptr->zvt[m];
					fwrite(pos, sizeof(float), 3, fp);
				}
				fwrite(ptr->rgb, sizeof(unsigned char), 3, fp);
			}
			/*
			else if (ptr->idx == -3) {
				// string
				pos[0] = (float)ptr->xvt[0];
				pos[1] = (float)ptr->yvt[0];
				pos[2] = (float)ptr->zvt[0];
				fwrite(pos, sizeof(float), 3, fp);
				float hgt = (float)ptr->hgt;
				fwrite(&hgt, sizeof(float), 1, fp);
				fwrite(ptr->rgb, sizeof(unsigned char), 3, fp);
				const int i4 = (int)strlen(ptr->str);
				fwrite(&i4, sizeof(int), 1, fp);
				//const int i5 = ptr->fontname;
				//fwrite(&i5, sizeof(int), 1, fp);
				fwrite(ptr->str, sizeof(char), i4, fp);
			}
			*/
			else if (ptr->idx == -4) {
				// title
				//float hgt = (float)ptr->hgt;
				//fwrite(&hgt, sizeof(float), 1, fp);
				//fwrite(ptr->rgb, sizeof(unsigned char), 3, fp);
				const int i4 = (int)strlen(ptr->str);
				fwrite(&i4, sizeof(int), 1, fp);
				//const int i5 = ptr->fontname;
				//fwrite(&i5, sizeof(int), 1, fp);
				fwrite(ptr->str, sizeof(char), i4, fp);
			}
		}
	}
}


// end (HTML)
static void ev3d_end_html(FILE *fp, const char fn[])
{
	//char fontname[][BUFSIZ] = {"sansserif", "serif", "monospace"};
	const char fmt[] = "[%g,%g,%g,%g,%g,%g,%d,%d,%d],\n";
	int *ndata = (int *)malloc((ev3d_page + 1) * sizeof(int));

	// header
	fprintf(fp, "<!doctype html>\n");
	fprintf(fp, "<html>\n");
	fprintf(fp, "<head>\n");
	fprintf(fp, "<meta charset=\"UTF-8\">\n");
	fprintf(fp, "<title>%s</title>\n", fn);
	fprintf(fp, "</head>\n");

	// <body>
	fprintf(fp, "<body onload=\"init()\">\n");

	// <script>
	fprintf(fp, "\n");
	fprintf(fp, "<script language=\"javascript\" type=\"text/javascript\">\n\n");

	fprintf(fp, "var Pline = [\n");
	for (int page = 0; page <= ev3d_page; page++) {
		ndata[page] = 0;
		if (!ev3d_ndata[page]) continue;
		fprintf(fp, "[\n");
		for (int n = 0; n < ev3d_ndata[page]; n++) {
			ev3d_t *ptr = &ev3d_data[page][n];
			//printf("%d %d %d\n", page, n, ptr->idx);
			if      (ptr->idx == 2) {
				// line -> 1 line
				fprintf(fp, fmt,
					ptr->xvt[0], ptr->yvt[0], ptr->zvt[0],
					ptr->xvt[1], ptr->yvt[1], ptr->zvt[1],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				ndata[page]++;
			}
			else if (ptr->idx == 3) {
				// triangle -> 3 lines
				fprintf(fp, fmt,
					ptr->xvt[0], ptr->yvt[0], ptr->zvt[0],
					ptr->xvt[1], ptr->yvt[1], ptr->zvt[1],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				fprintf(fp, fmt,
					ptr->xvt[1], ptr->yvt[1], ptr->zvt[1],
					ptr->xvt[2], ptr->yvt[2], ptr->zvt[2],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				fprintf(fp, fmt,
					ptr->xvt[2], ptr->yvt[2], ptr->zvt[2],
					ptr->xvt[0], ptr->yvt[0], ptr->zvt[0],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				ndata[page] += 3;
			}
			else if (ptr->idx == 4) {
				// quadrangle -> 4 lines
				fprintf(fp, fmt,
					ptr->xvt[0], ptr->yvt[0], ptr->zvt[0],
					ptr->xvt[1], ptr->yvt[1], ptr->zvt[1],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				fprintf(fp, fmt,
					ptr->xvt[1], ptr->yvt[1], ptr->zvt[1],
					ptr->xvt[2], ptr->yvt[2], ptr->zvt[2],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				fprintf(fp, fmt,
					ptr->xvt[2], ptr->yvt[2], ptr->zvt[2],
					ptr->xvt[3], ptr->yvt[3], ptr->zvt[3],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				fprintf(fp, fmt,
					ptr->xvt[3], ptr->yvt[3], ptr->zvt[3],
					ptr->xvt[0], ptr->yvt[0], ptr->zvt[0],
					ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				ndata[page] += 4;
			}
		}
		fprintf(fp, (page < ev3d_page) ? "],\n" : "]\n");
	}
	fprintf(fp, "];\n");

	int npage = 0;
	for (int page = 0; page <= ev3d_page; page++) {
		if (ndata[page] > 0) {
			npage++;
		}
	}

	fprintf(fp, "var Title = [\n");
	for (int page = 0; page <= ev3d_page; page++) {
		if (!ev3d_ndata[page]) continue;
		fprintf(fp, "[\n");
		for (int n = 0; n < ev3d_ndata[page]; n++) {
			ev3d_t *ptr = &ev3d_data[page][n];
			if (ptr->idx == -4) {
				//fprintf(fp, "[\"%s\",%g,%d,%d,%d,\"%s\"],\n", ptr->str, ptr->hgt, ptr->rgb[0], ptr->rgb[1], ptr->rgb[2], fontname[ptr->fontname]);
				fprintf(fp, "[\"%s\"],\n", ptr->str);
			}
		}
		fprintf(fp, (page < ev3d_page) ? "],\n" : "]\n");
	}
	fprintf(fp, "];\n");

	// bounding box
	fprintf(fp, "var boundingBox = [\n");
	for (int page = 0; page <= ev3d_page; page++) {
		if (!ev3d_ndata[page]) continue;
		double xmin = ev3d_data[page][0].xvt[0];
		double ymin = ev3d_data[page][0].yvt[0];
		double zmin = ev3d_data[page][0].zvt[0];
		double xmax = xmin;
		double ymax = ymin;
		double zmax = zmin;
		for (int n = 0; n < ev3d_ndata[page]; n++) {
			ev3d_t *ptr = &ev3d_data[page][n];
			int nvtx = (ptr->idx > 0) ? ptr->idx : 0;
			for (int m = 0; m < nvtx; m++) {
				if (ptr->xvt[m] < xmin) xmin = ptr->xvt[m];
				if (ptr->xvt[m] > xmax) xmax = ptr->xvt[m];
				if (ptr->yvt[m] < ymin) ymin = ptr->yvt[m];
				if (ptr->yvt[m] > ymax) ymax = ptr->yvt[m];
				if (ptr->zvt[m] < zmin) zmin = ptr->zvt[m];
				if (ptr->zvt[m] > zmax) zmax = ptr->zvt[m];
			}
		}
		fprintf(fp, "[%g,%g,%g,%g,%g,%g]%s\n",
			xmin, ymin, zmin, xmax, ymax, zmax, (page < ev3d_page ? "," : ""));
	}
	fprintf(fp, "];\n");
	fprintf(fp, "\n");

	// variables
	fprintf(fp, "var Width = %d;\n", ev3d_html_width);
	fprintf(fp, "var Height = %d;\n", ev3d_html_height);
	fprintf(fp, "var Mousedrag = 0;\n");
	fprintf(fp, "var Fzoom = 1;\n");
	fprintf(fp, "var Xoff = 0;\n");
	fprintf(fp, "var Yoff = 0;\n");
	fprintf(fp, "var Xold = 0;\n");
	fprintf(fp, "var Yold = 0;\n");
	fprintf(fp, "var Theta0 = %.3f;\n", ev3d_html_theta);
	fprintf(fp, "var Phi0 = %.3f;\n", ev3d_html_phi);
	fprintf(fp, "var Theta = Theta0;\n");
	fprintf(fp, "var Phi = Phi0;\n");
	fprintf(fp, "var Dtor = Math.PI / 180;\n");
	fprintf(fp, "var ipage = 0;\n");
	fprintf(fp, "var npage = %d;\n", npage);
	fprintf(fp, "var canvas, context;\n");
	fprintf(fp, "\n");

	// plot
	fprintf(fp, "function plot() {\n");
	fprintf(fp, "	document.getElementById(\"page\").innerHTML = ((npage > 0) ? (ipage + 1) : 0) + \"/\" + npage;\n");
	fprintf(fp, "	if ((ipage < 0) || (npage <= 0) || (ipage >= npage)) return;\n");
	fprintf(fp, "\n");
	fprintf(fp, "	context.fillStyle = \"rgb(255, 255, 255)\";\n");
	fprintf(fp, "	context.fillRect(0, 0, Width, Height);\n");
	fprintf(fp, "\n");

	fprintf(fp, "	var sint = Math.sin(Theta * Dtor);\n");
	fprintf(fp, "	var cost = Math.cos(Theta * Dtor);\n");
	fprintf(fp, "	var sinp = Math.sin(Phi   * Dtor);\n");
	fprintf(fp, "	var cosp = Math.cos(Phi   * Dtor);\n");
	fprintf(fp, "	var xu = [-sinp, +cosp, 0];\n");
	fprintf(fp, "	var zu = [-cost * cosp, -cost * sinp, +sint];\n");
	fprintf(fp, "	var xmin = boundingBox[ipage][0];\n");
	fprintf(fp, "	var ymin = boundingBox[ipage][1];\n");
	fprintf(fp, "	var zmin = boundingBox[ipage][2];\n");
	fprintf(fp, "	var xmax = boundingBox[ipage][3];\n");
	fprintf(fp, "	var ymax = boundingBox[ipage][4];\n");
	fprintf(fp, "	var zmax = boundingBox[ipage][5];\n");
	fprintf(fp, "	var xc = (xmin + xmax) / 2.0;\n");
	fprintf(fp, "	var yc = (ymin + ymax) / 2.0;\n");
	fprintf(fp, "	var zc = (zmin + zmax) / 2.0;\n");
	fprintf(fp, "	var dspan = Math.sqrt(\n");
	fprintf(fp, "		Math.pow(xmax - xmin, 2) +\n");
	fprintf(fp, "		Math.pow(ymax - ymin, 2) +\n");
	fprintf(fp, "		Math.pow(zmax - zmin, 2));\n");
	fprintf(fp, "	var xp0 = (xu[0] * xc) + (xu[1] * yc) + (xu[2] * zc);\n");
	fprintf(fp, "	var yp0 = (zu[0] * xc) + (zu[1] * yc) + (zu[2] * zc);\n");
	fprintf(fp, "	var fctr = Fzoom * Math.min(Width, Height) / dspan;\n");
	fprintf(fp, "	var x0 = Width / 2;\n");
	fprintf(fp, "	var y0 = Height / 2;\n");

	fprintf(fp, "	for (var i = 0; i < Pline[ipage].length; i++) {\n");
	fprintf(fp, "		context.beginPath();\n");
	fprintf(fp, "		var x1 = Pline[ipage][i][0];\n");
	fprintf(fp, "		var y1 = Pline[ipage][i][1];\n");
	fprintf(fp, "		var z1 = Pline[ipage][i][2];\n");
	fprintf(fp, "		var x2 = Pline[ipage][i][3];\n");
	fprintf(fp, "		var y2 = Pline[ipage][i][4];\n");
	fprintf(fp, "		var z2 = Pline[ipage][i][5];\n");
	fprintf(fp, "		var r  = Pline[ipage][i][6];\n");
	fprintf(fp, "		var g  = Pline[ipage][i][7];\n");
	fprintf(fp, "		var b  = Pline[ipage][i][8];\n");
	fprintf(fp, "		var xp1 = (xu[0] * x1) + (xu[1] * y1) + (xu[2] * z1);\n");
	fprintf(fp, "		var yp1 = (zu[0] * x1) + (zu[1] * y1) + (zu[2] * z1);\n");
	fprintf(fp, "		var xp2 = (xu[0] * x2) + (xu[1] * y2) + (xu[2] * z2);\n");
	fprintf(fp, "		var yp2 = (zu[0] * x2) + (zu[1] * y2) + (zu[2] * z2);\n");
	fprintf(fp, "		xp1 =        + (x0 + fctr * (xp1 - xp0) + Xoff);\n");
	fprintf(fp, "		yp1 = Height - (y0 + fctr * (yp1 - yp0) + Yoff);\n");
	fprintf(fp, "		xp2 =        + (x0 + fctr * (xp2 - xp0) + Xoff);\n");
	fprintf(fp, "		yp2 = Height - (y0 + fctr * (yp2 - yp0) + Yoff);\n");
	fprintf(fp, "		context.moveTo(xp1, yp1);\n");
	fprintf(fp, "		context.lineTo(xp2, yp2);\n");
	fprintf(fp, "		context.strokeStyle = \"rgb(\" + r + \",\" + g + \",\" + b + \")\";\n");
	fprintf(fp, "		context.stroke();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "\n");

	// XYZ arrows
	fprintf(fp, "	// XYZ arrows\n");
	fprintf(fp, "	context.font = \"%dpx %s\";\n", ev3d_html_fontsize, ev3d_html_fontname);
	fprintf(fp, "	context.strokeStyle = 'rgb(0, 0, 0)';\n");
	fprintf(fp, "	context.fillStyle = 'rgb(0, 0, 0)';\n");
	fprintf(fp, "	var da = (Width + Height) / 20;\n");
	fprintf(fp, "	var xa = da + 20;\n");
	fprintf(fp, "	var ya = Height - xa;\n");
	fprintf(fp, "	var dx = [-sinp, +cosp, 0];\n");
	fprintf(fp, "	var dy = [-cost * cosp, -cost * sinp, +sint];\n");
	fprintf(fp, "	var ca = ['X', 'Y', 'Z'];\n");
	fprintf(fp, "	for (var i = 0; i < 3; i++) {\n");
	fprintf(fp, "		context.beginPath();\n");
	fprintf(fp, "		context.moveTo(xa, ya);\n");
	fprintf(fp, "		context.lineTo(xa + da * dx[i], ya - da * dy[i]);\n");
	fprintf(fp, "		context.stroke();\n");
	fprintf(fp, "		context.fillText(ca[i], xa + 1.2 * da * dx[i] - 7, ya - 1.2 * da * dy[i] + 7);\n");
	fprintf(fp, "	}\n");

	// title
	fprintf(fp, "\n");
	fprintf(fp, "	// title\n");
	fprintf(fp, "	var xstr = 5;\n");
	fprintf(fp, "	var ystr = 5;\n");
	fprintf(fp, "	var hgt = %d;\n", ev3d_html_fontsize);
	fprintf(fp, "	context.font = \"%dpx %s\";\n", ev3d_html_fontsize, ev3d_html_fontname);
	fprintf(fp, "	context.fillStyle = \'rgb(0, 0, 0)\';\n");
	fprintf(fp, "	for (var i = 0; i < Title[ipage].length; i++) {\n");
	fprintf(fp, "		context.fillText(Title[ipage][i], xstr, ystr + 1.2 * hgt);\n");
	fprintf(fp, "		ystr += 1.2 * hgt;\n");
	fprintf(fp, "	}\n");

	// number of lines
	fprintf(fp, "\n");
	fprintf(fp, "	// number of lines\n");
	fprintf(fp, "	context.font = \"%dpx %s\";\n", ev3d_html_fontsize, ev3d_html_fontname);
	fprintf(fp, "	context.fillStyle = 'rgb(0, 0, 0)';\n");
	fprintf(fp, "	context.fillText(Pline[ipage].length + \" lines\", xstr, Height - 6);\n");

	fprintf(fp, "}\n");

	// init
	fprintf(fp, "function init() {\n");
	fprintf(fp, "	canvas = document.getElementById(\"ev3d\");\n");
	fprintf(fp, "	context = canvas.getContext(\"2d\");\n");
	fprintf(fp, "	Xoff = 0;\n");
	fprintf(fp, "	Yoff = 0;\n");
	fprintf(fp, "	Fzoom = 1;\n");
	fprintf(fp, "	Theta = Theta0;\n");
	fprintf(fp, "	Phi = Phi0;\n");
	fprintf(fp, "	first();\n");
	fprintf(fp, "}\n");

	// X
	fprintf(fp, "function xview() {\n");
	fprintf(fp, "	Theta = 90;\n");
	fprintf(fp, "	Phi = 0;\n");
	fprintf(fp, "	Mousedrag = 0;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// Y
	fprintf(fp, "function yview() {\n");
	fprintf(fp, "	Theta = 90;\n");
	fprintf(fp, "	Phi = 90;\n");
	fprintf(fp, "	Mousedrag = 0;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// Z
	fprintf(fp, "function zview() {\n");
	fprintf(fp, "	Theta = 0;\n");
	fprintf(fp, "	Phi = -90;\n");
	fprintf(fp, "	Mousedrag = 0;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// first
	fprintf(fp, "function first() {\n");
	fprintf(fp, "	ipage = 0;\n");
	fprintf(fp, "	Mousedrag = 0;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// prev
	fprintf(fp, "function prev() {\n");
	fprintf(fp, "	if (ipage > 0) {\n");
	fprintf(fp, "		ipage--;\n");
	fprintf(fp, "		Mousedrag = 0;\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");

	// next
	fprintf(fp, "function next() {\n");
	fprintf(fp, "	if (ipage < npage - 1) {\n");
	fprintf(fp, "		ipage++;\n");
	fprintf(fp, "		Mousedrag = 0;\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");

	// last
	fprintf(fp, "function last() {\n");
	fprintf(fp, "	ipage = npage - 1;\n");
	fprintf(fp, "	Mousedrag = 0;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// mouse
	fprintf(fp, "document.body.onmousedown = function(e) {\n");
	fprintf(fp, "	Mousedrag = 1;\n");
	fprintf(fp, "	Xold = e.clientX;\n");
	fprintf(fp, "	Yold = e.clientY;\n");
	fprintf(fp, "}\n");
	fprintf(fp, "document.body.onmouseup = function(e) {\n");
	fprintf(fp, "	Mousedrag = 0;\n");
	fprintf(fp, "}\n");
	fprintf(fp, "document.body.onmousemove = function(e) {\n");
	fprintf(fp, "	if (Mousedrag == 0) return;\n");
	fprintf(fp, "	if (e.shiftKey) {\n");
	fprintf(fp, "		Xoff += e.clientX - Xold;\n");
	fprintf(fp, "		Yoff -= e.clientY - Yold;\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "	else {\n");
	fprintf(fp, "		var fctr = 360.0 / (Width + Height);\n");
	fprintf(fp, "		Theta -= fctr * (e.clientY - Yold);\n");
	fprintf(fp, "		Phi   -= fctr * (e.clientX - Xold);\n");
	fprintf(fp, "		if (Theta < 0  ) Theta = 0;\n");
	fprintf(fp, "		if (Theta > 180) Theta = 180;\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "	Xold = e.clientX;\n");
	fprintf(fp, "	Yold = e.clientY;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// wheel
	fprintf(fp, "var mousewheelevt = (navigator.userAgent.indexOf(\"Firefox\") != -1) ? \"DOMMouseScroll\" : \"mousewheel\";\n");
	fprintf(fp, "document.body.addEventListener(mousewheelevt, wheel, false);\n");
	fprintf(fp, "function wheel(e) {\n");
	fprintf(fp, "	var evt = window.event || e;\n");
	fprintf(fp, "	var delta = evt.detail ? evt.detail * (-120) : evt.wheelDelta;\n");
	fprintf(fp, "	Fzoom *= (delta > 0) ? 1.2 : (1 / 1.2);\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// zoom
	fprintf(fp, "function zoomin() {\n");
	fprintf(fp, "	Fzoom *= 1.2;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function zoomout() {\n");
	fprintf(fp, "	Fzoom /= 1.2;\n");
	fprintf(fp, "	plot();\n");
	fprintf(fp, "}\n");

	// end of <script>
	fprintf(fp, "</script>\n");
	fprintf(fp, "\n");

	// button
	fprintf(fp, "<div id=\"page\">%d/%d</div>\n", (ev3d_page >= 0 ? 1 : 0), ev3d_page + 1);
	fprintf(fp, "<input type=\"button\" value=\"|<\" onClick=\"first()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"<\" onClick=\"prev()\">\n");
	fprintf(fp, "<input type=\"button\" value=\">\" onClick=\"next()\">\n");
	fprintf(fp, "<input type=\"button\" value=\">|\" onClick=\"last()\">&nbsp;&nbsp;\n");
	fprintf(fp, "<input type=\"button\" value=\"X\" onClick=\"xview()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"Y\" onClick=\"yview()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"Z\" onClick=\"zview()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"0\" onClick=\"init()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"+\" onClick=\"zoomin()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"-\" onClick=\"zoomout()\">\n");
	fprintf(fp, "<br>\n");

	// <canvas>
	fprintf(fp, "<canvas id=\"ev3d\" style=\"border: 1px solid;\" width=\"%d\" height=\"%d\">\n", ev3d_html_width, ev3d_html_height);
	fprintf(fp, "</canvas>\n");

	// date
	time_t now;
	time(&now);
	fprintf(fp, "<br>\n");
	fprintf(fp, "%s", ctime(&now));
	fprintf(fp, "<br>\n");
	fprintf(fp, "\n");
/*
	// [button]
	fprintf(fp, "<p>\n");
	fprintf(fp, "[button]<br>\n");
	fprintf(fp, "&lt;: previous page<br>\n");
	fprintf(fp, "&gt;: next page<br>\n");
	fprintf(fp, "0: initialize view<br>\n");
	fprintf(fp, "X: X-direction view<br>\n");
	fprintf(fp, "Y: Y-direction view<br>\n");
	fprintf(fp, "Z: Z-direction view<br>\n");
	fprintf(fp, "</p>\n");
	fprintf(fp, "\n");

	// [mouse]
	fprintf(fp, "<p>\n");
	fprintf(fp, "[mouse]<br>\n");
	fprintf(fp, "left-drag: rotation<br>\n");
	fprintf(fp, "wheel: zoom<br>\n");
	fprintf(fp, "right-click: save image<br>\n");
	fprintf(fp, "</p>\n");
	fprintf(fp, "\n");
*/
	// tailor
	fprintf(fp, "</body>\n");
	fprintf(fp, "</html>\n");
}


// file (option)
void ev3d_file(int type, const char fn[], int binary)
{
	ev3d_type = type;

	strcpy(ev3d_fn, fn);

	ev3d_binary = binary;
}


// output
void ev3d_output(void)
{
	// open file
	FILE *fp;
	if ((fp = fopen(ev3d_fn, (ev3d_binary ? "wb" : "w"))) == NULL) {
		fprintf(stderr, "*** %s file open error.\n", ev3d_fn);
		exit(1);
	}

	// output
	if (!ev3d_type) {
		ev3d_end_html(fp, ev3d_fn);
	}
	else if (ev3d_binary) {
		ev3d_end_ev3_binary(fp);
	}
	else {
		ev3d_end_ev3_ascii(fp);
	}

	// close file
	fclose(fp);
}


// set color
void ev3d_setColor(unsigned char r, unsigned char g, unsigned char b)
{
	ev3d_rgb[0] = r;
	ev3d_rgb[1] = g;
	ev3d_rgb[2] = b;
}


// set color
void ev3d_setColorA(const unsigned char rgb[])
{
	memcpy(ev3d_rgb, rgb, 3 * sizeof(unsigned char));
}


// f=0-1 -> (r, g, b)
// color = 0/1 : gray/color
void ev3d_setColorV(double f, int color)
{
	const unsigned char rgb[20][3] = {
		{  0,   0, 255},
		{  0, 110, 255},
		{  0, 140, 255},
		{  0, 170, 255},
		{  0, 200, 255},
		{  0, 230, 255},
		{  0, 255, 255},
		{  0, 255, 200},
		{  0, 255, 160},
		{  0, 255, 130},
		{  0, 255,   0},
		{130, 255,   0},
		{180, 255,   0},
		{220, 255,   0},
		{255, 255,   0},
		{255, 230,   0},
		{255, 190,   0},
		{255, 150,   0},
		{255, 100,   0},
		{255,   0,   0}
	};

	unsigned char r, g, b;
	if (color) {
		const int div = sizeof(rgb) / 3;
		int icolor = (int)(f * div);
		if (icolor < 0) icolor = 0;
		if (icolor >= div) icolor = div - 1;

		r = rgb[icolor][0];
		g = rgb[icolor][1];
		b = rgb[icolor][2];
	}
	else {
		if (f < 0) f = 0;
		if (f > 1) f = 1;
		r = g = b = (unsigned char)(255 * (1 - f));
	}
	ev3d_setColor(r, g, b);
}


// draw a line (primitive)
void ev3d_drawLine(double x1, double y1, double z1, double x2, double y2, double z2)
{
	// alloc
	alloc_inc();

	const int ndata = ev3d_ndata[ev3d_page];

	// index
	ev3d_data[ev3d_page][ndata].idx = 2;
	ev3d_data[ev3d_page][ndata].num = ev3d_num;

	// vertex
	ev3d_data[ev3d_page][ndata].xvt[0] = x1;
	ev3d_data[ev3d_page][ndata].yvt[0] = y1;
	ev3d_data[ev3d_page][ndata].zvt[0] = z1;
	ev3d_data[ev3d_page][ndata].xvt[1] = x2;
	ev3d_data[ev3d_page][ndata].yvt[1] = y2;
	ev3d_data[ev3d_page][ndata].zvt[1] = z2;

	// color
	memcpy(ev3d_data[ev3d_page][ndata].rgb, ev3d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev3d_ndata[ev3d_page]++;
}


// draw a triangle
void ev3d_drawTriangle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
	ev3d_drawLine(x1, y1, z1, x2, y2, z2);
	ev3d_drawLine(x2, y2, z2, x3, y3, z3);
	ev3d_drawLine(x3, y3, z3, x1, y1, z1);
}


// draw a quadrangle
void ev3d_drawQuadrangle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)
{
	ev3d_drawLine(x1, y1, z1, x2, y2, z2);
	ev3d_drawLine(x2, y2, z2, x3, y3, z3);
	ev3d_drawLine(x3, y3, z3, x4, y4, z4);
	ev3d_drawLine(x4, y4, z4, x1, y1, z1);
}


// fill a triangle (primitive)
void ev3d_fillTriangle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
	// alloc
	alloc_inc();

	const int ndata = ev3d_ndata[ev3d_page];

	// index
	ev3d_data[ev3d_page][ndata].idx = 3;
	ev3d_data[ev3d_page][ndata].num = ev3d_num;

	// vertex
	ev3d_data[ev3d_page][ndata].xvt[0] = x1;
	ev3d_data[ev3d_page][ndata].yvt[0] = y1;
	ev3d_data[ev3d_page][ndata].zvt[0] = z1;
	ev3d_data[ev3d_page][ndata].xvt[1] = x2;
	ev3d_data[ev3d_page][ndata].yvt[1] = y2;
	ev3d_data[ev3d_page][ndata].zvt[1] = z2;
	ev3d_data[ev3d_page][ndata].xvt[2] = x3;
	ev3d_data[ev3d_page][ndata].yvt[2] = y3;
	ev3d_data[ev3d_page][ndata].zvt[2] = z3;

	// color
	memcpy(ev3d_data[ev3d_page][ndata].rgb, ev3d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev3d_ndata[ev3d_page]++;
}


// draw a quadrangle (primitive)
void ev3d_fillQuadrangle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)
{
	// alloc
	alloc_inc();

	const int ndata = ev3d_ndata[ev3d_page];

	// index
	ev3d_data[ev3d_page][ndata].idx = 4;
	ev3d_data[ev3d_page][ndata].num = ev3d_num;

	// vertex
	ev3d_data[ev3d_page][ndata].xvt[0] = x1;
	ev3d_data[ev3d_page][ndata].yvt[0] = y1;
	ev3d_data[ev3d_page][ndata].zvt[0] = z1;
	ev3d_data[ev3d_page][ndata].xvt[1] = x2;
	ev3d_data[ev3d_page][ndata].yvt[1] = y2;
	ev3d_data[ev3d_page][ndata].zvt[1] = z2;
	ev3d_data[ev3d_page][ndata].xvt[2] = x3;
	ev3d_data[ev3d_page][ndata].yvt[2] = y3;
	ev3d_data[ev3d_page][ndata].zvt[2] = z3;
	ev3d_data[ev3d_page][ndata].xvt[3] = x4;
	ev3d_data[ev3d_page][ndata].yvt[3] = y4;
	ev3d_data[ev3d_page][ndata].zvt[3] = z4;

	// color
	memcpy(ev3d_data[ev3d_page][ndata].rgb, ev3d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev3d_ndata[ev3d_page]++;
}


// draw a rectangle
void ev3d_drawRectangle(char dir, double c0, double p1, double q1, double p2, double q2)
{
	dir = (char)toupper(dir);
	if      (dir == 'X') {
		ev3d_drawQuadrangle(c0, p1, q1, c0, p2, q1, c0, p2, q2, c0, p1, q2);
	}
	else if (dir == 'Y') {
		ev3d_drawQuadrangle(q1, c0, p1, q1, c0, p2, q2, c0, p2, q2, c0, p1);
	}
	else if (dir == 'Z') {
		ev3d_drawQuadrangle(p1, q1, c0, p2, q1, c0, p2, q2, c0, p1, q2, c0);
	}
}


// fill a rectangle
void ev3d_fillRectangle(char dir, double c0, double p1, double q1, double p2, double q2)
{
	dir = (char)toupper(dir);
	if      (dir == 'X') {
		ev3d_fillQuadrangle(c0, p1, q1, c0, p2, q1, c0, p2, q2, c0, p1, q2);
	}
	else if (dir == 'Y') {
		ev3d_fillQuadrangle(q1, c0, p1, q1, c0, p2, q2, c0, p2, q2, c0, p1);
	}
	else if (dir == 'Z') {
		ev3d_fillQuadrangle(p1, q1, c0, p2, q1, c0, p2, q2, c0, p1, q2, c0);
	}
}


// draw a box
void ev3d_drawBox(double x1, double y1, double z1, double x2, double y2, double z2)
{
	// floor/ceil
	ev3d_drawRectangle('Z', z1, x1, y1, x2, y2);
	ev3d_drawRectangle('Z', z2, x1, y1, x2, y2);

	// pillar
	ev3d_drawLine(x1, y1, z1, x1, y1, z2);
	ev3d_drawLine(x2, y1, z1, x2, y1, z2);
	ev3d_drawLine(x2, y2, z1, x2, y2, z2);
	ev3d_drawLine(x1, y2, z1, x1, y2, z2);
}


// fill a box
void ev3d_fillBox(double x1, double y1, double z1, double x2, double y2, double z2)
{
	ev3d_fillRectangle('X', x1, y1, z1, y2, z2);
	ev3d_fillRectangle('X', x2, y1, z1, y2, z2);
	ev3d_fillRectangle('Y', y1, z1, x1, z2, x2);
	ev3d_fillRectangle('Y', y2, z1, x1, z2, x2);
	ev3d_fillRectangle('Z', z1, x1, y1, x2, y2);
	ev3d_fillRectangle('Z', z2, x1, y1, x2, y2);
}


// ellipse
static void ellipse(int fill, char dir, double c0, double p1, double q1, double p2, double q2, int div)
{
	const double pi = 4 * atan(1);
	dir = (char)toupper(dir);

	for (int i = 0; i < div; i++) {
		double p0 = (p1 + p2) / 2;
		double q0 = (q1 + q2) / 2;
		double pr = fabs(p1 - p2) / 2;
		double qr = fabs(q1 - q2) / 2;
		double a1 = (i    ) * (2 * pi) / div;
		double a2 = (i + 1) * (2 * pi) / div;
		double x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0;
		if      (dir == 'X') {
			x0 = c0;
			y0 = p0;
			z0 = q0;
			x1 = c0;
			y1 = p0 + pr * cos(a1);
			z1 = q0 + qr * sin(a1);
			x2 = c0;
			y2 = p0 + pr * cos(a2);
			z2 = q0 + qr * sin(a2);
		}
		else if (dir == 'Y') {
			y0 = c0;
			z0 = p0;
			x0 = q0;
			y1 = c0;
			z1 = p0 + pr * cos(a1);
			x1 = q0 + qr * sin(a1);
			y2 = c0;
			z2 = p0 + pr * cos(a2);
			x2 = q0 + qr * sin(a2);
		}
		else if (dir == 'Z') {
			z0 = c0;
			x0 = p0;
			y0 = q0;
			z1 = c0;
			x1 = p0 + pr * cos(a1);
			y1 = q0 + qr * sin(a1);
			z2 = c0;
			x2 = p0 + pr * cos(a2);
			y2 = q0 + qr * sin(a2);
		}
		else {
			return;
		}

		if (fill) {
			ev3d_fillTriangle(x0, y0, z0, x1, y1, z1, x2, y2, z2);
		}
		else {
			ev3d_drawLine(x1, y1, z1, x2, y2, z2);
		}
	}
}


// draw an ellipse
void ev3d_drawEllipse(char dir, double c0, double p1, double q1, double p2, double q2, int div)
{
	ellipse(0, dir, c0, p1, q1, p2, q2, div);
}


// fill an ellipse
void ev3d_fillEllipse(char dir, double c0, double p1, double q1, double p2, double q2, int div)
{
	ellipse(1, dir, c0, p1, q1, p2, q2, div);
}

/*
// string (primitive)
static void string(int idx, double x, double y, double z, double h, const char str[])
{
	// alloc
	alloc_inc();

	const int ndata = ev3d_ndata[ev3d_page];

	// index
	ev3d_data[ev3d_page][ndata].idx = idx;
	ev3d_data[ev3d_page][ndata].num = ev3d_num;

	// position (string)
	if (idx == -3) {
		ev3d_data[ev3d_page][ndata].xvt[0] = x;
		ev3d_data[ev3d_page][ndata].yvt[0] = y;
		ev3d_data[ev3d_page][ndata].zvt[0] = z;
	}

	// height
	ev3d_data[ev3d_page][ndata].hgt = h;

	// color
	memcpy(ev3d_data[ev3d_page][ndata].rgb, ev3d_rgb, 3 * sizeof(unsigned char));

	// string
	ev3d_data[ev3d_page][ndata].str = (char *)malloc((strlen(str) + 1) * sizeof(char));
	strcpy(ev3d_data[ev3d_page][ndata].str, str);

	// fontname
	//ev3d_data[ev3d_page][ndata].fontname = ev3d_fontname;

	// data++
	ev3d_ndata[ev3d_page]++;
}

// draw string (3D) (Ver.1.2, not supported)
void ev3d_drawString(double x, double y, double z, double h, const char str[])
{
	string(-3, x, y, z, h, str);
}
*/


// string (primitive)
static void string(int idx, const char str[])
{
	// alloc
	alloc_inc();

	const int ndata = ev3d_ndata[ev3d_page];

	// index
	ev3d_data[ev3d_page][ndata].idx = idx;
	//ev3d_data[ev3d_page][ndata].num = ev3d_num;

	// height
	//ev3d_data[ev3d_page][ndata].hgt = h;

	// string
	ev3d_data[ev3d_page][ndata].str = (char *)malloc((strlen(str) + 1) * sizeof(char));
	strcpy(ev3d_data[ev3d_page][ndata].str, str);

	// data++
	ev3d_ndata[ev3d_page]++;
}


// draw title (2D)
void ev3d_drawTitle(const char str[])
{
	string(-4, str);
}


// window size (HTML only)
void ev3d_html_size(int width, int height)
{
	ev3d_html_width  = width;
	ev3d_html_height = height;
}


// set angle (HTML only)
void ev3d_html_angle(double theta, double phi)
{
	ev3d_html_theta = theta;
	ev3d_html_phi = phi;
}


// set fontname (HTML only)
void ev3d_html_font(int fontname, int fontsize)
{
	if      (fontname == 0) {
		strcpy(ev3d_html_fontname, "sansserif");
	}
	else if (fontname == 1) {
		strcpy(ev3d_html_fontname, "serif");
	}
	else if (fontname == 2) {
		strcpy(ev3d_html_fontname, "monospace");
	}
	ev3d_html_fontsize = fontsize;
}


// index
void ev3d_index(int num)
{
	ev3d_num = num;
}
