/*
ev2d.c

output: ev2 or HTML(2D)

ev2 format:
-1 W H                    new page, width and height of the page
-2 R G B                  set color, RGB(0-255)
-3 X Y H                  draw a string: start position, height
characters                string
2 X1 Y1 X2 Y2             draw a line
3 X1 Y1 X2 Y2 X3 Y3       fill a triangle
4 X1 Y1 X2 Y2 X3 Y3 X4 Y4 fill a quadrangle
21 X1 Y1 X2 Y2            draw an ellipse
22 X1 Y1 X2 Y2            fill an ellipse
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct {
	int           idx;  // 2:line, 3:triangle(fill), 4:quadrangle(fill), 21:ellipse(line), 22:ellipse(fill), -3:characters
	double        xvt[4];
	double        yvt[4];
	unsigned char rgb[3];
	double        hgt;
	char          *str;
} ev2d_t;

int ev2d_type = 1;
char ev2d_fn[BUFSIZ] = "ev.ev2";
int ev2d_page = -1;
double ev2d_width = 0;
double ev2d_height = 0;
unsigned char ev2d_rgb[3] = {0, 0, 0};
int *ev2d_ndata;
ev2d_t **ev2d_data;
const int ev2d_increment = 10000;


// alloc increment
static void alloc_inc(void)
{
	if ((ev2d_ndata[ev2d_page]) > 0 && (ev2d_ndata[ev2d_page] % ev2d_increment == 0)) {
		ev2d_data[ev2d_page] = (ev2d_t *)realloc(ev2d_data[ev2d_page], (ev2d_ndata[ev2d_page] + ev2d_increment) * sizeof(ev2d_t));
	}
}


// initialize
void ev2d_init(double width, double height)
{
	// page
	ev2d_page = -1;

	// alloc
	ev2d_ndata = (int *)    malloc(1 * sizeof(int));
	ev2d_data  = (ev2d_t **)malloc(1 * sizeof(ev2d_t *));

	// window size
	ev2d_width  = width;
	ev2d_height = height;
}


// new page
void ev2d_newPage(void)
{
	// skip no data page
	if ((ev2d_page >= 0) && !ev2d_ndata[ev2d_page]) return;

	// page++
	ev2d_page++;

	// alloc
	if (ev2d_page > 0) {
		ev2d_ndata = (int *)    realloc(ev2d_ndata, (ev2d_page + 1) * sizeof(int));
		ev2d_data  = (ev2d_t **)realloc(ev2d_data,  (ev2d_page + 1) * sizeof(ev2d_t *));
	}
	ev2d_ndata[ev2d_page] = 0;
	ev2d_data[ev2d_page] = (ev2d_t *)malloc(ev2d_increment * sizeof(ev2d_t));

	// clear color
	ev2d_rgb[0] = ev2d_rgb[1] = ev2d_rgb[2] = 0;
}


// end (ev.ev2)
static void ev2d_end_data(FILE *fp)
{
	for (int page = 0; page <= ev2d_page; page++) {
		if (!ev2d_ndata[page]) continue;
		unsigned char rgb[] = {0, 0, 0};
		fprintf(fp, "-1 %g %g\n", ev2d_width, ev2d_height);
		fprintf(fp, "-2 %u %u %u\n", rgb[0], rgb[1], rgb[2]);
		for (int n = 0; n < ev2d_ndata[page]; n++) {
			ev2d_t *ptr = &ev2d_data[page][n];
			if ((ptr->rgb[0] != rgb[0]) || (ptr->rgb[1] != rgb[1]) || (ptr->rgb[2] != rgb[2])) {
				// change color
				fprintf(fp, "-2 %u %u %u\n", ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				rgb[0] = ptr->rgb[0];
				rgb[1] = ptr->rgb[1];
				rgb[2] = ptr->rgb[2];
			}

			fprintf(fp, "%d ", ptr->idx);
			if      (ptr->idx == 2) {
				fprintf(fp, "%g %g %g %g\n", ptr->xvt[0], ptr->yvt[0], ptr->xvt[1], ptr->yvt[1]);
			}
			else if (ptr->idx == 3) {
				fprintf(fp, "%g %g %g %g %g %g\n", ptr->xvt[0], ptr->yvt[0], ptr->xvt[1], ptr->yvt[1], ptr->xvt[2], ptr->yvt[2]);
			}
			else if (ptr->idx == 4) {
				fprintf(fp, "%g %g %g %g %g %g %g %g\n", ptr->xvt[0], ptr->yvt[0], ptr->xvt[1], ptr->yvt[1], ptr->xvt[2], ptr->yvt[2], ptr->xvt[3], ptr->yvt[3]);
			}
			else if ((ptr->idx == 21) || (ptr->idx == 22)) {
				fprintf(fp, "%g %g %g %g\n", ptr->xvt[0], ptr->yvt[0], ptr->xvt[1], ptr->yvt[1]);
			}
			else if (ptr->idx == -3) {
				fprintf(fp, "%g %g %g\n", ptr->xvt[0], ptr->yvt[0], ptr->hgt);
				fprintf(fp, "%s\n", ptr->str);
			}
		}
	}
}


// end (HTML)
static void ev2d_end_html(FILE *fp, const char fn[])
{
	// number of page
	int npage = 0;
	for (int page = 0; page <= ev2d_page; page++) {
		if (ev2d_ndata[page] > 0) {
			npage++;
		}
	}

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

	fprintf(fp, "var Data = [\n");
	for (int page = 0; page <= ev2d_page; page++) {
		if (!ev2d_ndata[page]) continue;
		fprintf(fp, "[\n");
		for (int n = 0; n < ev2d_ndata[page]; n++) {
			ev2d_t *ptr = &ev2d_data[page][n];
			fprintf(fp, "[%d", ptr->idx);
			if ((ptr->idx == 2) || (ptr->idx == 3) || (ptr->idx == 4)) {
				// line, triangle, quadrangle
				for (int i = 0; i < ptr->idx; i++) {
					fprintf(fp, ",%g,%g", ptr->xvt[i], ptr->yvt[i]);
				}
				fprintf(fp, ",%d,%d,%d]", ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
			}
			else if (ptr->idx == -3) {
				// characters
				fprintf(fp, ",%g,%g", ptr->xvt[0], ptr->yvt[0]);
				fprintf(fp, ",%d,%d,%d", ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
				//fprintf(fp, ",\"%s\",\"%dpx %s\"]", ptr->str, ptr->hgt, ptr->fontname);
				fprintf(fp, ",\"%s\",\"%gpx monospace\"]", ptr->str, ptr->hgt);
			}
			else if ((ptr->idx == 21) || (ptr->idx == 22)) {
				// ellipse
				fprintf(fp, ",%g,%g,%g,%g", ptr->xvt[0], ptr->yvt[0], ptr->xvt[1], ptr->yvt[1]);
				fprintf(fp, ",%d,%d,%d]", ptr->rgb[0], ptr->rgb[1], ptr->rgb[2]);
			}
			fprintf(fp, (n < ev2d_ndata[page] - 1) ? ",\n" : "\n");
		}
		fprintf(fp, (page < ev2d_page) ? "],\n" : "]\n");
	}
	fprintf(fp, "];\n\n");

	fprintf(fp, "var Width = %g;\n", ev2d_width);
	fprintf(fp, "var Height = %g;\n", ev2d_height);
	fprintf(fp, "var ipage = 0;\n");
	fprintf(fp, "var npage = %d;\n", npage);
	fprintf(fp, "var canvas, context;\n");
	fprintf(fp, "var timer = null;\n");
	fprintf(fp, "var interval = 500;\n");
	fprintf(fp, "\n");
	fprintf(fp, "function plot() {\n");
	fprintf(fp, "	document.getElementById(\"page\").innerHTML = ((npage > 0) ? (ipage + 1) : 0) + \"/\" + npage;\n");
	fprintf(fp, "	if ((ipage < 0) || (npage <= 0) || (ipage >= npage)) return;\n");
	fprintf(fp, "\n");
	fprintf(fp, "	for (var n = 0; n < Data[ipage].length; n++) {\n");
	fprintf(fp, "		var idx = Data[ipage][n][0];\n");
	fprintf(fp, "		if ((idx == 2) || (idx == 3) || (idx == 4)) {\n");
	fprintf(fp, "			context.beginPath();\n");
	fprintf(fp, "			context.moveTo(Data[ipage][n][1], Height - Data[ipage][n][2]);\n");
	fprintf(fp, "			for (var i = 1; i < idx; i++) {\n");
	fprintf(fp, "				context.lineTo(Data[ipage][n][2 * i + 1], Height - Data[ipage][n][2 * i + 2]);\n");
	fprintf(fp, "			}\n");
	fprintf(fp, "			context.closePath();\n");
	fprintf(fp, "			var r = Data[ipage][n][2 * idx + 1];\n");
	fprintf(fp, "			var g = Data[ipage][n][2 * idx + 2];\n");
	fprintf(fp, "			var b = Data[ipage][n][2 * idx + 3];\n");
	fprintf(fp, "			context.strokeStyle = \"rgb(\" + r + \",\" + g + \",\" + b + \")\";\n");
	fprintf(fp, "			context.fillStyle = \"rgb(\" + r + \",\" + g + \",\" + b + \")\";\n");
	fprintf(fp, "			if ((idx == 3) || (idx == 4)) {\n");
	fprintf(fp, "				context.fill();\n");
	fprintf(fp, "			}\n");
	fprintf(fp, "			context.stroke();\n");
	fprintf(fp, "		}\n");
	fprintf(fp, "		else if ((idx == 21) || (idx == 22)) {\n");
	fprintf(fp, "			var x1 = Data[ipage][n][1];\n");
	fprintf(fp, "			var y1 = Height - Data[ipage][n][2];\n");
	fprintf(fp, "			var x2 = Data[ipage][n][3];\n");
	fprintf(fp, "			var y2 = Height - Data[ipage][n][4];\n");
	fprintf(fp, "			var x0 = (x1 + x2) / 2;\n");
	fprintf(fp, "			var y0 = (y1 + y2) / 2;\n");
	fprintf(fp, "			var rad = Math.abs((x2 - x1) / 2);\n");
	fprintf(fp, "			var r = Data[ipage][n][5];\n");
	fprintf(fp, "			var g = Data[ipage][n][6];\n");
	fprintf(fp, "			var b = Data[ipage][n][7];\n");
	fprintf(fp, "			context.save();\n");
	fprintf(fp, "			context.translate(x0, y0);\n");
	fprintf(fp, "			context.scale(1, Math.abs((y2 - y1) / (x2 - x1)));\n");
	fprintf(fp, "			context.beginPath();\n");
	fprintf(fp, "			context.arc(0, 0, rad, 0, 2 * Math.PI, true);\n");
	fprintf(fp, "			context.strokeStyle = \"rgb(\" + r + \",\" + g + \",\" + b + \")\";\n");
	fprintf(fp, "			context.fillStyle = \"rgb(\" + r + \",\" + g + \",\" + b + \")\";\n");
	fprintf(fp, "			if (idx == 22) {\n");
	fprintf(fp, "				context.fill();\n");
	fprintf(fp, "			}\n");
	fprintf(fp, "			context.stroke();\n");
	fprintf(fp, "			context.restore();\n");
	fprintf(fp, "		}\n");
	fprintf(fp, "		else if (idx == -3) {\n");
	fprintf(fp, "			var x = Data[ipage][n][1];\n");
	fprintf(fp, "			var y = Height - Data[ipage][n][2];\n");
	fprintf(fp, "			var r = Data[ipage][n][3];\n");
	fprintf(fp, "			var g = Data[ipage][n][4];\n");
	fprintf(fp, "			var b = Data[ipage][n][5];\n");
	fprintf(fp, "			var str = Data[ipage][n][6];\n");
	fprintf(fp, "			var font = Data[ipage][n][7];\n");
	fprintf(fp, "			context.fillStyle = \"rgb(\" + r + \",\" + g + \",\" + b + \")\";\n");
	fprintf(fp, "			context.font = font;\n");
	fprintf(fp, "			context.fillText(str, x, y);\n");
	fprintf(fp, "		}\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function init() {\n");
	fprintf(fp, "	canvas = document.getElementById(\"ev2d\");\n");
	fprintf(fp, "	context = canvas.getContext(\"2d\");\n");
	fprintf(fp, "	first();\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function clear() {\n");
	fprintf(fp, "	context.fillStyle = \"rgb(255, 255, 255)\";\n");
	fprintf(fp, "	context.fillRect(0, 0, Width, Height);\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function first() {\n");
	fprintf(fp, "	if (npage > 0) {\n");
	fprintf(fp, "		ipage = 0;\n");
	fprintf(fp, "		clear();\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function prev() {\n");
	fprintf(fp, "	if (ipage > 0) {\n");
	fprintf(fp, "		ipage--;\n");
	fprintf(fp, "		clear();\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function next() {\n");
	fprintf(fp, "	if (ipage < npage - 1) {\n");
	fprintf(fp, "		ipage++;\n");
	fprintf(fp, "		clear();\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function last() {\n");
	fprintf(fp, "	if (npage > 0) {\n");
	fprintf(fp, "		ipage = npage - 1;\n");
	fprintf(fp, "		clear();\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function allpage() {\n");
	fprintf(fp, "	clear();\n");
	fprintf(fp, "	for (var i = 0; i < npage; i++) {\n");
	fprintf(fp, "		ipage = i;\n");
	fprintf(fp, "		plot();\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "function play() {\n");
	fprintf(fp, "	if (timer == null) {\n");
	fprintf(fp, "		timer = setInterval(function() {\n");
	fprintf(fp, "			clear();\n");
	fprintf(fp, "			ipage++;\n");
	fprintf(fp, "			if (ipage == npage) ipage = 0;\n");
	fprintf(fp, "			plot();\n");
	fprintf(fp, "		}, interval);\n");
	fprintf(fp, "		document.getElementById(\"play\").value = \"pause\";\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "	else {\n");
	fprintf(fp, "		clearInterval(timer);\n");
	fprintf(fp, "		timer = null;\n");
	fprintf(fp, "		document.getElementById(\"play\").value = \"play\";\n");
	fprintf(fp, "	}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "</script>\n");
	fprintf(fp, "\n");

	// button
	fprintf(fp, "<div id=\"page\">%d/%d</div>\n", (ev2d_page > 0 ? 1 : 0), (ev2d_page > 0 ? ev2d_page : 0));
	fprintf(fp, "<input type=\"button\" value=\"|<\" onClick=\"first()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"<\" onClick=\"prev()\">\n");
	fprintf(fp, "<input type=\"button\" value=\">\" onClick=\"next()\">\n");
	fprintf(fp, "<input type=\"button\" value=\">|\" onClick=\"last()\">&nbsp;&nbsp;\n");
	fprintf(fp, "<input type=\"button\" value=\"all\" onClick=\"allpage()\">\n");
	fprintf(fp, "<input type=\"button\" value=\"play\" onClick=\"play()\" id=\"play\">\n");
	fprintf(fp, "<br>\n");

	// <canvas>
	fprintf(fp, "<canvas id=\"ev2d\" style=\"border: 1px solid;\" width=\"%g\" height=\"%g\">\n", ev2d_width, ev2d_height);
	fprintf(fp, "</canvas>\n");

	// date
	time_t now;
	time(&now);
	fprintf(fp, "<br>\n");
	fprintf(fp, "%s", ctime(&now));
	fprintf(fp, "<br>\n");

	// tailor
	fprintf(fp, "</body>\n");
	fprintf(fp, "</html>\n");
}


// file (option)
void ev2d_file(int type, const char fn[])
{
	ev2d_type = type;

	strcpy(ev2d_fn, fn);
}


// output
void ev2d_output(void)
{
	// open file
	FILE *fp;
	if ((fp = fopen(ev2d_fn, "w")) == NULL) {
		fprintf(stderr, "*** %s file open error.\n", ev2d_fn);
		exit(1);
	}

	// output
	if (!ev2d_type) {
		ev2d_end_html(fp, ev2d_fn);
	}
	else {
		ev2d_end_data(fp);
	}

	// close file
	fclose(fp);
}


// set color
void ev2d_setColor(unsigned char r, unsigned char g, unsigned char b)
{
	ev2d_rgb[0] = r;
	ev2d_rgb[1] = g;
	ev2d_rgb[2] = b;
}


// set color
void ev2d_setColorA(const unsigned char rgb[])
{
	memcpy(ev2d_rgb, rgb, 3 * sizeof(unsigned char));
}


// f=0-1 -> (r, g, b)
// color = 0/1 : gray/color
void ev2d_setColorV(double f, int color)
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
	ev2d_setColor(r, g, b);
}


// draw a line (primitive)
void ev2d_drawLine(double x1, double y1, double x2, double y2)
{
	// alloc
	alloc_inc();

	const int ndata = ev2d_ndata[ev2d_page];

	// index
	ev2d_data[ev2d_page][ndata].idx = 2;

	// vertex
	ev2d_data[ev2d_page][ndata].xvt[0] = x1;
	ev2d_data[ev2d_page][ndata].yvt[0] = y1;
	ev2d_data[ev2d_page][ndata].xvt[1] = x2;
	ev2d_data[ev2d_page][ndata].yvt[1] = y2;

	// color
	memcpy(ev2d_data[ev2d_page][ndata].rgb, ev2d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev2d_ndata[ev2d_page]++;
}


// draw a triangle
void ev2d_drawTriangle(double x1, double y1, double x2, double y2, double x3, double y3)
{
	ev2d_drawLine(x1, y1, x2, y2);
	ev2d_drawLine(x2, y2, x3, y3);
	ev2d_drawLine(x3, y3, x1, y1);
}


// fill a triangle (primitive)
void ev2d_fillTriangle(double x1, double y1, double x2, double y2, double x3, double y3)
{
	// alloc
	alloc_inc();

	const int ndata = ev2d_ndata[ev2d_page];

	// index
	ev2d_data[ev2d_page][ndata].idx = 3;

	// vertex
	ev2d_data[ev2d_page][ndata].xvt[0] = x1;
	ev2d_data[ev2d_page][ndata].yvt[0] = y1;
	ev2d_data[ev2d_page][ndata].xvt[1] = x2;
	ev2d_data[ev2d_page][ndata].yvt[1] = y2;
	ev2d_data[ev2d_page][ndata].xvt[2] = x3;
	ev2d_data[ev2d_page][ndata].yvt[2] = y3;

	// color
	memcpy(ev2d_data[ev2d_page][ndata].rgb, ev2d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev2d_ndata[ev2d_page]++;
}


// draw a quadangle
void ev2d_drawQuadrangle(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
	ev2d_drawLine(x1, y1, x2, y2);
	ev2d_drawLine(x2, y2, x3, y3);
	ev2d_drawLine(x3, y3, x4, y4);
	ev2d_drawLine(x4, y4, x1, y1);
}


// fill a quadangle (primitive)
void ev2d_fillQuadrangle(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
	// alloc
	alloc_inc();

	const int ndata = ev2d_ndata[ev2d_page];

	// index
	ev2d_data[ev2d_page][ndata].idx = 4;

	// vertex
	ev2d_data[ev2d_page][ndata].xvt[0] = x1;
	ev2d_data[ev2d_page][ndata].yvt[0] = y1;
	ev2d_data[ev2d_page][ndata].xvt[1] = x2;
	ev2d_data[ev2d_page][ndata].yvt[1] = y2;
	ev2d_data[ev2d_page][ndata].xvt[2] = x3;
	ev2d_data[ev2d_page][ndata].yvt[2] = y3;
	ev2d_data[ev2d_page][ndata].xvt[3] = x4;
	ev2d_data[ev2d_page][ndata].yvt[3] = y4;

	// color
	memcpy(ev2d_data[ev2d_page][ndata].rgb, ev2d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev2d_ndata[ev2d_page]++;
}


// draw a rectangle
void ev2d_drawRectangle(double x1, double y1, double x2, double y2)
{
	ev2d_drawLine(x1, y1, x2, y1);
	ev2d_drawLine(x2, y1, x2, y2);
	ev2d_drawLine(x2, y2, x1, y2);
	ev2d_drawLine(x1, y2, x1, y1);
}


// fill a rectangle
void ev2d_fillRectangle(double x1, double y1, double x2, double y2)
{
	ev2d_fillQuadrangle(x1, y1, x2, y1, x2, y2, x1, y2);
}


// ellipse
static void ellipse(double x1, double y1, double x2, double y2, unsigned char fil)
{
	// alloc
	alloc_inc();

	const int ndata = ev2d_ndata[ev2d_page];

	// index
	ev2d_data[ev2d_page][ndata].idx = fil ? 22 : 21;

	// vertex
	ev2d_data[ev2d_page][ndata].xvt[0] = x1;
	ev2d_data[ev2d_page][ndata].yvt[0] = y1;
	ev2d_data[ev2d_page][ndata].xvt[1] = x2;
	ev2d_data[ev2d_page][ndata].yvt[1] = y2;

	// color
	memcpy(ev2d_data[ev2d_page][ndata].rgb, ev2d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev2d_ndata[ev2d_page]++;
}


// draw an ellipse
void ev2d_drawEllipse(double x1, double y1, double x2, double y2)
{
	ellipse(x1, y1, x2, y2, 0);
}


// fill an ellipse
void ev2d_fillEllipse(double x1, double y1, double x2, double y2)
{
	ellipse(x1, y1, x2, y2, 1);
}


// draw a polyline (open)
void ev2d_drawPolyline(int n, const double *x, const double *y)
{
	if (n <= 0) return;

	for (int i = 0; i < n - 1; i++) {
		ev2d_drawLine(x[i], y[i], x[i + 1], y[i + 1]);
	}
}


// draw a polygon (closed)
void ev2d_drawPolygon(int n, const double *x, const double *y)
{
	if (n <= 0) return;

	ev2d_drawPolyline(n, x, y);
	ev2d_drawLine(x[n - 1], y[n - 1], x[0], y[0]);
}


// draw a string (primitive)
void ev2d_drawString(double x, double y, double h, const char str[])
{
	// alloc
	alloc_inc();

	const int ndata = ev2d_ndata[ev2d_page];

	// index
	ev2d_data[ev2d_page][ndata].idx = -3;

	// position
	ev2d_data[ev2d_page][ndata].xvt[0] = x;
	ev2d_data[ev2d_page][ndata].yvt[0] = y;

	// string
	ev2d_data[ev2d_page][ndata].str = (char *)malloc((strlen(str) + 1) * sizeof(char));
	strcpy(ev2d_data[ev2d_page][ndata].str, str);

	// height
	ev2d_data[ev2d_page][ndata].hgt = h;

	// color
	memcpy(ev2d_data[ev2d_page][ndata].rgb, ev2d_rgb, 3 * sizeof(unsigned char));

	// data++
	ev2d_ndata[ev2d_page]++;
}
