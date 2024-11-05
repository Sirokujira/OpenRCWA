// ev.h
#ifndef _EV_H_
#define _EV_H_

#ifdef __cplusplus
extern "C" {
#endif

// ev2d
extern void ev2d_init(double, double);
extern void ev2d_newPage(void);
extern void ev2d_file(int, const char []);
extern void ev2d_output(void);
extern void ev2d_setColor(unsigned char, unsigned char, unsigned char);
extern void ev2d_setColorA(const unsigned char []);
extern void ev2d_setColorV(double, int);
extern void ev2d_drawLine(double, double, double, double);
extern void ev2d_drawTriangle(double, double, double, double, double, double);
extern void ev2d_fillTriangle(double, double, double, double, double, double);
extern void ev2d_drawQuadrangle(double, double, double, double, double, double, double, double);
extern void ev2d_fillQuadrangle(double, double, double, double, double, double, double, double);
extern void ev2d_drawRectangle(double, double, double, double);
extern void ev2d_fillRectangle(double, double, double, double);
extern void ev2d_drawPolyline(int, const double *, const double *);
extern void ev2d_drawPolygon(int, const double *, const double *);
extern void ev2d_drawEllipse(double, double, double, double);
extern void ev2d_fillEllipse(double, double, double, double);
extern void ev2d_drawString(double, double, double, const char []);

// ev3d
extern void ev3d_init(void);
extern void ev3d_newPage(void);
extern void ev3d_file(int, const char [], int);
extern void ev3d_output(void);
extern void ev3d_setColor(unsigned char, unsigned char, unsigned char);
extern void ev3d_setColorA(const unsigned char []);
extern void ev3d_setColorV(double, int);
extern void ev3d_drawLine(double, double, double, double, double, double);
extern void ev3d_drawTriangle(double, double, double, double, double, double, double, double, double);
extern void ev3d_fillTriangle(double, double, double, double, double, double, double, double, double);
extern void ev3d_drawQuadrangle(double, double, double, double, double, double, double, double, double, double, double, double);
extern void ev3d_fillQuadrangle(double, double, double, double, double, double, double, double, double, double, double, double);
extern void ev3d_drawRectangle(char, double, double, double, double, double);
extern void ev3d_fillRectangle(char, double, double, double, double, double);
extern void ev3d_drawBox(double, double, double, double, double, double);
extern void ev3d_fillBox(double, double, double, double, double, double);
extern void ev3d_drawEllipse(char, double, double, double, double, double, int);
extern void ev3d_fillEllipse(char, double, double, double, double, double, int);
extern void ev3d_drawTitle(const char []);
extern void ev3d_html_size(int, int);
extern void ev3d_html_angle(double, double);
extern void ev3d_html_font(int, int);
extern void ev3d_index(int);

// ev2dlib
extern void ev2dlib_grid(double, double, double, double, int, int);
extern void ev2dlib_func2(int, const double [], const double [], double, double, double, double, double, double);
extern void ev2dlib_func1(int, const double [], double, double, double, double, double, double);
extern void ev2dlib_sample(double, double, double, double, int);
extern void ev2dlib_contour(int, int, const double [], const double [], double **, double, double, int);
extern void ev2dlib_arrow(double, double, double, double);
extern void ev2dlib_Xaxis(double, double, double, double, const char *, const char *, const char *);
extern void ev2dlib_Yaxis(double, double, double, double, const char *, const char *, const char *);
extern void ev2dlib_CircleFunc(double, double, double, double, const double [], int, double, double, double, double);
extern void ev2dlib_CircleMesh(double, double, double, double, int, int);
extern void ev2dlib_drawLineA(double [2][2]);
extern void ev2dlib_drawTriangleA(double [3][2]);
extern void ev2dlib_fillTriangleA(double [3][2]);
extern void ev2dlib_drawQuadrangleA(double [4][2]);
extern void ev2dlib_fillQuadrangleA(double [4][2]);

// ev3dlib
extern void ev3dlib_axis(double, double, double, double, int);
extern void ev3dlib_func(int, double, double, int, double, double, double **);
extern void ev3dlib_drawLineA(double [2][3]);
extern void ev3dlib_drawTriangleA(double [3][3]);
extern void ev3dlib_fillTriangleA(double [3][3]);
extern void ev3dlib_drawQuadrangleA(double [4][3]);
extern void ev3dlib_fillQuadrangleA(double [4][3]);

#ifdef __cplusplus
}
#endif

#endif		// _EV_H_
