#ifndef _ORCWA_PROTOTYPE_H_
#define _ORCWA_PROTOTYPE_H_

#ifdef __CUDACC__
#include <cuComplex.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// C
extern void        alloc_farfield(void);
extern d_complex_t calcdft(int, const double [], double, double, double);
extern void        calcFar1d(const double [], d_complex_t ***, d_complex_t ***);
extern void        calcFar2d(const double [], d_complex_t ***, d_complex_t ***);
extern void        calcNear1d(void);
extern void        calcNear2d(void);
extern void        calcPin(void);
extern void        calcSpara(void);
extern void        calcZin(void);
extern d_complex_t coupling(int, int, int);
extern double      cputime(void);
extern double      factorMur(double, id_t);
extern void        farComponent(d_complex_t, d_complex_t, double []);
extern double      farfactor(int);
extern void        farfield(int, double, double, double, d_complex_t *, d_complex_t *);
extern void        fitgeometry(void);
extern int         geomlines(int, int, int *, int *, double (*)[8], double (*)[2][3], int *, double);
extern void        getspan(const double [], int, int, int, double, double, int *, int *, double);
extern int         ingeometry(double, double, double, int, double *, double);
extern int         inout3(double, double, double [3][2], double);
extern int         input_data(FILE *);
extern void        memalloc1(void);
extern void        memalloc2(void);
extern void        memalloc3(void);
extern void        memfree1(void);
extern void        memfree2(void);
extern void        memfree3(void);
extern void        monitor1(FILE *, const char []);
extern void        monitor2(FILE *, int, int);
extern void        monitor3(FILE *, const char [], const char []);
extern void        monitor4(FILE *, const double []);
extern int         nearest(double, int, int, const double *);
extern void        NodeE_c(int, int, int, int, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        NodeH_c(int, int, int, int, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        outputChars(FILE *);
extern void        outputCoupling(FILE *);
extern void        outputCross(FILE *);
extern void        outputFar0d(void);
extern void        outputFar1d(void);
extern void        outputFar2d(void);
extern void        outputNear1d(void);
extern void        outputNear2d(void);
extern void        outputSpara(FILE *);
extern void        outputZin(FILE *);
extern void        planewave(double, double, double, double, d_complex_t [], d_complex_t []);
extern void        plot2dCoupling(int, int, int, double ***, scale_t, int, const double [], const char [], double, double, double);
extern void        plot2dFar0d0(int, double *[7], int, double, double, int, int, double, double, const char [], const char [], const char [], double, double, double);
extern void        plot2dFar1d0(int, double (*)[7], int [], char, double, int, int, int, double, double, int, const char [], const char [], double, double, double, double);
extern void        plot2dFchar(int, const double [], scale_t, int, const double [], const char [], const char [], const char [], int, int, int);
extern void        plot2dFeed(void);
extern void        plot2dFreq(void);
extern void        plot2dIter(void);
extern void        plot2dNear1d0(const char [], int, double *[3], double *[3], const double [], int, int, double, double, int, const char [], double, double [2][3], double, double, double);
extern void        plot2dNear2d0(int, int, double **, const double [], const double [], int, int, int, double, double, int, const char [], double, char, double, const char [], int, double (*)[2][3], double, double, double);
extern void        plot2dPoint(void);
extern void        plot2dRef(int, int, const d_complex_t [], const double [], scale_t, int, const double [], const char [], double, double, double);
extern void        plot2dSmith(int, int, const d_complex_t [], const double [], const double [], const char [], double, double, double);
extern void        plot2dSpara(int, int, const d_complex_t [], scale_t, int, const double [], const char [], double, double, double);
extern void        plot2dYin(int, int, const d_complex_t [], scale_t, int, const double [], const char [], double, double, double);
extern void        plot2dZin(int, int, const d_complex_t [], scale_t, int, const double [], const char [], double, double, double);
extern void        plot3dFar2d(d_complex_t ***, d_complex_t ***);
extern void        plot3dFar2d0(int, int, double **, int, int, double, double, int, double (*)[2][3], double, int, char **);
extern void        plot3dGeom(void);
extern void        plot3dNear2d0(int, int, double **, const double [], const double [], int, int, int, double, double, int, const char [], double, char, double, const char [], int, double (*)[2][3]);
extern void        post(void);
extern int         post_data(FILE *);
// post_data() が設定する: RCWA モードの入力 (rcwalayer を含む) なら 1
extern int         PostRcwaInput;
extern void        readout(FILE *);
extern void        rectangleContour(double [4][3], double, double, int);
extern void        rgbColor(double, unsigned char *, unsigned char *, unsigned char *, int);
extern void        setup(void);
extern void        setupABCsize(void);
extern void        setupDft(void);
extern void        setupDispersion(void);
extern void        setupId(void);
extern void        setupMurHx(int);
extern void        setupMurHy(int);
extern void        setupMurHz(int);
extern void        setupPml(void);
extern void        setupPmlEx(int);
extern void        setupPmlEy(int);
extern void        setupPmlEz(int);
extern void        setupPmlHx(int);
extern void        setupPmlHy(int);
extern void        setupPmlHz(int);
extern void        setupSize(int, int, int, int);
extern void        setup_cells(int, int, int, int *, int *, int *);
extern void        setup_center(void);
extern void        setup_geomline3d(void);
extern void        setup_farfield(void);
extern void        setup_feed(const double *, const double *, const double *);
extern void        setup_load(int, char *, double *, double *, double *, char *, double *, int);
extern void        setup_node(int, int, int, double *, double *, double *, int *, int *, int *);
extern void        setup_planewave(void);
extern void        setup_point(const double *, const double *, const double *, const char []);
extern int         tokenize(char *, const char *, char *[], int);
extern double      vfeed(double, double, double);
extern void        writeout(FILE *);

// MPI
extern void        setup_mpi(void);
extern void        comm_average(double []);
extern void        comm_broadcast(void);
extern void        comm_check(int, int, int);
extern double      comm_cputime(void);
extern void        comm_feed(void);
extern int         comm_inproc(int, int, int);
extern void        comm_near3d(void);
extern void        comm_point(void);
extern void        comm_X(int);
extern void        comm_Y(int);
extern void        comm_Z(int);
extern void        mpi_init(int, char **);
extern void        mpi_close(void);

/* 発熱量密度 (sol/powerloss.c)。CUDA 版は無く C の実装を全ビルドで共有する
   ので、C++ (.cu) から呼べるよう extern "C" の中に置く。 */
extern void        calculatePowerLoss(double *);

#ifdef __cplusplus
}
#endif

// C + CUDA
extern void        average(double []);
extern void        dftNear1d(int, int *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear1dX(int, int, int64_t, int64_t, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear1dY(int, int, int64_t, int64_t, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear1dZ(int, int, int64_t, int64_t, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear2d(int, int *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear2dX(int, int64_t, int64_t, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear2dY(int, int64_t, int64_t, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear2dZ(int, int64_t, int64_t, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *, d_complex_t *);
extern void        dftNear3d(int);
extern void        efeed(int);
extern void        eload(void);
extern void        initfield(void);
extern void        murH(int64_t, mur_t *, real_t *);
extern void        pmlEx(void);
extern void        pmlEy(void);
extern void        pmlEz(void);
extern void        pmlHx(void);
extern void        pmlHy(void);
extern void        pmlHz(void);
extern void        pbcx(void);
extern void        pbcy(void);
extern void        pbcz(void);
extern void        solve(int, double *, FILE *);
extern void        updateEx(double);
extern void        updateEy(double);
extern void        updateEz(double);
extern void        updateHx(double);
extern void        updateHy(double);
extern void        updateHz(double);
extern void        dispersionEx(double);
extern void        dispersionEy(double);
extern void        dispersionEz(double);
extern void        vpoint(int);

// CUDA
#ifdef __CUDACC__
extern void        comm_cuda_X(int);
extern void        comm_cuda_Y(int);
extern void        comm_cuda_Z(int);
extern void        cuda_free(int, void *);
extern void        cuda_malloc(int, int, void **, size_t);
extern void        cuda_memcpy(int, void *, const void *, size_t, cudaMemcpyKind);
extern void        cuda_memset(int, void *, int, size_t);
extern void        info_gpu(FILE *, int, int, int);
extern void        info_gpu_mpi(FILE *, int, const int [], int, int, int, int, int);
extern void        memalloc2_gpu();
extern void        memalloc3_gpu();
extern void        memcopy3_gpu();
extern void        memfree2_gpu();
extern void        memfree3_gpu();
extern void        setup_gpu();
extern void        setup_host();
#endif

#endif
