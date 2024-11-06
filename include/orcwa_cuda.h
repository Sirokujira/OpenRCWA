// orcwa_cuda.h
#ifndef _ORCWA_CUDA_H_
#define _ORCWA_CUDA_H_

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

#define LA(p,i,j,k) ((i)*((p)->Ni)+(j)*((p)->Nj)+(k)*((p)->Nk)+((p)->N0))
#define CEIL(n,d) (((n) + (d) - 1) / (d))

// parameter
typedef struct {
	int64_t Ni, Nj, Nk, N0;
	int     Nx, Ny, Nz;
	int     iMin, iMax, jMin, jMax, kMin, kMax;
	int     NFeed, IPlanewave;
	real_t  ei[3], hi[3], r0[3], ri[3], ai, dt;
} param_t;

EXTERN int          GPU;
EXTERN int          UM;

// constant memory
EXTERN param_t h_Param;
__constant__ param_t d_Param;
__constant__ bid_t d_Bid;

// execution configuration
EXTERN dim3         updateBlock;
EXTERN int          dispersionBlock;
EXTERN dim3         sumGrid, sumBlock;
EXTERN int          murBlock;
EXTERN int          pmlBlock;
EXTERN int          pbcBlock;
EXTERN int          near1dBlock;
EXTERN dim3         near2dBlock;
EXTERN dim3         bufBlock;

// host memory
EXTERN real_t       *h_Xn, *h_Yn, *h_Zn;
EXTERN real_t       *h_Xc, *h_Yc, *h_Zc;

// device memory
EXTERN real_t       *d_Xn, *d_Yn, *d_Zn;
EXTERN real_t       *d_Xc, *d_Yc, *d_Zc;
EXTERN real_t       *d_RXn, *d_RYn, *d_RZn;
EXTERN real_t       *d_RXc, *d_RYc, *d_RZc;
EXTERN id_t         *d_iEx, *d_iEy, *d_iEz;
EXTERN id_t         *d_iHx, *d_iHy, *d_iHz;
EXTERN real_t       *d_C1, *d_C2;
EXTERN real_t       *d_D1, *d_D2;
EXTERN real_t       *d_K1Ex, *d_K2Ex, *d_K1Ey, *d_K2Ey, *d_K1Ez, *d_K2Ez;
EXTERN real_t       *d_K1Hx, *d_K2Hx, *d_K1Hy, *d_K2Hy, *d_K1Hz, *d_K2Hz;
EXTERN real_t       *d_DispersionEx, *d_DispersionEy, *d_DispersionEz;
EXTERN dispersion_t *d_mDispersionEx, *d_mDispersionEy, *d_mDispersionEz;
EXTERN mur_t        *d_fMurHx, *d_fMurHy, *d_fMurHz;
EXTERN pml_t        *d_fPmlEx, *d_fPmlEy, *d_fPmlEz;
EXTERN pml_t        *d_fPmlHx, *d_fPmlHy, *d_fPmlHz;
EXTERN real_t       *d_gPmlXn, *d_gPmlYn, *d_gPmlZn;
EXTERN real_t       *d_gPmlXc, *d_gPmlYc, *d_gPmlZc;
EXTERN real_t       *d_rPmlE, *d_rPmlH, *d_rPml;
EXTERN feed_t       *d_Feed;
EXTERN double       *d_VFeed, *d_IFeed;
EXTERN inductor_t   *d_Inductor;
EXTERN point_t      *d_Point;
EXTERN double       *d_VPoint;
EXTERN real_t       *h_sumE, *h_sumH;
EXTERN real_t       *d_sumE, *d_sumH;
EXTERN float        *d_cEx_r, *d_cEx_i, *d_cEy_r, *d_cEy_i, *d_cEz_r, *d_cEz_i;
EXTERN float        *d_cHx_r, *d_cHx_i, *d_cHy_r, *d_cHy_i, *d_cHz_r, *d_cHz_i;
EXTERN real_t       *d_Sendbuf_hy_x, *d_Sendbuf_hz_x, *d_Sendbuf_hz_y, *d_Sendbuf_hx_y, *d_Sendbuf_hx_z, *d_Sendbuf_hy_z;
EXTERN real_t       *d_Recvbuf_hy_x, *d_Recvbuf_hz_x, *d_Recvbuf_hz_y, *d_Recvbuf_hx_y, *d_Recvbuf_hx_z, *d_Recvbuf_hy_z;

#endif  // _ORCWA_CUDA_H_
