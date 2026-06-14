#pragma once

#if defined(__CUDACC__) || __has_include(<cuComplex.h>)
#include <cuComplex.h>
#define ACACIA_HAS_CUDA 1
#endif

#if !defined(NAMESPACE_BEGIN)
#   define NAMESPACE_BEGIN(name) namespace name {
#endif
#if !defined(NAMESPACE_END)
#   define NAMESPACE_END(name) }
#endif

// #define ACA_USE_FLOAT32
NAMESPACE_BEGIN(acacia::gpu)
#ifdef ACACIA_HAS_CUDA
#define HAVE_CUBLAS
#ifdef ACA_USE_FLOAT32
  using complex_t = cuFloatComplex;
  using Real = float;
#else
  using complex_t = cuDoubleComplex;
  using Real = double;
#endif

#ifdef ACA_USE_FLOAT32
  #define acacia_gpu_em_Complexgetrf_ cusolverDnCgetrf
  #define acacia_gpu_em_Complexgetrs_ cusolverDnCgetrs
  #define acacia_gpu_em_Complexgetrf_bufferSize cusolverDnCgetrf_bufferSize
  #define acacia_gpu_em_Complexgemm3m_ cublasCgemm3m
  #define acacia_gpu_em_Complexaxpy_ cublasCaxpy
  #define acacia_gpu_em_Complexscal_ cublasCscal
  #define acacia_gpu_em_Complexamax_ cublasIcamax
#else
  #define acacia_gpu_em_Complexgetrf_ cusolverDnZgetrf
  #define acacia_gpu_em_Complexgetrs_ cusolverDnZgetrs
  #define acacia_gpu_em_Complexgemm3m_ cublasZgemm3m
  #define acacia_gpu_em_Complexgetrf_bufferSize cusolverDnZgetrf_bufferSize
  #define acacia_gpu_em_Complexaxpy_ cublasZaxpy
  #define acacia_gpu_em_Complexscal_ cublasZscal
  #define acacia_gpu_em_Complexamax_ cublasIzamax
#endif
#else // !ACACIA_HAS_CUDA
// CUDA ツールキットがない CPU ビルド用に cuComplex とレイアウト互換の型を定義する
struct acacia_float2  { float  x, y; };
struct acacia_double2 { double x, y; };
#ifdef ACA_USE_FLOAT32
  using complex_t = acacia_float2;
  using Real = float;
#else
  using complex_t = acacia_double2;
  using Real = double;
#endif
#endif // ACACIA_HAS_CUDA
NAMESPACE_END(acacia::gpu)
