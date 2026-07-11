#pragma once
#include "defns.h"
#include <Eigen/Core>

#if !defined(USE_MKL) && !defined(USE_LAPACKE)
// MKL/LAPACKE のいずれも利用できない環境では Eigen の固有値ソルバにフォールバックする
#include <Eigen/Eigenvalues>

template <class MatrixType>
class MKLEigenSolver
{
private:
    typedef typename MatrixType::Scalar ScalarType;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorType;

    Eigen::ComplexEigenSolver<MatrixType> es_;

public:
    MKLEigenSolver() {}
    void compute(const MatrixType& A) { es_.compute(A, true); }
    const MatrixType& eigenvectors() const { return es_.eigenvectors(); }
    const VectorType& eigenvalues() const { return es_.eigenvalues(); }
};

#elif defined(USE_LAPACKE)
// MKL が無いがシステムの LAPACKE (zgeev/cgeev) が利用できる環境向けの実装。
// std::complex<double>/<float> を lapack_complex_double/float としてそのまま
// LAPACKE に渡せるようにし、Eigen::ComplexEigenSolver より高速な経路を使う。
// zgeev の失敗 (info != 0) や残差の大きい結果は Eigen ソルバで解き直す。
#include <complex>
#include <iostream>
#include <vector>
#include <Eigen/Eigenvalues>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
extern "C" {
#include <lapacke.h>
}

template <class MatrixType>
class MKLEigenSolver
{
private:
    typedef typename MatrixType::Scalar ScalarType;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorType;

    MatrixType V_;
    VectorType d_;

    // ‖A・V − V・diag(d)‖∞ / ‖A‖∞ による検算
    bool verified(const MatrixType& A) const
    {
        const auto normA = A.cwiseAbs().maxCoeff();
        if (!(normA > 0)) return true;
        const auto res =
            (A * V_ - V_ * d_.asDiagonal()).cwiseAbs().maxCoeff();
        return res <= 1e-8 * normA;
    }

    void computeEigenFallback(const MatrixType& A)
    {
        Eigen::ComplexEigenSolver<MatrixType> es(A, true);
        V_ = es.eigenvectors();
        d_ = es.eigenvalues();
    }

public:
    MKLEigenSolver() {}
    void compute(const MatrixType& A);
    const MatrixType& eigenvectors() const { return V_; }
    const VectorType& eigenvalues() const { return d_; }
};

template <class MatrixType>
void MKLEigenSolver<MatrixType>::compute(const MatrixType& A)
{
    int n = A.rows();
    lapack_int info = 0;
    if constexpr (std::is_same_v<ScalarType, std::complex<double>>) {
        std::vector<std::complex<double>> a(n * n), w(n), vr(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[i * n + j] = A(i, j);

        info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n,
            a.data(), n, w.data(), nullptr, n, vr.data(), n);

        V_.resize(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                V_(i, j) = vr[i * n + j];

        d_.resize(n);
        for (int i = 0; i < n; ++i)
            d_(i) = w[i];
    } else if constexpr (std::is_same_v<ScalarType, std::complex<float>>) {
        std::vector<std::complex<float>> a(n * n), w(n), vr(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[i * n + j] = A(i, j);

        info = LAPACKE_cgeev(LAPACK_ROW_MAJOR, 'N', 'V', n,
            a.data(), n, w.data(), nullptr, n, vr.data(), n);

        V_.resize(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                V_(i, j) = vr[i * n + j];

        d_.resize(n);
        for (int i = 0; i < n; ++i)
            d_(i) = w[i];
    } else {
        ScalarType::unimplemented;
    }

    if (info != 0 || !verified(A)) {
        std::clog << "[MKLEigenSolver] LAPACKE geev unreliable (info="
                  << info << "), falling back to Eigen solver\n";
        computeEigenFallback(A);
        if (!verified(A)) {
            std::clog << "[MKLEigenSolver] warning: Eigen solver residual"
                         " also exceeds tolerance\n";
        }
    }
}

#else // USE_MKL
#include <mkl.h>

template <class MatrixType>
class MKLEigenSolver
{
private:
    typedef typename MatrixType::Scalar ScalarType;
    typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorType;

    MatrixType V_;
	VectorType d_;

public:
	MKLEigenSolver() {}
	void compute(const MatrixType& A);
	const MatrixType& eigenvectors() const { return V_; }
	const VectorType& eigenvalues() const { return d_; }
};


template <class MatrixType>
void MKLEigenSolver<MatrixType>::compute(const MatrixType& A)
{
    int n = A.rows();
    if constexpr (std::is_same_v<ScalarType, std::complex<double>>) {
        lapack_complex_double * a = nullptr;
        lapack_complex_double * w = nullptr;
        lapack_complex_double * vr = nullptr;

        a = new lapack_complex_double[n*n];
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                a[i * n + j].real = A(i, j).real();
                a[i * n + j].imag = A(i, j).imag();
            }
        }
        w = new lapack_complex_double[n];
        vr = new lapack_complex_double[n*n];

        LAPACKE_zgeev(LAPACK_ROW_MAJOR,
            'N',
            'V',
            n,
            a,
            n,
            w,
            nullptr,
            n,
            vr,
            n);

        V_.resize(n , n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                V_(i, j) = scalex(vr[i * n + j].real, vr[i * n + j].imag);
            }
        }

        d_.resize(n);
        for (int i = 0; i < n; ++i)
        {
            d_(i) = scalex(w[i].real, w[i].imag);
        }

        delete[] a;
        delete[] w;
        delete[] vr;
    } else if constexpr (std::is_same_v<ScalarType, std::complex<float>>) {
        
        lapack_complex_float * a = nullptr;
        lapack_complex_float * w = nullptr;
        lapack_complex_float * vr = nullptr;

        a = new lapack_complex_float[n*n];
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                a[i * n + j].real = A(i, j).real();
                a[i * n + j].imag = A(i, j).imag();
            }
        }
        w = new lapack_complex_float[n];
        vr = new lapack_complex_float[n*n];

        LAPACKE_cgeev(LAPACK_ROW_MAJOR,
            'N',
            'V',
            n,
            a,
            n,
            w,
            nullptr,
            n,
            vr,
            n);

        V_.resize(n , n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                V_(i, j) = scalex(vr[i * n + j].real, vr[i * n + j].imag);
            }
        }

        d_.resize(n);
        for (int i = 0; i < n; ++i)
        {
            d_(i) = scalex(w[i].real, w[i].imag);
        }

        delete[] a;
        delete[] w;
        delete[] vr;
    } else {
        ScalarType::unimplemented;
    }
}
#endif // USE_MKL