#pragma once
#include "defns.h"
#include <Eigen/Core>
#include <atomic>

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
    void compute(const MatrixType& A)
    {
        // RCWA の波動方程式行列は ±m 次数の構造縮退を持ち、縮退固有値の
        // 固有ベクトル基底が処理系依存でほぼ線形従属になることがある
        // (MSVC ビルドで回折効率が発散する事象を確認)。対角に微小な
        // 非一様摂動を加えて縮退を決定的に分離し、良条件の基底を得る。
        // 固有値の変化は ~1e-11・‖A‖ で回折効率への影響は無視できる。
        MatrixType Ap = A;
        const double scale = (double)A.cwiseAbs().maxCoeff();
        if (scale > 0) {
            for (Eigen::Index i = 0; i < Ap.rows(); ++i)
                Ap(i, i) += ScalarType(scale * 1e-11 * (double)(i + 1));
        }
        es_.compute(Ap, true);
    }
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

    // ‖A・V − V・diag(d)‖∞ / ‖A‖∞ (相対残差)
    double relResidual(const MatrixType& A) const
    {
        const auto normA = A.cwiseAbs().maxCoeff();
        if (!(normA > 0)) return 0.0;
        const auto res =
            (A * V_ - V_ * d_.asDiagonal()).cwiseAbs().maxCoeff();
        return static_cast<double>(res) / static_cast<double>(normA);
    }

    // 相対残差による検算
    bool verified(const MatrixType& A) const
    {
        return relResidual(A) <= 1e-8;
    }

    void computeEigenFallback(const MatrixType& A)
    {
        // 縮退固有値の基底不良条件対策の微小対角摂動 (Eigen フォールバック側と同じ)
        MatrixType Ap = A;
        const double scale = (double)A.cwiseAbs().maxCoeff();
        if (scale > 0) {
            for (Eigen::Index i = 0; i < Ap.rows(); ++i)
                Ap(i, i) += ScalarType(scale * 1e-11 * (double)(i + 1));
        }
        Eigen::ComplexEigenSolver<MatrixType> es(Ap, true);
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
    // 検算に落ちたら以後のソルブでは LAPACKE を試さず直接 Eigen で解く
    // (プロセス内で粘着)。macOS/arm64 の Homebrew LAPACKE で全ソルブが
    // 相対残差 O(1) (実測 1.96〜3.60) の壊れた固有ベクトルを返す事象があり、
    // 毎回「壊れた zgeev + 検算 + Eigen 解き直し」を払うのは無駄なため。
    // Eigen 経路は縮退対策の対角摂動つきで Windows が常用しており、
    // どのソルブに対しても正しい結果を返すので、粘着させても
    // 正しさは損なわれない (性能の判断だけを固定する)。
    static std::atomic<bool> lapackeBroken{false};
    if (lapackeBroken.load(std::memory_order_relaxed)) {
        computeEigenFallback(A);
        return;
    }

    int n = A.rows();
    lapack_int info = 0;
    // LAPACK_COL_MAJOR で呼ぶ: LAPACKE の row-major ラッパは内部で
    // 行列全体の転置コピーとワーク領域を作る。コピーはどのみち自前で
    // 行っているので、詰め替え時に転置して直接 col-major で渡す方が
    // 余分な転置が消え、ラッパ層の実装差 (プラットフォーム依存の疑い)
    // も踏まない。
    if constexpr (std::is_same_v<ScalarType, std::complex<double>>) {
        std::vector<std::complex<double>> a(n * n), w(n), vr(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[j * n + i] = A(i, j);   // col-major 詰め替え

        info = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'V', n,
            a.data(), n, w.data(), nullptr, n, vr.data(), n);

        V_.resize(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                V_(i, j) = vr[j * n + i];   // col-major 読み出し

        d_.resize(n);
        for (int i = 0; i < n; ++i)
            d_(i) = w[i];
    } else if constexpr (std::is_same_v<ScalarType, std::complex<float>>) {
        std::vector<std::complex<float>> a(n * n), w(n), vr(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[j * n + i] = A(i, j);   // col-major 詰め替え

        info = LAPACKE_cgeev(LAPACK_COL_MAJOR, 'N', 'V', n,
            a.data(), n, w.data(), nullptr, n, vr.data(), n);

        V_.resize(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                V_(i, j) = vr[j * n + i];   // col-major 読み出し

        d_.resize(n);
        for (int i = 0; i < n; ++i)
            d_(i) = w[i];
    } else {
        ScalarType::unimplemented;
    }

    if (info != 0 || !verified(A)) {
        // 相対残差を併記する: 「わずかに閾値 (1e-8) を超えた」のか
        // 「桁違いに壊れている」のかをログだけで切り分けられるようにする。
        // 実測では macOS/arm64 で O(1) (= 出力が完全に壊れている) だった。
        std::clog << "[MKLEigenSolver] LAPACKE geev unreliable (info="
                  << info << ", relative residual=" << relResidual(A)
                  << " > 1e-8); using the Eigen solver for the rest of"
                     " this process\n";
        lapackeBroken.store(true, std::memory_order_relaxed);
        computeEigenFallback(A);
        if (!verified(A)) {
            std::clog << "[MKLEigenSolver] warning: Eigen solver residual"
                         " also exceeds tolerance (relative residual="
                      << relResidual(A) << ")\n";
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