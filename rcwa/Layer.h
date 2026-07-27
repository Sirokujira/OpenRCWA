#pragma once
#include "defns.h"
#include <Eigen/Core>
#include <vector>
#include "Parameterization.h"


class Layer
{
private:
	Eigen::VectorXs coordX_;
	Eigen::VectorXs coordY_;
	Eigen::MatrixXcs eps_;

	int nCx_;
	int nCy_;


	bool isSolved_;
	Eigen::MatrixXcs eigvecE_;
	Eigen::MatrixXcs eigvecH_;
	Eigen::VectorXcs gamma_;

	scalar lambda_ = 0;  // 直近の solve() で用いた波長 (k0 の復元用)

public:
	Layer(const Eigen::VectorXs& coordX,
		const Eigen::VectorXs& coordY,
		const Eigen::MatrixXcs& eps);

	const Eigen::MatrixXcs& eps() const { return eps_; }
	scalex eps(scalar x, scalar y) const;

	// solve() 時に設定された波長と対応する真空波数 k0 = 2π/λ。
	scalar lambda() const { return lambda_; }
	scalar k0() const { return 2 * Pi / lambda_; }

	// 誘電率の畳み込み行列 [[ε]] (ContinuousXY 規則, サイズ nx*ny × nx*ny)。
	// 縦方向成分 Ez = i·[[ε]]⁻¹(kx·Hy − ky·Hx) の評価に用いる。
	Eigen::MatrixXcs epsConvolution(int nx, int ny) const;
	
	scalar Lx() const { return coordX_[nCx_-1] - coordX_[0]; }
	scalar Ly() const { return coordY_[nCy_-1] - coordY_[0]; }

	scalar b0() const { return coordY_[0]; }
	scalar b1() const { return coordY_[1]; }
	scalar u0() const { return coordY_[nCy_-2]; }
	scalar u1() const { return coordY_[nCy_-1]; }

	scalar l0() const { return coordX_[0]; }
	scalar l1() const { return coordX_[1]; }
	scalar r0() const { return coordX_[nCx_-2]; }
	scalar r1() const { return coordX_[nCx_-1]; }

	scalar coordX(int layout) const { return coordX_[layout]; }
	scalar coordY(int layout) const { return coordY_[layout]; }
	const Eigen::VectorXs& coordX() const { return coordX_; }
	const Eigen::VectorXs& coordY() const { return coordY_; }
	void setCoordX(const Eigen::VectorXs& coordX) { coordX_ = coordX; isSolved_ = false; }
	void setCoordY(const Eigen::VectorXs& coordY) { coordY_ = coordY; isSolved_ = false; }

	void setCoordX(int layout, scalar coordX) { coordX_[layout] = coordX; isSolved_ = false; }
	void setCoordY(int layout, scalar coordY) { coordY_[layout] = coordY; isSolved_ = false; }
	
	void setEps(int row, const Eigen::VectorXcs& reps) { eps_.row(row) = reps; isSolved_ = false; }
	int nCx() const { return nCx_; }
	int nCy() const { return nCy_; }

	bool isSolved() const { return isSolved_; }

    void waveEqnCoeff(const Eigen::MatrixXcs& Kx, 
        const Eigen::MatrixXcs& Ky, 
        scalar lambda, 
        int nx, 
        int ny,
        Eigen::MatrixXcs& P,
        Eigen::MatrixXcs& Q,
        std::vector<Eigen::MatrixXcs>* dP = nullptr,
        std::vector<Eigen::MatrixXcs>* dQ = nullptr,
        const Parameterization& para_x = {},
        const Parameterization& para_y = {},
        Eigen::MatrixXcs * dPl = nullptr,
        Eigen::MatrixXcs * dQl = nullptr);
        
	void solve(const Eigen::MatrixXcs& Kx, 
        const Eigen::MatrixXcs& Ky, 
        scalar lambda, 
        int nx, 
        int ny);

	void permuteEigVecX(
		Eigen::MatrixXcs& eigvecE,
		Eigen::MatrixXcs& eigvecH,
		scalar dx, 
		int nx, 
		int ny);

	void dpermuteEigVecX(
		Eigen::MatrixXcs& deigvecE,
		Eigen::MatrixXcs& deigvecH,
		scalar dx,
		int nx, 
		int ny);
    
	const Eigen::VectorXcs& gamma() const { return gamma_; }
	const Eigen::MatrixXcs& eigvecE() const { return eigvecE_; }
	const Eigen::MatrixXcs& eigvecH() const { return eigvecH_; }

	void setSolutions(
		const Eigen::MatrixXcs& eigvecE,
		const Eigen::MatrixXcs& eigvecH,
		const Eigen::VectorXcs& gamma);
	
	// px/py は複素振幅 (円偏波・任意位相の偏波に対応)。scalar からは暗黙変換される。
	void generatePlaneWave(
		const Eigen::VectorXcs& delta,
		scalex px,
		scalex py,
		Eigen::VectorXcs& c);
	
	void getHarmonics(
		const Eigen::VectorXcs& c,
		Eigen::VectorXcs& harmonics);
};
