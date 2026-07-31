#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Layer.h"
#include "FourierSolver2D.h"
#include "WaveEquationSolver.h"

#include <vector>
#include <algorithm>
#include <iomanip>
#include "utils/timer.hpp"

using namespace Eigen;
using namespace std;

Layer::Layer(const Eigen::VectorXs& coordX,
	const Eigen::VectorXs& coordY,
	const Eigen::MatrixXcs& eps):
	coordX_(coordX),
	coordY_(coordY), 
	eps_(eps),
	isSolved_(false)
{
	nCx_ = coordX.size();
	nCy_ = coordY.size();

	if (eps.rows() != nCy_ - 1 ||
		eps.cols() != nCx_ - 1)
	{
		std::cerr << "The configuration of permittivity does not match coordinates!\n";
	}

	// 非磁性 (μr = 1) が既定
	mu_ = Eigen::MatrixXcs::Ones(eps_.rows(), eps_.cols());
	isMagnetic_ = false;
}

Layer::Layer(const Eigen::VectorXs& coordX,
	const Eigen::VectorXs& coordY,
	const Eigen::MatrixXcs& eps,
	const Eigen::MatrixXcs& mu):
	coordX_(coordX),
	coordY_(coordY),
	eps_(eps),
	mu_(mu),
	isSolved_(false)
{
	nCx_ = coordX.size();
	nCy_ = coordY.size();

	if (eps.rows() != nCy_ - 1 ||
		eps.cols() != nCx_ - 1)
	{
		std::cerr << "The configuration of permittivity does not match coordinates!\n";
	}
	if (mu.rows() != eps.rows() || mu.cols() != eps.cols())
	{
		std::cerr << "The configuration of permeability does not match permittivity!\n";
		mu_ = Eigen::MatrixXcs::Ones(eps_.rows(), eps_.cols());
	}

	// μ が全セル 1 なら非磁性として扱い、P/Q の高速経路を使う。
	isMagnetic_ = false;
	for (int i = 0; i < mu_.rows() && !isMagnetic_; ++i)
		for (int j = 0; j < mu_.cols(); ++j)
			if (std::abs(mu_(i, j) - scalex(1.0, 0.0)) > 1e-12)
			{
				isMagnetic_ = true;
				break;
			}
}

void Layer::waveEqnCoeff(const Eigen::MatrixXcs& Kx, 
    const Eigen::MatrixXcs& Ky, 
    scalar lambda, 
    int nx, 
    int ny,
    Eigen::MatrixXcs& P,
    Eigen::MatrixXcs& Q,
    std::vector<Eigen::MatrixXcs>* dP,
    std::vector<Eigen::MatrixXcs>* dQ,
    const Parameterization& para_x,
    const Parameterization& para_y,
    Eigen::MatrixXcs * dPl,
    Eigen::MatrixXcs * dQl)
{

    // auto t1 = sploosh::now();
    if (para_y.size() != 0)
    {
        cerr << "[ERROR]: currently not support variation of y\n";
    }

    using std::max;
    int nPara = max(count(para_x), count(para_y));
    // cout << "number of parameters: " << nPara << endl;

    vector<MatrixXcs> dfe11, dfe22, dfe33;
  
	FourierSolver2D e11Solver(coordX_, coordY_);
	MatrixXcs fe11;
	e11Solver.solve(eps_, nx, ny, DiscontinuousX, fe11,
                &dfe11, para_x, para_y);
  // ofstream fout1("DiscontinousX.txt");
  // fout1 << fe11 << endl;

	FourierSolver2D e22Solver(coordX_, coordY_);
	MatrixXcs fe22;
	e22Solver.solve(eps_, nx, ny, DiscontinuousY, fe22,
                &dfe22, para_x, para_y);
  // ofstream fout2("DiscontinousY.txt");
  // fout2 << fe22 << endl;

	FourierSolver2D e33Solver(coordX_, coordY_);
	MatrixXcs fe33;
	e33Solver.solve(eps_, nx, ny, ContinuousXY, fe33,
                &dfe33, para_x, para_y);
  int blockDim = nx * ny;
	//std::clog << "begin to solve linear system with e33...";
	BDCSVD<MatrixXcs> SVDSolver(fe33, ComputeThinU | ComputeThinV);
	const MatrixXcs& F1h = SVDSolver.solve(Kx);
	const MatrixXcs& F2h = SVDSolver.solve(Ky);

	// 透磁率の畳み込み行列。非磁性層 (μ=1) では bt_mu_11 / tb_mu_22 /
	// inv_mu_33 がすべて単位行列になるので、追加の Fourier 変換と SVD を
	// 省略できる (既定の高速経路)。
	MatrixXcs fmu11, fmu22, G1h, G2h;
	if (isMagnetic_)
	{
		FourierSolver2D m11Solver(coordX_, coordY_);
		m11Solver.solve(mu_, nx, ny, DiscontinuousX, fmu11);

		FourierSolver2D m22Solver(coordX_, coordY_);
		m22Solver.solve(mu_, nx, ny, DiscontinuousY, fmu22);

		FourierSolver2D m33Solver(coordX_, coordY_);
		MatrixXcs fmu33;
		m33Solver.solve(mu_, nx, ny, ContinuousXY, fmu33);

		BDCSVD<MatrixXcs> muSolver(fmu33, ComputeThinU | ComputeThinV);
		G1h = muSolver.solve(Kx);   // μzz⁻¹ Kx
		G2h = muSolver.solve(Ky);   // μzz⁻¹ Ky
	}

    // auto t3 = sploosh::now();
	//std::clog << "done!\n";

	P.resize(2*blockDim, 2*blockDim);
	scalar k0 = 2 * Pi / lambda;
	scalex k02 = scalex(k0 * k0);

	// P = [[ Kx εzz⁻¹ Ky,          k0² μyy − Kx εzz⁻¹ Kx ],
	//      [ Ky εzz⁻¹ Ky − k0² μxx, −Ky εzz⁻¹ Kx          ]]
	P.topLeftCorner(blockDim, blockDim) = Kx * F2h;
	P.topRightCorner(blockDim, blockDim) = - Kx * F1h;
	P.bottomLeftCorner(blockDim, blockDim) = Ky * F2h;
	P.bottomRightCorner(blockDim, blockDim) = -Ky * F1h;
	if (isMagnetic_)
	{
		P.topRightCorner(blockDim, blockDim) += k02 * fmu22;
		P.bottomLeftCorner(blockDim, blockDim) -= k02 * fmu11;
	}
	else
	{
		P.topRightCorner(blockDim, blockDim).diagonal().array() += k02;
		P.bottomLeftCorner(blockDim, blockDim).diagonal().array() -= k02;
	}

    // auto t4 = sploosh::now();

    // change of this parameter can get higher-order derivatives
    scalex dk02 = scalex(-8.*Pi*Pi/(lambda * lambda * lambda));

    if (dPl != nullptr)
    {
        dPl->resize(2*blockDim, 2*blockDim);
        dPl->setZero();

        // TODO may need optimization, now use a dense matrix
        dPl->topRightCorner(blockDim, blockDim).diagonal().array() += dk02;
        dPl->bottomLeftCorner(blockDim, blockDim).diagonal().array() -= dk02; 
    }

    if (dP != nullptr && nPara > 0)
    {
        dP->resize(nPara);
        for (int p = 0; p < nPara; ++p)
        {
            const MatrixXcs& temp = SVDSolver.solve(dfe33[p]);
            MatrixXcs& dp = (*dP)[p];
            dp.resize(2*blockDim, 2*blockDim);
            dp.topLeftCorner(blockDim, blockDim) = Kx * (-temp * F2h);
            dp.topRightCorner(blockDim, blockDim) = - Kx * (-temp * F1h);
            dp.bottomLeftCorner(blockDim, blockDim) = Ky * (-temp * F2h);
            dp.bottomRightCorner(blockDim, blockDim) = -Ky * (-temp * F1h);
        }
    }

    // auto t5 = sploosh::now();

	// Q = [[ −Kx μzz⁻¹ Ky,          Kx μzz⁻¹ Kx − k0² εyy ],
	//      [ −Ky μzz⁻¹ Ky + k0² εxx, Ky μzz⁻¹ Kx           ]]
	// (P と ε↔μ で対称。μ=1 では μzz⁻¹Kx → Kx となり従来式に一致する)
	Q.resize(2*blockDim, 2*blockDim);
	if (isMagnetic_)
	{
		Q.topLeftCorner(blockDim, blockDim) = -Kx * G2h;
		Q.topRightCorner(blockDim, blockDim) = Kx * G1h - k02 * fe22;
		Q.bottomLeftCorner(blockDim, blockDim) = -Ky * G2h + k02 * fe11;
		Q.bottomRightCorner(blockDim, blockDim) = Ky * G1h;
	}
	else
	{
		Q.topLeftCorner(blockDim, blockDim) = -Kx * Ky;
		Q.topRightCorner(blockDim, blockDim) = Kx * Kx - k02 * fe22;
		Q.bottomLeftCorner(blockDim, blockDim) = -Ky * Ky + k02 * fe11;
		Q.bottomRightCorner(blockDim, blockDim) = Ky * Kx;
	}

    // auto t6 = sploosh::now();

    if (dQl != nullptr)
    {
        dQl->resize(2*blockDim, 2*blockDim);
        dQl->setZero();

        dQl->topRightCorner(blockDim, blockDim) = -dk02 * fe22;
        dQl->bottomLeftCorner(blockDim, blockDim) = dk02 * fe11;
    }

    if (dQ != nullptr && nPara > 0)
    {
        dQ->resize(nPara);
        for (int p = 0; p < nPara; ++p)
        {
            MatrixXcs& dq = (*dQ)[p];
            dq.resize(2*blockDim, 2*blockDim);
            dq.setZero();
            dq.topRightCorner(blockDim, blockDim) = -k02 * dfe22[p];
            dq.bottomLeftCorner(blockDim, blockDim) = k02 * dfe11[p];
        }
    }
    // auto t7 = sploosh::now();

    // cout << "wavecoeff timing: " << sploosh::duration_milli_d(t1, t2) << ", "
    // << sploosh::duration_milli_d(t2, t3) << ", "
    // << sploosh::duration_milli_d(t3, t4) << ", "
    // << sploosh::duration_milli_d(t4, t5) << ", "
    // << sploosh::duration_milli_d(t5, t6) << ", "
    // << sploosh::duration_milli_d(t6, t7) << endl;
}

void Layer::solve(const Eigen::MatrixXcs& Kx, 
	const Eigen::MatrixXcs& Ky,
	scalar lambda,
	int nx, int ny)
{
	scalar k0 = 2 * Pi / lambda;

    MatrixXcs F, G;
    waveEqnCoeff(Kx, Ky, lambda, nx, ny, F, G);


    WaveEquationSolver solver;
    MatrixXcs FG;
    solver.solve(lambda, F, G, FG, eigvecE_, eigvecH_, gamma_);
    gamma_ *= 1./ k0;

	lambda_ = lambda;
	isSolved_ = true;
}

Eigen::MatrixXcs Layer::epsConvolution(int nx, int ny) const
{
	// waveEqnCoeff の fe33 と同じ ContinuousXY 規則による誘電率の畳み込み行列。
	// Ez = i·[[ε]]⁻¹(kx·Hy − ky·Hx) で用いるため再計算する。
	FourierSolver2D e33Solver(coordX_, coordY_);
	MatrixXcs fe33;
	e33Solver.solve(eps_, nx, ny, ContinuousXY, fe33);
	return fe33;
}


void Layer::permuteEigVecX(
	Eigen::MatrixXcs& eigvecE,
	Eigen::MatrixXcs& eigvecH,
	scalar dx, 
	int nx, 
	int ny)
{
	if (!isSolved())
	{
		std::cerr << "The layer has not been solved!\n";
		return;
	}


	MatrixXcs eigvecEtemp = eigvecE_;
	MatrixXcs eigvecHtemp = eigvecH_;


	// eigvecE_ and eigvecH_: both (2*nx*ny) x (2*nx*ny) matrices
	int maxX = (nx - 1) / 2;
	for (int m = -maxX; m <= maxX; ++m)
	{
		scalex perturb = exp(scalex(0, -2.* Pi * m / Lx() * dx));
		int i = m + maxX;
		for (int j = 0; j < ny; ++j)
		{
			int row1 = i * ny + j;
			int row2 = row1 + nx * ny;

			eigvecEtemp.row(row1) *= perturb;
			eigvecEtemp.row(row2) *= perturb;

			eigvecHtemp.row(row1) *= perturb;
			eigvecHtemp.row(row2) *= perturb;
		}
	}

	int nRow = eigvecEtemp.rows();
	int nCol = eigvecEtemp.cols();


	eigvecE.resize(nRow, nCol);
	eigvecH.resize(nRow, nCol);

	eigvecE = eigvecEtemp;
	eigvecH = eigvecHtemp;
}

void Layer::dpermuteEigVecX(
	Eigen::MatrixXcs& deigvecE,
	Eigen::MatrixXcs& deigvecH,
	scalar dx,
	int nx, 
	int ny)
{
	if (!isSolved())
	{
		std::cerr << "The layer has not been solved!\n";
		return;
	}


	MatrixXcs eigvecEtemp = eigvecE_;
	MatrixXcs eigvecHtemp = eigvecH_;


	// eigvecE_ and eigvecH_: both (2*nx*ny) x (2*nx*ny) matrices
	int maxX = (nx - 1) / 2;
	for (int m = -maxX; m <= maxX; ++m)
	{
		scalex perturb = scalex(0, -2 * Pi * m / Lx())
		 * exp(scalex(0, -2.* Pi * m / Lx() * dx));
		int i = m + maxX;
		for (int j = 0; j < ny; ++j)
		{
			int row1 = i * ny + j;
			int row2 = row1 + nx * ny;

			eigvecEtemp.row(row1) *= perturb;
			eigvecEtemp.row(row2) *= perturb;

			eigvecHtemp.row(row1) *= perturb;
			eigvecHtemp.row(row2) *= perturb;
		}
	}

	int nRow = eigvecEtemp.rows();
	int nCol = eigvecEtemp.cols();


	deigvecE.resize(nRow, nCol);
	deigvecH.resize(nRow, nCol);

	deigvecE = eigvecEtemp;
	deigvecH = eigvecHtemp;	
}

scalex Layer::eps(scalar x, scalar y) const
{
	// consider periodc case along x axis
	while(x < l0())
	{
		x += Lx();
	}

	while(x > r1())
	{
		x -= Lx();
	}
	
	// if (x < l0() || x > r1())
	// {
	// 	std::cerr << "The coordinate x is out of range!\n";
	// 	return scalex(.0);
	// }

	if (y < b0() || y > u1())
	{
		std::cerr << "The coordinate y is out of range!\n";
		return scalex(.0);
	}

	int nx = -1, ny = -1;

	for (int i = 0; i < coordX_.size()-1; ++i)
	{
		if (x <= coordX_(i+1))
		{
			nx = i;
			break;
		}
	}

	for (int j = 0; j < coordY_.size()-1; ++j)
	{
		if (y <= coordY_(j+1))
		{
			ny = j;
			break;
		}
	}

	return eps_(ny, nx);
}

scalex Layer::mu(scalar x, scalar y) const
{
	// eps(x, y) と同じセル探索。mu_ は eps_ と同形状。
	while (x < l0()) x += Lx();
	while (x > r1()) x -= Lx();

	if (y < b0() || y > u1())
	{
		std::cerr << "The coordinate y is out of range!\n";
		return scalex(.0);
	}

	int nx = -1, ny = -1;
	for (int i = 0; i < coordX_.size()-1; ++i)
		if (x <= coordX_(i+1)) { nx = i; break; }
	for (int j = 0; j < coordY_.size()-1; ++j)
		if (y <= coordY_(j+1)) { ny = j; break; }

	return mu_(ny, nx);
}


void Layer::generatePlaneWave(
	const Eigen::VectorXcs& delta,
	scalex px,
	scalex py,
	Eigen::VectorXcs& c)
{
	if (!isSolved_)
	{
		std::cerr << "The layer has not been solved!\n";
		return;
	}

	int nDim = gamma_.size() / 2;

	if (delta.size() != nDim)
	{
		std::cerr << "The size of harmonics does not match!\n";
		return;
	}

	VectorXcs s(2*nDim);
	s.head(nDim) = px * delta;
	s.tail(nDim) = py * delta;

	MatrixXcs Wref = eigvecE_;
	BDCSVD<MatrixXcs> csolver(Wref, ComputeThinU | ComputeThinV);
	c = csolver.solve(s);
}


void Layer::getHarmonics(
	const Eigen::VectorXcs& c,
	Eigen::VectorXcs& harmonics)
{
	harmonics = eigvecE_ * c;
}

void Layer::setSolutions(
	const Eigen::MatrixXcs& eigvecE,
	const Eigen::MatrixXcs& eigvecH,
	const Eigen::VectorXcs& gamma)
{
	eigvecE_ = eigvecE;
	eigvecH_ = eigvecH;
	gamma_ = gamma;

	isSolved_ = true;
}
