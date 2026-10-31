/*
rcwa_bridge.cpp

.ofd の rcwa/rcwalayer キーで指定された格子スタック (1D または 2D) を
RCWA コア (rcwa/RCWASolver) で解き、frequency1 掃引の回折効率を
rcwa_efficiency.csv に出力する。

規約 (rcwa/ コアの実装に合わせる):
  - 斜め入射に対応 (planewave = theta phi ... の theta/phi を使う)。
    偏波基底と Bloch 波数の決め方は rcwa/RCWADriver.cpp と同一にすること。
    **正規化 (1/kzinc) と偏波ベクトルは対で決まっており、片方だけ変えると
    斜め入射で R+T≠1 になる。**
  - planewave の pol は使わない。RCWA モードは常に TE/TM 両方を出力する
    (CSV の列構成を固定するため)。
  - Layer::gamma() は k0*neff [rad/m] を返すため、solve() の thickness は
    物理厚 [m] をそのまま渡す
  - scatterPlaneWave() の REF/TRN は入射波数 k0*sqrt(eps_ref) で
    正規化された次数別の効率 (合計が全反射率/全透過率)
  - 誘電率は複素数。時間規約 exp(-iwt) なので損失は「正」の虚部。
  - rcwadispersion を指定した材料は周波数ごとに Lorentz 式で誘電率を求める
    (rcwa/RCWADriver.cpp の epsAt と同じ式・同じ規約)。
*/

#include "orcwa_rcwa.h"
#include "rcwa_hdf5.h"

#include "rcwa/RCWASolver.h"
#include "rcwa/Layer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

extern "C" {
extern int    NFreq1;
extern double *Freq1;
extern char   Title[];
extern void   monitor1(FILE *, const char []);
}

namespace {

const double C0 = 2.99792458e8;  /* 真空中の光速 [m/s] */

/* 全層で共通のセル境界 [0, ..., period]。各層の fill*period を全部含める。
   RCWASolver::addLayer() は全層のセル境界 (先頭/末尾のセル) が一致することを
   要求するため、層ごとに自分の fill だけで境界を作ると、fill の違う層が
   混ざった時点で "The new layer does not match layers before!" で層が捨てられ、
   その後の solve() で落ちていた (fill が全層同じ既存サンプルでは顕在化しない)。
   区分一定の誘電率を同じ値のセルに細分しても Fourier 係数は変わらない。 */
std::vector<double> commonEdges(double period, bool useY)
{
	std::vector<double> e;
	e.push_back(0.0);
	for (int n = 0; n < NRcwaLayer; n++) {
		const double f = useY ? RcwaLayer[n].fill_y : RcwaLayer[n].fill;
		if (f < 1.0) e.push_back(f * period);
	}
	e.push_back(period);
	std::sort(e.begin(), e.end());
	std::vector<double> u;
	for (double v : e) {
		if (u.empty() || (v - u.back()) > 1e-12 * period) u.push_back(v);
	}
	return u;
}

/* 1層 = 1 周期に eps1 の矩形 [0, fill*Λx] × [0, fill_y*Λy] が 1 つ、残りが eps2。
   fill_y = 1 なら y 方向に一様 (1D 線格子)。xe / ye は commonEdges() の境界。
   eps は (y セル行 × x セル列) */
Layer makeLayer(const rcwalayer_t &l, scalex e1, scalex e2,
                const std::vector<double> &xe, const std::vector<double> &ye)
{
	const double px = xe.back(), py = ye.back();
	Eigen::VectorXs cx(xe.size()), cy(ye.size());
	for (size_t i = 0; i < xe.size(); i++) cx[i] = xe[i];
	for (size_t j = 0; j < ye.size(); j++) cy[j] = ye[j];
	Eigen::MatrixXcs eps(ye.size() - 1, xe.size() - 1);
	for (int j = 0; j < eps.rows(); j++) {
		const double yc = 0.5 * (ye[j] + ye[j + 1]);
		for (int i = 0; i < eps.cols(); i++) {
			const double xc = 0.5 * (xe[i] + xe[i + 1]);
			const bool in1 = (xc < l.fill * px) && (yc < l.fill_y * py);
			eps(j, i) = in1 ? e1 : e2;
		}
	}
	return Layer(cx, cy, eps);
}

/* (層 n, side) に有効な分散 (ae ≠ 0 の極) があるか */
bool hasDispersion(int n, int side)
{
	for (int p = 0; p < NRcwaPole; p++) {
		const rcwapole_t &q = RcwaPole[p];
		if ((q.layer == n) && (q.side == side) && (q.ae != 0.0)) return true;
	}
	return false;
}

/* 層 n の eps1 (side=1) / eps2 (side=2) の角周波数 omega [rad/s] での値。
   分散が指定されていれば Lorentz 式 (rcwa/RCWADriver.cpp の epsAt と同じ):
     eps = einf + Σ ae² / (ce² − ω² − i·be·ω)   (exp(-iωt) で be > 0 が損失)
   einf は最後に書いた行の値、ae = 0 の行は極として数えない。 */
scalex layerEps(int n, int side, double omega)
{
	const rcwalayer_t &l = RcwaLayer[n];
	if (!hasDispersion(n, side)) {
		return (side == 1) ? scalex(l.eps1, l.eps1i) : scalex(l.eps2, l.eps2i);
	}
	double einf = 0.0;
	for (int p = 0; p < NRcwaPole; p++) {
		const rcwapole_t &q = RcwaPole[p];
		if ((q.layer == n) && (q.side == side)) einf = q.einf;
	}
	scalex e(einf, 0.0);
	for (int p = 0; p < NRcwaPole; p++) {
		const rcwapole_t &q = RcwaPole[p];
		if ((q.layer != n) || (q.side != side) || (q.ae == 0.0)) continue;
		e += scalex(q.ae * q.ae, 0.0) / scalex(q.ce * q.ce - omega * omega, -q.be * omega);
	}
	return e;
}

/* 層が無損失か。静的な値は虚部ゼロ、分散は全極 be = 0 (Drude の衝突なし等) */
bool isLossless(int n)
{
	for (int side = 1; side <= 2; side++) {
		if (hasDispersion(n, side)) {
			for (int p = 0; p < NRcwaPole; p++) {
				const rcwapole_t &q = RcwaPole[p];
				if ((q.layer == n) && (q.side == side) && (q.ae != 0.0) && (q.be != 0.0))
					return false;
			}
		} else {
			const double im = (side == 1) ? RcwaLayer[n].eps1i : RcwaLayer[n].eps2i;
			if (im != 0.0) return false;
		}
	}
	return true;
}

} // namespace

extern "C" int rcwa_run(FILE *fp_log)
{
	char str[BUFSIZ];

	if ((NRcwaLayer < 2) || (RcwaPeriod <= 0) || (NFreq1 < 1)) {
		printf("%s\n", "*** invalid rcwa setup");
		return 1;
	}

	if (NRcwaHarmonicsY > 0) {
		sprintf(str, "RCWA : 2D Nx=%d Ny=%d (orders %dx%d), period=%.4e x %.4e[m],"
			" %d layers, %d frequencies",
			NRcwaHarmonics, NRcwaHarmonicsY, 2 * NRcwaHarmonics + 1,
			2 * NRcwaHarmonicsY + 1, RcwaPeriod, RcwaPeriodY, NRcwaLayer, NFreq1);
	} else {
		sprintf(str, "RCWA : N=%d (order %d), period=%.4e[m], %d layers, %d frequencies",
			NRcwaHarmonics, 2 * NRcwaHarmonics + 1, RcwaPeriod, NRcwaLayer, NFreq1);
	}
	monitor1(fp_log, str);
	sprintf(str, "RCWA : incidence theta=%.4g[deg] phi=%.4g[deg]", RcwaTheta, RcwaPhi);
	monitor1(fp_log, str);

	/* 損失層があると R+T<1 が正常なので、エネルギー保存の警告を切り替える */
	bool lossless = true;
	int ndisp = 0;
	for (int n = 0; n < NRcwaLayer; n++) {
		if (!isLossless(n)) lossless = false;
		ndisp += hasDispersion(n, 1) + hasDispersion(n, 2);
	}
	if (ndisp > 0) {
		sprintf(str, "RCWA : dispersive (Lorentz) materials : %d (layer, eps1/eps2) pairs", ndisp);
		monitor1(fp_log, str);
	}

	const double thetaRad = RcwaTheta * M_PI / 180.0;
	const double phiRad   = RcwaPhi   * M_PI / 180.0;

	/* 入射 E 場の横 (xy) 成分。単位振幅 (|E_inc_3D|²=1) の直交基底:
	     TM (p): 入射面内     → ( cosθ·cosφ,  cosθ·sinφ)
	     TE (s): 入射面に直交 → (−sinφ,       cosφ)
	   法線入射では両者が縮退するので x/y 方向をそのまま使う。
	   ここは rcwa/RCWADriver.cpp と同一でなければならない。 */
	double tmX, tmY, teX, teY;
	if (std::fabs(std::sin(thetaRad)) < 1e-9) {
		tmX = 1.0; tmY = 0.0;
		teX = 0.0; teY = 1.0;
	} else {
		tmX = std::cos(thetaRad) * std::cos(phiRad);
		tmY = std::cos(thetaRad) * std::sin(phiRad);
		teX = -std::sin(phiRad);
		teY =  std::cos(phiRad);
	}


	const std::vector<double> xedges = commonEdges(RcwaPeriod, false);
	const std::vector<double> yedges = commonEdges(RcwaPeriodY, true);

	std::vector<int> stack(NRcwaLayer);
	std::vector<scalar> thickness(NRcwaLayer);
	for (int n = 0; n < NRcwaLayer; n++) {
		stack[n] = n;
		/* 先頭/末尾は半無限層: 厚さは位相に使われない */
		thickness[n] = ((n == 0) || (n == NRcwaLayer - 1))
			? 0.0 : RcwaLayer[n].thickness;
	}

	const char csvname[] = "rcwa_efficiency.csv";
	FILE *csv = fopen(csvname, "w");
	if (csv == NULL) {
		printf("*** file %s open error.\n", csvname);
		return 1;
	}
	fprintf(csv, "frequency[Hz],lambda[m],R_TE,T_TE,R_TM,T_TM\n");

	/* HDF5 出力用に全点を溜める (CSV は従来どおり逐次書き出す)。
	   レイアウトは [npol][nfreq] の行優先: 先に TE の全周波数、次に TM。 */
	std::vector<rcwa_spectrum_row_t> spectrum(2 * (size_t)NFreq1);
	int nwritten = 0;

	int ierr = 0;
	for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
		const double freq = Freq1[ifreq];
		if (freq <= 0) {
			printf("%s\n", "*** invalid frequency1 data");
			ierr = 1;
			break;
		}
		const double lambda = C0 / freq;
		const double k0 = 2.0 * M_PI / lambda;
		const double omega = 2.0 * M_PI * freq;

		/* 入射側半無限層の屈折率 (Bloch 波数に使う)。両端は均質前提なので eps1 を見る。
		   複素平方根の実部を取ることで eps<1 (プラズマ等) にも対応する。
		   分散があると周波数ごとに変わるので、ここで求める。 */
		const scalex eps_inc = layerEps(0, 1, omega);
		const double n_inc = std::max(std::sqrt(eps_inc).real(), 1e-6);
		if (ndisp > 0) {
			static bool warned = false;
			const scalex eps_trn = layerEps(NRcwaLayer - 1, 1, omega);
			if (!warned && ((eps_inc.imag() > 1e-6 * std::fabs(eps_inc.real())) ||
			                (eps_trn.imag() > 1e-6 * std::fabs(eps_trn.real())))) {
				sprintf(str, "*** warning : a semi-infinite layer is lossy at f=%.4e[Hz];"
					" R/T normalization assumes lossless outer media", freq);
				monitor1(fp_log, str);
				warned = true;
			}
		}

		/* 層の固有値問題は波長依存のため周波数毎に解き直す */
		RCWASolver solver(NRcwaHarmonics, NRcwaHarmonicsY);
		solver.disablePML();
		for (int n = 0; n < NRcwaLayer; n++) {
			solver.addLayer(makeLayer(RcwaLayer[n], layerEps(n, 1, omega),
			                          layerEps(n, 2, omega), xedges, yedges));
		}
		/* 斜め入射の Bloch 波数。solve() より前に設定する必要がある */
		solver.setBlochWavevector(k0 * n_inc * std::sin(thetaRad) * std::cos(phiRad),
		                          k0 * n_inc * std::sin(thetaRad) * std::sin(phiRad));
		solver.solve(lambda, stack, thickness);

		/* 法線入射では TE = E//y (格子溝方向)、TM = E//x に退化する */
		double R[2], T[2];
		for (int pol = 0; pol < 2; pol++) {
			const double px = (pol == 0) ? teX : tmX;
			const double py = (pol == 0) ? teY : tmY;
			Eigen::VectorXcs cInc;
			solver.generateHorizontalPlaneWave(px, py, 0, cInc);
			Eigen::VectorXs REF, TRN;
			solver.scatterPlaneWave(cInc, 0, NRcwaLayer - 1, k0, REF, TRN);
			R[pol] = REF.sum();
			T[pol] = TRN.sum();
		}

		fprintf(csv, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
			freq, lambda, R[0], T[0], R[1], T[1]);

		for (int pol = 0; pol < 2; pol++) {
			rcwa_spectrum_row_t &row =
				spectrum[(size_t)pol * (size_t)NFreq1 + (size_t)ifreq];
			row.frequency = freq;
			row.lambda    = lambda;
			row.R = R[pol];
			row.T = T[pol];
			row.A = 1.0 - R[pol] - T[pol];
		}
		nwritten++;

		sprintf(str, "  f=%.4e[Hz] lambda=%.4e[m] R/T(TE)=%.5f/%.5f R/T(TM)=%.5f/%.5f",
			freq, lambda, R[0], T[0], R[1], T[1]);
		monitor1(fp_log, str);

		/* 無損失入力に対するエネルギー保存の破れは実装異常。
		   損失層があると R+T<1 が正常なので、その場合は R+T>1 (利得) のみ疑う。 */
		for (int pol = 0; pol < 2; pol++) {
			const double sum = R[pol] + T[pol];
			if (lossless) {
				if (std::fabs(sum - 1.0) > 1e-3) {
					sprintf(str, "*** warning : R+T=%.6f (energy not conserved)", sum);
					monitor1(fp_log, str);
				}
			} else if (sum > 1.0 + 1e-3) {
				sprintf(str, "*** warning : R+T=%.6f > 1 with lossy layers"
					" (check the sign of the eps imaginary part)", sum);
				monitor1(fp_log, str);
			}
		}
	}

	fclose(csv);

	if (!ierr) {
		sprintf(str, "output : %s", csvname);
		monitor1(fp_log, str);

		/* GUI 表示用の HDF5。CSV と並行して出力する (既存の流れは変えない)。
		   書けなくても計算結果は CSV にあるので致命的にはしない。 */
		const char h5name[] = "time_series_data.h5";
		rcwa_meta_t meta;
		meta.title      = Title;
		meta.nharmonics   = NRcwaHarmonics;
		meta.period       = RcwaPeriod;
		meta.nharmonics_y = NRcwaHarmonicsY;
		meta.period_y     = (NRcwaHarmonicsY > 0) ? RcwaPeriodY : 0.0;
		meta.nlayer     = NRcwaLayer;
		meta.theta      = RcwaTheta;
		meta.phi        = RcwaPhi;
		const char *pol_labels[2] = { "TE", "TM" };
		/* 途中で break した場合に備え、実際に埋まった点数だけ詰め直す */
		if (nwritten < NFreq1) {
			for (int i = 0; i < nwritten; i++)
				spectrum[(size_t)nwritten + (size_t)i] =
					spectrum[(size_t)NFreq1 + (size_t)i];
		}
		if (rcwa_write_hdf5(h5name, &meta, spectrum.data(),
		                    2, nwritten, pol_labels) == 0) {
			sprintf(str, "output : %s", h5name);
		} else {
			sprintf(str, "*** warning : %s write failed (CSV is still valid)", h5name);
		}
		monitor1(fp_log, str);
	}

	return ierr;
}
