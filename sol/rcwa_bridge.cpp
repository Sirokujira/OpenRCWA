/*
rcwa_bridge.cpp

.ofd の rcwa/rcwalayer キーで指定された 1D 格子スタックを
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
*/

#include "orcwa_rcwa.h"

#include "rcwa/RCWASolver.h"
#include "rcwa/Layer.h"

#include <cmath>
#include <cstdio>
#include <vector>

extern "C" {
extern int    NFreq1;
extern double *Freq1;
extern void   monitor1(FILE *, const char []);
}

namespace {

const double C0 = 2.99792458e8;  /* 真空中の光速 [m/s] */

/* 1層 = x 方向 2 セル (eps1: 0..fill*P, eps2: fill*P..P)、y 方向一様 */
Layer makeLayer(const rcwalayer_t &l, double period)
{
	Eigen::VectorXs cx(3);
	cx << 0.0, l.fill * period, period;
	Eigen::VectorXs cy(2);
	cy << 0.0, period;
	Eigen::MatrixXcs eps(1, 2);
	eps << scalex(l.eps1, l.eps1i), scalex(l.eps2, l.eps2i);
	return Layer(cx, cy, eps);
}

/* 層が無損失か (虚部がゼロか) */
bool isLossless(const rcwalayer_t &l)
{
	return (l.eps1i == 0.0) && (l.eps2i == 0.0);
}

} // namespace

extern "C" int rcwa_run(FILE *fp_log)
{
	char str[BUFSIZ];

	if ((NRcwaLayer < 2) || (RcwaPeriod <= 0) || (NFreq1 < 1)) {
		printf("%s\n", "*** invalid rcwa setup");
		return 1;
	}

	sprintf(str, "RCWA : N=%d (order %d), period=%.4e[m], %d layers, %d frequencies",
		NRcwaHarmonics, 2 * NRcwaHarmonics + 1, RcwaPeriod, NRcwaLayer, NFreq1);
	monitor1(fp_log, str);
	sprintf(str, "RCWA : incidence theta=%.4g[deg] phi=%.4g[deg]", RcwaTheta, RcwaPhi);
	monitor1(fp_log, str);

	/* 損失層があると R+T<1 が正常なので、エネルギー保存の警告を切り替える */
	bool lossless = true;
	for (int n = 0; n < NRcwaLayer; n++) {
		if (!isLossless(RcwaLayer[n])) lossless = false;
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

	/* 入射側半無限層の屈折率 (Bloch 波数に使う)。両端は均質前提なので eps1 を見る。
	   複素平方根の実部を取ることで eps<1 (プラズマ等) にも対応する。 */
	const scalex eps_inc = scalex(RcwaLayer[0].eps1, RcwaLayer[0].eps1i);
	const double n_inc = std::max(std::sqrt(eps_inc).real(), 1e-6);

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

		/* 層の固有値問題は波長依存のため周波数毎に解き直す */
		RCWASolver solver(NRcwaHarmonics, 0);
		solver.disablePML();
		for (int n = 0; n < NRcwaLayer; n++) {
			solver.addLayer(makeLayer(RcwaLayer[n], RcwaPeriod));
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
	}

	return ierr;
}
