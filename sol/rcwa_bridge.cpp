/*
rcwa_bridge.cpp

.ofd の rcwa/rcwalayer キーで指定された 1D 格子スタックを
RCWA コア (rcwa/RCWASolver) で解き、frequency1 掃引の回折効率を
rcwa_efficiency.csv に出力する。

規約 (rcwa/ コアの実装に合わせる):
  - 垂直入射のみ (generateHorizontalPlaneWave)
  - Layer::gamma() は k0*neff [rad/m] を返すため、solve() の thickness は
    物理厚 [m] をそのまま渡す
  - scatterPlaneWave() の REF/TRN は入射波数 k0*sqrt(eps_ref) で
    正規化された次数別の効率 (合計が全反射率/全透過率)
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
	eps << scalex(l.eps1), scalex(l.eps2);
	return Layer(cx, cy, eps);
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
		solver.solve(lambda, stack, thickness);

		/* TE (E//y: 格子溝方向) と TM (E//x) の両偏波 */
		double R[2], T[2];
		for (int pol = 0; pol < 2; pol++) {
			Eigen::VectorXcs cInc;
			solver.generateHorizontalPlaneWave(
				(pol == 0) ? 0.0 : 1.0,
				(pol == 0) ? 1.0 : 0.0, 0, cInc);
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

		/* 無損失入力に対するエネルギー保存の破れは実装異常 */
		for (int pol = 0; pol < 2; pol++) {
			if (std::fabs(R[pol] + T[pol] - 1.0) > 1e-3) {
				sprintf(str, "*** warning : R+T=%.6f (energy not conserved)",
					R[pol] + T[pol]);
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
