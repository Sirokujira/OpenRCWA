/*
outputSpara.c

S-parameters
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"
#define GROUP_NAME "/metadata"
#define SPARA_DATASET_NAME "s_parameters"
typedef struct {
    double magnitude_dB;
    double phase_deg;
} spara_data_t;

void outputTouchstone(FILE *fp, int num_ports);

void calcSpara(void)
{
	if ((NPoint <= 0) || (NFreq1 <= 0)) return;

	// alloc
	d_complex_t *cv = (d_complex_t *)malloc((NPoint + 2) * NFreq1 * sizeof(d_complex_t));

	// DFT
	for (int ipoint = 0; ipoint < NPoint + 2; ipoint++) {
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			const int id = (ipoint * NFreq1) + ifreq;
			cv[id] = calcdft(Ntime, &VPoint[ipoint * (Solver.maxiter + 1)], Freq1[ifreq], Dt, 0);
		}
	}

	// S-parameters
	for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
		d_complex_t cv0 = cv[ 0           * NFreq1 + ifreq];  // V1
		d_complex_t cvp = cv[ NPoint      * NFreq1 + ifreq];  // V1+
		d_complex_t cvm = cv[(NPoint + 1) * NFreq1 + ifreq];  // V1-
		//printf("%d %f %f %f\n", ifreq, d_abs(cv0), d_abs(cvp), d_abs(cvm));
		//printf("%d %f %f %f %f %f %f\n", ifreq, cv0.r, cv0.i, cvp.r, cvp.i, cvm.r, cvm.i);
		d_complex_t c1 = d_div(d_add(cvp, cvm), cv0);
		d_complex_t c2 = d_sqrt(d_sub(d_mul(c1, c1), d_complex(4, 0)));
		d_complex_t c3 = d_add(c1, c2);
		if (c3.i < 0) c3 = d_sub(c1, c2);  // Im > 0
		c3 = d_rmul(0.5, c3);              // exp(+gd)
		d_complex_t c4 = d_div(d_complex(1, 0), c3);   // exp(-gd)
		d_complex_t c5 = d_sub(d_mul(c4, c4), d_mul(c3, c3));  // exp(-2gd) - exp(2gd)
		d_complex_t c6 = d_div(d_sub(d_mul(cvp, c4), d_mul(cvm, c3)), c5);  // V+
		d_complex_t c7 = d_div(d_sub(d_mul(cvm, c4), d_mul(cvp, c3)), c5);  // V-
		// S11 = E- / E+
		Spara[ifreq] = d_div(c7, c6);
		// Sn1 (n > 1)
		for (int ipoint = 1; ipoint < NPoint; ipoint++) {
			const int id = (ipoint * NFreq1) + ifreq;
			Spara[id] = d_div(cv[id], c6);  // Sn1 = Vn / V+
		}
	}

	// free
	free(cv);
}

static void _outputSpara(FILE *fp)
{
	fprintf(fp, "=== S-parameters ===\n");

	fprintf(fp, "  frequency[Hz]");
	for (int ipoint = 0; ipoint < NPoint; ipoint++) {
		fprintf(fp, "  S%d1[dB] S%d1[deg]", ipoint + 1, ipoint + 1);
	}
	fprintf(fp, "\n");

	for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
		fprintf(fp, "  %13.5e", Freq1[ifreq]);
		for (int ipoint = 0; ipoint < NPoint; ipoint++) {
			const int id = (ipoint * NFreq1) + ifreq;
			fprintf(fp, "%9.3f%9.3f", 20 * log10(MAX(d_abs(Spara[id]), EPS2)), d_deg(Spara[id]));
		}
		fprintf(fp, "\n");
	}

	fflush(fp);
}

static void write_spara_data_to_hdf5()
{
    hid_t file_id, group_id, dataset_id, dataspace_id, datatype_id;
    herr_t status;

    // spara_data_t 型のデータを保持する配列を作成 (MSVC は C99 VLA 非対応のため malloc で確保)
    spara_data_t *data = (spara_data_t *)malloc((size_t)NFreq1 * NPoint * sizeof(spara_data_t));
    if (data == NULL) {
        fprintf(stderr, "Error allocating memory for S-parameters data\n");
        return;
    }

    // S-parameters データの計算
    for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
        for (int ipoint = 0; ipoint < NPoint; ipoint++) {
            const int id = (ipoint * NFreq1) + ifreq;
            const size_t did = (size_t)ifreq * NPoint + ipoint;
            data[did].magnitude_dB = 20 * log10(MAX(d_abs(Spara[id]), EPS2));
            data[did].phase_deg = d_deg(Spara[id]);
        }
    }

    // HDF5 ファイルを開く
    file_id = H5Fopen(FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", FILE_NAME);
        free(data);
        return;
    }

    // グループを開く
    group_id = H5Gopen(file_id, GROUP_NAME, H5P_DEFAULT);
    if (group_id < 0) {
        fprintf(stderr, "Error opening group: %s\n", GROUP_NAME);
        H5Fclose(file_id);
        free(data);
        return;
    }

    // データタイプを作成
    datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(spara_data_t));
    H5Tinsert(datatype_id, "magnitude_dB", HOFFSET(spara_data_t, magnitude_dB), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "phase_deg", HOFFSET(spara_data_t, phase_deg), H5T_NATIVE_DOUBLE);

    // データ空間の次元を定義
    hsize_t dims[2] = {NFreq1, NPoint};
    dataspace_id = H5Screate_simple(2, dims, NULL);

    // データセットを作成
    dataset_id = H5Dcreate(group_id, SPARA_DATASET_NAME, datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset: %s\n", SPARA_DATASET_NAME);
    	H5Tclose(datatype_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        free(data);
        return;
    }

    // S-parameters データを書き込み
    // メモリ側も複合型で渡す (H5T_NATIVE_DOUBLE では型不一致になる)
    status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) {
        fprintf(stderr, "Error writing dataset: %s\n", SPARA_DATASET_NAME);
    }

    // リソースを解放
    H5Dclose(dataset_id);
	H5Tclose(datatype_id);
    H5Sclose(dataspace_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
    free(data);
}

// .s(NPoint)p
#define FN_snp      "test.snp"
void outputSpara(FILE *fp)
{
	_outputSpara(stdout);
	_outputSpara(fp);

	write_spara_data_to_hdf5();

	// input data
	FILE *fp_snp = NULL;
	int ierr = 0;
	if ((fp_snp = fopen(FN_snp, "w")) == NULL) {
		ierr = 1;
	}
	if (!ierr) {
		// Touchstone形式で出力
		// input  : NFeed?
		// output : NPoint
	    outputTouchstone(fp_snp, NPoint);

		fclose(fp_snp);
	}
	// fopen 失敗時 (fp_snp == NULL) は何も閉じない
}

void outputTouchstone(FILE *fp, int num_ports) {
    // ヘッダー情報の出力
    fprintf(fp, "# Hz S MA R 50\n");  // Sパラメータ、マグニチュード・角度表記、50Ω基準

    // 各周波数でのSパラメータを出力
    for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
        // 周波数の出力
        fprintf(fp, "%13.5e", Freq1[ifreq]);
        // 各ポート間のSパラメータを出力（異なる入力/出力ポートに対応）
        for (int input_port = 0; input_port < 1; input_port++) {
            for (int output_port = 0; output_port < num_ports; output_port++) {
                int index = (input_port * num_ports + output_port) * NFreq1 + ifreq;
                double magnitude = d_abs(Spara[index]);
                double phase = d_deg(Spara[index]);
                fprintf(fp, " %13.5e %13.5e", magnitude, phase);
            }
        }
        fprintf(fp, "\n");
    }

    fflush(fp);
}
