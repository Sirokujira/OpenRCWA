/*
outputCoupling.c

coupling
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"
#define GROUP_NAME "/metadata"
#define COUPLING_DATASET_NAME "coupling"

d_complex_t coupling(int ifeed, int ipoint, int ifreq)
{
	const d_complex_t cvf = calcdft(Ntime, &VFeed[ifeed   * (Solver.maxiter + 1)], Freq1[ifreq], Dt, 0);
	const d_complex_t cvp = calcdft(Ntime, &VPoint[ipoint * (Solver.maxiter + 1)], Freq1[ifreq], Dt, 0);

	return d_div(cvp, cvf);
}

typedef struct {
    double magnitude_dB;
    double phase_deg;
} coupling_data_t;

void write_coupling_data_to_hdf5()
{
    hid_t file_id, group_id, dataset_id, dataspace_id, datatype_id;
    herr_t status;

    // coupling_data_t 型のデータを保持する配列を作成
    coupling_data_t data[NFreq1][NFeed][NPoint];

    // coupling データの計算
    for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
        for (int ifeed = 0; ifeed < NFeed; ifeed++) {
            for (int ipoint = 0; ipoint < NPoint; ipoint++) {
                const d_complex_t couple = coupling(ifeed, ipoint, ifreq);
                data[ifreq][ifeed][ipoint].magnitude_dB = 20 * log10(fmax(d_abs(couple), EPS2));
                data[ifreq][ifeed][ipoint].phase_deg = d_deg(couple);
            }
        }
    }

    // HDF5 ファイルを開く
    file_id = H5Fopen(FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", FILE_NAME);
        return;
    }

    // グループを開く
    group_id = H5Gopen(file_id, GROUP_NAME, H5P_DEFAULT);
    if (group_id < 0) {
        fprintf(stderr, "Error opening group: %s\n", GROUP_NAME);
        H5Fclose(file_id);
        return;
    }

    // データ積の次元を定義
    hsize_t dims[3] = {NFreq1, NFeed, NPoint};
    dataspace_id = H5Screate_simple(3, dims, NULL);

    // coupling_data_t のデータ型を定義
    datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(coupling_data_t));
    H5Tinsert(datatype_id, "magnitude_dB", HOFFSET(coupling_data_t, magnitude_dB), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "phase_deg", HOFFSET(coupling_data_t, phase_deg), H5T_NATIVE_DOUBLE);

    // データセットを作成
    dataset_id = H5Dcreate(group_id, COUPLING_DATASET_NAME, datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset: %s\n", COUPLING_DATASET_NAME);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // coupling データを書き込み
    status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) {
        fprintf(stderr, "Error writing dataset: %s\n", COUPLING_DATASET_NAME);
    }

    // リソースを解放
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(datatype_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

static void _outputCoupling(FILE *fp)
{
	char str[BUFSIZ];

	fprintf(fp, "=== coupling ===\n");

	fprintf(fp, "  frequency[Hz]");
	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		for (int ipoint = 0; ipoint < NPoint; ipoint++) {
			sprintf(str, "C%d%d", ipoint + 1, ifeed + 1);
			fprintf(fp, "  %s[dB] %s[deg]", str, str);
		}
	}
	fprintf(fp, "\n");

	for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
		fprintf(fp, "  %13.5e", Freq1[ifreq]);
		for (int ifeed = 0; ifeed < NFeed; ifeed++) {
			for (int ipoint = 0; ipoint < NPoint; ipoint++) {
				const d_complex_t couple = coupling(ifeed, ipoint, ifreq);
				fprintf(fp, "%9.3f", 20 * log10(MAX(d_abs(couple), EPS2)));
				fprintf(fp, "%9.3f", d_deg(couple));
			}
		}
		fprintf(fp, "\n");
	}

	fflush(fp);
}


void outputCoupling(FILE *fp)
{
	_outputCoupling(stdout);
	_outputCoupling(fp);

	write_coupling_data_to_hdf5();
}

