/*
outputZin.c

input impedance
*/

#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

#include "hdf5.h"

#define FILE_NAME "time_series_data.h5"
#define GROUP_NAME "/metadata"
#define ZIN_DATASET_NAME "input_impedance"
typedef struct {
    double frequency;
    double rin;
    double xin;
    double gin;
    double bin;
    double ref;
    double vswr;
} input_impedance_data_t;

// calculate input impedance and reflection
void calcZin(void)
{
	if ((NFeed <= 0) || (NFreq1 <= 0)) return;

	Zin = (d_complex_t *)malloc(NFeed * NFreq1 * sizeof(d_complex_t));
	Ref =      (double *)malloc(NFeed * NFreq1 * sizeof(double));

	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		double *fv = &VFeed[ifeed * (Solver.maxiter + 1)];
		double *fi = &IFeed[ifeed * (Solver.maxiter + 1)];
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			const int id = (ifeed * NFreq1) + ifreq;

			// Zin
			const d_complex_t vin = calcdft(Ntime, fv, Freq1[ifreq], Dt, 0);
			const d_complex_t iin = calcdft(Ntime, fi, Freq1[ifreq], Dt, -0.5);
			Zin[id] = d_div(vin, iin);

			// Reflection = (Zin - Z0) / (Zin + Z0)
			const d_complex_t z0 = d_complex(Feed[ifeed].z0, 0);
			const d_complex_t ref = d_div(d_sub(Zin[id], z0), d_add(Zin[id], z0));
			Ref[id] = 10 * log10(d_norm(ref));
		}
	}
}


// calculate input Power (for far field gain : for post)
void calcPin(void)
{
	if ((NFeed <= 0) || (NFreq2 <= 0)) return;

	for (int i = 0; i < 2; i++) {
		Pin[i] = (double *)malloc(NFeed * NFreq2 * sizeof(double));
	}

	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		double *fv = &VFeed[ifeed * (Solver.maxiter + 1)];
		double *fi = &IFeed[ifeed * (Solver.maxiter + 1)];
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			const d_complex_t vin = d_div(calcdft(Ntime, fv, Freq2[ifreq], Dt, 0),    cFdft[ifreq]);
			const d_complex_t iin = d_div(calcdft(Ntime, fi, Freq2[ifreq], Dt, -0.5), cFdft[ifreq]);
			const d_complex_t zin = d_div(vin, iin);
			const double rin = zin.r;
			const double xin = zin.i;
			const double z0 = Feed[ifeed].z0;
			const double denom = 1
				 - ((rin - z0) * (rin - z0) + (xin * xin))
				 / ((rin + z0) * (rin + z0) + (xin * xin));
			const int id = (ifeed * NFreq2) + ifreq;
			Pin[0][id] = (vin.r * iin.r) + (vin.i * iin.i);
			Pin[1][id] = Pin[0][id] / MAX(denom, EPS);
		}
	}
}

static void write_input_impedance_data_to_hdf5()
{
    hid_t file_id, group_id, dataset_id, dataspace_id, memtype_id;
    herr_t status;

    // input_impedance_data_t 型のデータを保持する配列を作成
    input_impedance_data_t data[NFeed][NFreq1];

    // input impedance データの計算
    for (int ifeed = 0; ifeed < NFeed; ifeed++) {
        for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
            const int id = (ifeed * NFreq1) + ifreq;
            const d_complex_t yin = d_inv(Zin[id]);
            const double gamma = pow(10, Ref[id] / 20);
            const double vswr = (fabs(1 - gamma) > EPS) ? (1 + gamma) / (1 - gamma) : 1000;
            data[ifeed][ifreq].frequency = Freq1[ifreq];
            data[ifeed][ifreq].rin = Zin[id].r;
            data[ifeed][ifreq].xin = Zin[id].i;
            data[ifeed][ifreq].gin = yin.r * 1e3;
            data[ifeed][ifreq].bin = yin.i * 1e3;
            data[ifeed][ifreq].ref = Ref[id];
            data[ifeed][ifreq].vswr = vswr;
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

    // データ空間の次元を定義
    hsize_t dims[2] = {NFeed, NFreq1};
    dataspace_id = H5Screate_simple(2, dims, NULL);

    // メモリタイプを定義
    memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(input_impedance_data_t));
    H5Tinsert(memtype_id, "frequency", HOFFSET(input_impedance_data_t, frequency), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "rin", HOFFSET(input_impedance_data_t, rin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "xin", HOFFSET(input_impedance_data_t, xin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "gin", HOFFSET(input_impedance_data_t, gin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "bin", HOFFSET(input_impedance_data_t, bin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "ref", HOFFSET(input_impedance_data_t, ref), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "vswr", HOFFSET(input_impedance_data_t, vswr), H5T_NATIVE_DOUBLE);

    // データセットを作成
    dataset_id = H5Dcreate(group_id, ZIN_DATASET_NAME, memtype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset: %s\n", ZIN_DATASET_NAME);
        H5Sclose(dataspace_id);
        H5Tclose(memtype_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // input impedance データを書き込み
    status = H5Dwrite(dataset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) {
        fprintf(stderr, "Error writing dataset: %s\n", ZIN_DATASET_NAME);
    }

    // リソースを解放
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(memtype_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

static void _outputZin(FILE *fp)
{
	fprintf(fp, "=== input impedance ===\n");

	for (int ifeed = 0; ifeed < NFeed; ifeed++) {
		fprintf(fp, "feed #%d (Z0[ohm] = %.2f)\n", ifeed + 1, Feed[ifeed].z0);
		fprintf(fp, "  %s\n", "frequency[Hz] Rin[ohm]   Xin[ohm]    Gin[mS]    Bin[mS]    Ref[dB]       VSWR");
		for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
			const int id = (ifeed * NFreq1) + ifreq;
			const d_complex_t yin = d_inv(Zin[id]);
			const double gamma = pow(10, Ref[id] / 20);
			const double vswr = (fabs(1 - gamma) > EPS) ? (1 + gamma) / (1 - gamma) : 1000;
			fprintf(fp, "%13.5e%11.3f%11.3f%11.3f%11.3f%11.3f%11.3f\n",
				Freq1[ifreq], Zin[id].r, Zin[id].i, yin.r * 1e3, yin.i * 1e3, Ref[id], vswr);
		}
	}

	fflush(fp);
}


void outputZin(FILE *fp)
{
	_outputZin(stdout);
	_outputZin(fp);
	
	write_input_impedance_data_to_hdf5();
}

