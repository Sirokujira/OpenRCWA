#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"
#define GROUP_NAME "/metadata"
#define CROSS_SECTION_DATASET_NAME "cross_section"

typedef struct {
    double frequency;
    double backward;
    double forward;
} cross_section_data_t;

static void _outputCross(FILE *fp, const double bcs[], const double fcs[])
{
    fprintf(fp, "=== cross section ===\n");
    fprintf(fp, "  %s\n", "frequency[Hz] backward[m*m]  forward[m*m]");
    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
        fprintf(fp, "  %13.5e  %12.4e  %12.4e\n", Freq2[ifreq], bcs[ifreq], fcs[ifreq]);
    }
}

void write_cross_section_data_to_hdf5(const double bcs[], const double fcs[])
{
    hid_t file_id, group_id, dataset_id, dataspace_id, datatype_id;
    herr_t status;

    // cross_section_data_t 型のデータを保持する配列を作成
    cross_section_data_t data[NFreq2];

    // cross section データの計算
    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
        data[ifreq].frequency = Freq2[ifreq];
        data[ifreq].backward = bcs[ifreq];
        data[ifreq].forward = fcs[ifreq];
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
    hsize_t dims[1] = {NFreq2};
    dataspace_id = H5Screate_simple(1, dims, NULL);

    // cross_section_data_t 型のデータ型を定義
    datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(cross_section_data_t));
    H5Tinsert(datatype_id, "frequency", HOFFSET(cross_section_data_t, frequency), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "backward", HOFFSET(cross_section_data_t, backward), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "forward", HOFFSET(cross_section_data_t, forward), H5T_NATIVE_DOUBLE);

    // データセットを作成
    dataset_id = H5Dcreate(group_id, CROSS_SECTION_DATASET_NAME, datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset: %s\n", CROSS_SECTION_DATASET_NAME);
        H5Tclose(datatype_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // cross section データを書き込み
    status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) {
        fprintf(stderr, "Error writing dataset: %s\n", CROSS_SECTION_DATASET_NAME);
    }

    // リソースを解放
    H5Dclose(dataset_id);
    H5Tclose(datatype_id);
    H5Sclose(dataspace_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

void outputCross(FILE *fp)
{
    // alloc
    double *bcs = (double *)malloc(NFreq2 * sizeof(double));
    double *fcs = (double *)malloc(NFreq2 * sizeof(double));

    // calculation
    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
        const double ffctr = farfactor(ifreq);
        d_complex_t etheta, ephi;
        // BCS
        farfield(ifreq,       Planewave.theta,       Planewave.phi, ffctr, &etheta, &ephi);
        bcs[ifreq] = d_norm(etheta) + d_norm(ephi);
        // FCS
        farfield(ifreq, 180 - Planewave.theta, 180 + Planewave.phi, ffctr, &etheta, &ephi);
        fcs[ifreq] = d_norm(etheta) + d_norm(ephi);
    }

    // output
    _outputCross(stdout, bcs, fcs);
    _outputCross(fp,     bcs, fcs);

    write_cross_section_data_to_hdf5(bcs, fcs);

    // free
    free(bcs);
    free(fcs);
}

