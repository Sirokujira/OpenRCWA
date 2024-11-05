#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"
#include "ev.h"

#include "hdf5.h"

#define FILE_NAME "time_series_data.h5"

#define GROUP_NAME "/metadata"
#define DATA_GROUP_PREFIX "/data"
//#define TIME_STEP 200 // Solver.nout
//#define MAX_STEP 1000 // Solver.maxiter

void read_metadata_from_hdf5() {
    // HDF5ファイルのオープン
    hid_t file_id = H5Fopen(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    herr_t status;
    hid_t group_id;
    hid_t dataset_id;

    //// /metadata グループのオープン
    //hid_t metadata_group_id = H5Gopen(file_id, "/metadata", H5P_DEFAULT);
    //
    //// Titleの読み込み
    //hid_t dataset_id = H5Dopen(metadata_group_id, "Title", H5P_DEFAULT);
    //status = H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
    //H5Dclose(dataset_id);
    //
    //// 各種整数型メタデータの読み込み
    //struct {
    //    const char *name;
    //    void *value;
    //    hid_t type;
    //} metadata[] = {
    //    {"Nx", &Nx, H5T_NATIVE_INT},
    //    {"Ny", &Ny, H5T_NATIVE_INT},
    //    {"Nz", &Nz, H5T_NATIVE_INT},
    //    {"Ni", &Ni, H5T_NATIVE_INT},
    //    {"Nj", &Nj, H5T_NATIVE_INT},
    //    {"Nk", &Nk, H5T_NATIVE_INT},
    //    {"N0", &N0, H5T_NATIVE_INT},
    //    {"NN", &NN, H5T_NATIVE_INT64},
    //    {"NFreq1", &NFreq1, H5T_NATIVE_INT},
    //    {"NFreq2", &NFreq2, H5T_NATIVE_INT},
    //    {"NFeed", &NFeed, H5T_NATIVE_INT},
    //    {"NPoint", &NPoint, H5T_NATIVE_INT},
    //    {"Niter", &Niter, H5T_NATIVE_INT},
    //    {"Ntime", &Ntime, H5T_NATIVE_INT},
    //    {"Solver_maxiter", &Solver.maxiter, H5T_NATIVE_INT},
    //    {"Solver_nout", &Solver.nout, H5T_NATIVE_INT},
    //    {"NGline", &NGline, H5T_NATIVE_INT},
    //    {"IPlanewave", &IPlanewave, H5T_NATIVE_INT}
    //};
    //
    //for (int i = 0; i < sizeof(metadata) / sizeof(metadata[0]); i++) {
    //    dataset_id = H5Dopen(metadata_group_id, metadata[i].name, H5P_DEFAULT);
    //    status = H5Dread(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
    //    H5Dclose(dataset_id);
    //}
    //
    //// Dtの読み込み
    //dataset_id = H5Dopen(metadata_group_id, "Dt", H5P_DEFAULT);
    //status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    //H5Dclose(dataset_id);
    //
    //// Planewaveの読み込み
    //dataset_id = H5Dopen(metadata_group_id, "Planewave", H5P_DEFAULT);
    //status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
    //H5Dclose(dataset_id);
    //
    //// 配列データの読み込み
    //struct {
    //    const char *name;
    //    double **data;
    //    size_t size;
    //} arrays[] = {
    //    {"Xn", &Xn, Nx + 1},
    //    {"Yn", &Yn, Ny + 1},
    //    {"Zn", &Zn, Nz + 1},
    //    {"Xc", &Xc, Nx},
    //    {"Yc", &Yc, Ny},
    //    {"Zc", &Zc, Nz},
    //    {"Eiter", &Eiter, Niter},
    //    {"Hiter", &Hiter, Niter},
    //    {"VFeed", &VFeed, NFeed * (Solver.maxiter + 1)},
    //    {"IFeed", &IFeed, NFeed * (Solver.maxiter + 1)},
    //    {"VPoint", &VPoint, NPoint * (Solver.maxiter + 1)},
    //    {"Freq1", &Freq1, NFreq1},
    //    {"Freq2", &Freq2, NFreq2},
    //    {"Gline", &Gline, NGline * 2 * 3}
    //};
    //
    //for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
    //    *arrays[i].data = (double *)malloc(sizeof(double) * arrays[i].size);
    //    dataset_id = H5Dopen(metadata_group_id, arrays[i].name, H5P_DEFAULT);
    //    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *arrays[i].data);
    //    H5Dclose(dataset_id);
    //}

    //// Surfaceデータの読み込み
    //dataset_id = H5Dopen(metadata_group_id, "NSurface", H5P_DEFAULT);
    //status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NSurface);
    //H5Dclose(dataset_id);
    //
    //Surface = (double *)malloc(sizeof(double) * NSurface);
    //dataset_id = H5Dopen(metadata_group_id, "Surface", H5P_DEFAULT);
    //status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
    //H5Dclose(dataset_id);

    // 電場データ用のメモリ割り当て
    //int NFreq2 = 10; // 仮の周波数数
    //int NN = 100; // 仮のデータポイント数
    //float *cEx_r = (float *)malloc(NFreq2 * NN * sizeof(float));
    //float *cEx_i = (float *)malloc(NFreq2 * NN * sizeof(float));
    //float *cEy_r = (float *)malloc(NFreq2 * NN * sizeof(float));
    //float *cEy_i = (float *)malloc(NFreq2 * NN * sizeof(float));
    //float *cEz_r = (float *)malloc(NFreq2 * NN * sizeof(float));
    //float *cEz_i = (float *)malloc(NFreq2 * NN * sizeof(float));
    //if (!cEx_r || !cEx_i || !cEy_r || !cEy_i || !cEz_r || !cEz_i) {
    //    fprintf(stderr, "Memory allocation failed for electric field components\n");
    //    H5Gclose(group_id);
    //    H5Fclose(file_id);
    //    return 1;
    //}
    // メタデータグループのクローズ
    //H5Gclose(metadata_group_id);

    // 時系列データの読み込み(Niter)
    char data_group_name[256];
    for (int step = 0; step <= Ntime; step += Solver.nout) {
        snprintf(data_group_name, sizeof(data_group_name), "%s%06d", DATA_GROUP_PREFIX, step);
        group_id = H5Gopen(file_id, data_group_name, H5P_DEFAULT);
        if (group_id < 0) {
            fprintf(stderr, "Failed to open data group %s\n", data_group_name);
            continue;
        }

        // 電場データの読み込み
        dataset_id = H5Dopen(group_id, "E", H5P_DEFAULT);
        if (dataset_id >= 0) {
            hsize_t e_dims[4] = {1, NFreq2, NN, 6};
            hid_t dataspace_id = H5Dget_space(dataset_id);

            hsize_t mem_dims[1] = {6};
            hid_t memspace_id = H5Screate_simple(1, mem_dims, NULL);

            for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                int64_t n0 = ifreq * NN;
                for (int nn = 0; nn < NN; nn++) {
                    float e_value[6];

                    hsize_t e_offset[4] = {0, ifreq, nn, 0};
                    hsize_t e_count[4] = {1, 1, 1, 6};
                    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, e_offset, NULL, e_count, NULL);

                    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, e_value);
                    if (status < 0) {
                        fprintf(stderr, "Failed to read E data at step=%d, ifreq=%d, nn=%d\n", step, ifreq, nn);
                    } else {
                        cEx_r[n0 + nn] = e_value[0];
                        cEy_r[n0 + nn] = e_value[1];
                        cEz_r[n0 + nn] = e_value[2];
                        cEx_i[n0 + nn] = e_value[3];
                        cEy_i[n0 + nn] = e_value[4];
                        cEz_i[n0 + nn] = e_value[5];
                    }
                }
            }
            H5Sclose(memspace_id);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
        }
        // 磁場データの読み込み
        dataset_id = H5Dopen(group_id, "H", H5P_DEFAULT);
        if (dataset_id >= 0) {
            hsize_t e_dims[4] = {1, NFreq2, NN, 6};
            hid_t dataspace_id = H5Dget_space(dataset_id);

            hsize_t mem_dims[1] = {6};
            hid_t memspace_id = H5Screate_simple(1, mem_dims, NULL);

            for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                int64_t n0 = ifreq * NN;
                for (int nn = 0; nn < NN; nn++) {
                    float e_value[6];

                    hsize_t e_offset[4] = {0, ifreq, nn, 0};
                    hsize_t e_count[4] = {1, 1, 1, 6};
                    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, e_offset, NULL, e_count, NULL);

                    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, e_value);
                    if (status < 0) {
                        fprintf(stderr, "Failed to read E data at step=%d, ifreq=%d, nn=%d\n", step, ifreq, nn);
                    } else {
                        cHx_r[n0 + nn] = e_value[0];
                        cHy_r[n0 + nn] = e_value[1];
                        cHz_r[n0 + nn] = e_value[2];
                        cHx_i[n0 + nn] = e_value[3];
                        cHy_i[n0 + nn] = e_value[4];
                        cHz_i[n0 + nn] = e_value[5];
                    }
                }
            }
            H5Sclose(memspace_id);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
        }

        {
            // 時系列毎の計算・出力対応
            printf("step(%d)(start).\n", step);
            if (NNear2d) {
                //setup_near2d();
                calcNear2d();
                outputNear2d();
            }
            printf("step(%d)(end).\n", step);
        }


        H5Gclose(group_id);
    }

    printf("ev2d_output(start).\n");
    ev2d_file(!HTML, (!HTML ? FN_ev2d_1 : FN_ev2d_0));
    ev2d_output();
    printf("ev2d_output(end).\n");

    // HDF5ファイルのクローズ
    H5Fclose(file_id);
}

