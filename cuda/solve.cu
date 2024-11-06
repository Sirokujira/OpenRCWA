/*
solve.cu (CUDA)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"

#include "hdf5.h"
#define FILE_NAME "time_series_data.h5"

static void copy_to_host();

void solve(int io, double *tdft, FILE *fp)
{
    // HDF5ファイルの作成
    // 関数から?(fp の入替え?)
    hid_t file_id;
    // local
    hid_t group_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;

    double fmax[] = {0, 0};
    char   str[BUFSIZ];
    int    converged = 0;

    // setup host memory
    setup_host();

    // setup (GPU)
    if (GPU) {
        setup_gpu();
    }

    // initial field
    initfield();

    // HDF5ファイルの作成
    file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // time step iteration
    int itime;
    double t = 0;

    int lx = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
    int ly = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
    int lz = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;

    for (itime = 0; itime <= Solver.maxiter; itime++) {

        // update H
        t += 0.5 * Dt;
        updateHx(t);
        updateHy(t);
        updateHz(t);

        // ABC H
        if      (iABC == 0) {
            murH(numMurHx, (GPU ? d_fMurHx : fMurHx), Hx);
            murH(numMurHy, (GPU ? d_fMurHy : fMurHy), Hy);
            murH(numMurHz, (GPU ? d_fMurHz : fMurHz), Hz);
        }
        else if (iABC == 1) {
            pmlHx();
            pmlHy();
            pmlHz();
        }

        // PBC H
        if (PBCx) {
            pbcx();
        }
        if (PBCy) {
            pbcy();
        }
        if (PBCz) {
            pbcz();
        }

        // update E
        t += 0.5 * Dt;
        updateEx(t);
        updateEy(t);
        updateEz(t);

        // dispersion E
        if (numDispersionEx) {
            dispersionEx(t);
        }
        if (numDispersionEy) {
            dispersionEy(t);
        }
        if (numDispersionEz) {
            dispersionEz(t);
        }

        // ABC E
        if      (iABC == 1) {
            pmlEx();
            pmlEy();
            pmlEz();
        }

        // feed
        if (NFeed) {
            efeed(itime);
        }

        // inductor
        if (NInductor) {
            eload();
        }

        // point
        if (NPoint) {
            vpoint(itime);
        }

        // DFT
        if (GPU) cudaDeviceSynchronize();
        const double t0 = cputime();
        dftNear3d(itime);
        if (GPU) cudaDeviceSynchronize();
        *tdft += cputime() - t0;

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // average (plot)
            Eiter[Niter] = fsum[0];
            Hiter[Niter] = fsum[1];
            //Niter++;

            // monitor
            if (io) {
                sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
                fprintf(fp,     "%s\n", str);
                fprintf(stdout, "%s\n", str);
                fflush(fp);
                fflush(stdout);

                // copy near3d from device to host
                memcopy3_gpu();

                // 各時間ステップごとにグループを作成
                char group_name[32];
                snprintf(group_name, sizeof(group_name), "/data%06d", itime);
                group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                //sprintf(str, "group_name : %s", group_name);
                //fprintf(stdout, "%s\n", str);

                // Eフィールドデータセットの作成と書き込み
                hsize_t e_dims[4] = {1, NFreq2, NN, 6};
                //hsize_t e_dims[4] = {1, NFreq2, Nx*Ny*Nz, 6};
                dataspace_id = H5Screate_simple(4, e_dims, NULL);
                dataset_id = H5Dcreate(group_id, "E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims[1] = {6};
                memspace_id = H5Screate_simple(1, mem_dims, NULL);

                //sprintf(str, "NFreq2 : %d", NFreq2);
                //fprintf(stdout, "%s\n", str);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        //fprintf(stdout, "set e_value.\n");

                        double e_value[6] = {
                            cEx_r[n0 + nn], cEy_r[n0 + nn], cEz_r[n0 + nn],
                            cEx_i[n0 + nn], cEy_i[n0 + nn], cEz_i[n0 + nn]
                        };

                        hsize_t e_offset[4] = {0, ifreq, nn, 0};
                        hsize_t e_count[4] = {1, 1, 1, 6};
                        //fprintf(stdout, "H5Sselect_hyperslab.\n");
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, e_offset, NULL, e_count, NULL);

                        // 書き込み
                        //fprintf(stdout, "H5Dwrite.\n");
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, e_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing E data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }

                    /*
                    //Nx
                    int datacount = 0;
                    //for (int x = iMin + lx; x <= iMax - lx; x++)
                    for (int x = iMin; x < iMax; x++)
                    {
                        //Ny
                        //for (int y = jMin - ly; y <= jMax + ly; y++)
                        for (int y = jMin; y < jMax; y++)
                        {
                            //Ny
                            //for (int z = kMin + lz; z <= kMax - lz; z++)
                            for (int z = kMin; z < kMax; z++)
                            {
                                //境界線の対応は間引く
                                int nn = NA(x, y, z);
                                //printf("%d.\", nn)
                                //sprintf(str, "%7d", nn);
                                //fprintf(stdout, "%s\n", str);
                                double e_value[6] = {
                                    cEx_r[n0 + nn], cEy_r[n0 + nn], cEz_r[n0 + nn],
                                    cEx_i[n0 + nn], cEy_i[n0 + nn], cEz_i[n0 + nn]
                                };

                                hsize_t e_offset[4] = {0, ifreq, datacount, 0};
                                hsize_t e_count[4] = {1, 1, 1, 6};
                                //fprintf(stdout, "H5Sselect_hyperslab.\n");
                                H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, e_offset, NULL, e_count, NULL);

                                // 書き込み
                                //fprintf(stdout, "H5Dwrite.\n");
                                status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, e_value);
                                if (status < 0) {
                                    fprintf(stderr, "Error writing E data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                                }
                                datacount++;
                            }
                        }
                    }
                    */
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Hフィールドデータセットの作成と書き込み
                hsize_t h_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, h_dims, NULL);
                dataset_id = H5Dcreate(group_id, "H", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double h_value[6] = {
                            cHx_r[n0 + nn], cHy_r[n0 + nn], cHz_r[n0 + nn],
                            cHx_i[n0 + nn], cHy_i[n0 + nn], cHz_i[n0 + nn]
                        };

                        hsize_t h_offset[4] = {0, ifreq, nn, 0};
                        hsize_t h_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, h_offset, NULL, h_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, h_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // 複素数用のHDF5データ型を定義
                hid_t complex_datatype = H5Tcreate(H5T_COMPOUND, sizeof(d_complex_t));
                H5Tinsert(complex_datatype, "real", HOFFSET(d_complex_t, r), H5T_NATIVE_DOUBLE);
                H5Tinsert(complex_datatype, "imag", HOFFSET(d_complex_t, i), H5T_NATIVE_DOUBLE);

                // Surfaceフィールドデータセットの作成と書き込み
                hsize_t surf_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, surf_dims, NULL);
                dataset_id = H5Dcreate(group_id, "Surface", complex_datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    //int64_t surf0 = ifreq * NSurface;
                    for (int surf = 0; surf < NSurface; surf++) {
                        d_complex_t surf_value[6] = {
                            SurfaceEx[ifreq][surf], SurfaceEy[ifreq][surf], SurfaceEz[ifreq][surf],
                            SurfaceHx[ifreq][surf], SurfaceHy[ifreq][surf], SurfaceHz[ifreq][surf]
                        };

                        hsize_t surf_offset[4] = {0, ifreq, surf, 0};
                        hsize_t surf_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, surf_offset, NULL, surf_count, NULL);
                        status = H5Dwrite(dataset_id, complex_datatype, memspace_id, dataspace_id, H5P_DEFAULT, surf_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, surf=%d\n", itime, ifreq, surf);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Pフィールドデータセットの作成と書き込み（仮の例）
                hsize_t p_dims[4] = {1, NFreq2, NN, 3};
                dataspace_id = H5Screate_simple(4, p_dims, NULL);
                dataset_id = H5Dcreate(group_id, "P", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims2[1] = {3};
                memspace_id = H5Screate_simple(1, mem_dims2, NULL);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double p_value[3] = {
                            cEx_r[n0 + nn] * cHy_r[n0 + nn] - cEy_r[n0 + nn] * cHx_r[n0 + nn],
                            cEy_r[n0 + nn] * cEz_r[n0 + nn] - cEz_r[n0 + nn] * cHy_r[n0 + nn],
                            cEz_r[n0 + nn] * cHx_r[n0 + nn] - cEx_r[n0 + nn] * cEz_r[n0 + nn]
                        };

                        hsize_t p_offset[4] = {0, ifreq, nn, 0};
                        hsize_t p_count[4] = {1, 1, 1, 3};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, p_offset, NULL, p_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, p_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing P data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);
                // グループのクローズ
                H5Gclose(group_id);
            }

            // check convergence
            fmax[0] = MAX(fmax[0], fsum[0]);
            fmax[1] = MAX(fmax[1], fsum[1]);
            if ((fsum[0] < fmax[0] * Solver.converg) &&
                (fsum[1] < fmax[1] * Solver.converg)) {
                converged = 1;
                break;
            }
            
            // Niterを増加
            Niter++;
        }
    }

    // メモリスペース、データセットとデータスペースのクローズ
    status = H5Sclose(memspace_id);

    // result
    if (io) {
        sprintf(str, "    --- %s ---", (converged ? "converged" : "max steps"));
        fprintf(fp,     "%s\n", str);
        fprintf(stdout, "%s\n", str);
        fflush(fp);
        fflush(stdout);
    }

    // time steps
    Ntime = itime + converged;

    // copy point from device to host
    if (GPU) {
        copy_to_host();
    }

    // メタデータの作成
    hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //sprintf(str, "group_name : %s", group_name);
    fprintf(stdout, "meta1.\n");

    // 時間に関するメタデータの書き込み
    double time_metadata[1] = {Solver.maxiter * Dt};
    //dataspace_id = H5Screate_simple(1, count, NULL);
    hsize_t time_count[1] = {1};
    dataspace_id = H5Screate_simple(1, time_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_metadata);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // グリッドに関するメタデータの書き込み
    //double grid_metadata[3] = {Xn, Yn, Zn};
    //hsize_t grid_count[1] = {3};
    //dataspace_id = H5Screate_simple(1, grid_count, NULL);
    //dataset_id = H5Dcreate(metadata_group_id, "grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_metadata);
    //H5Dclose(dataset_id);
    //H5Sclose(dataspace_id);

    // Title
    fprintf(stdout, "meta2.\n");
    hsize_t title_dims[1] = {256};
    dataspace_id = H5Screate_simple(1, title_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    fprintf(stdout, "meta3.\n");
    // 各種整数型メタデータの書き込み
    struct {
        const char *name;
        void *value;
        hid_t type;
    } metadata[] = {
        {"Nx", &Nx, H5T_NATIVE_INT},
        {"Ny", &Ny, H5T_NATIVE_INT},
        {"Nz", &Nz, H5T_NATIVE_INT},
        {"Ni", &Ni, H5T_NATIVE_INT},
        {"Nj", &Nj, H5T_NATIVE_INT},
        {"Nk", &Nk, H5T_NATIVE_INT},
        {"N0", &N0, H5T_NATIVE_INT},
        {"NN", &NN, H5T_NATIVE_INT64},
        //境界線領域対応
        {"lx", &lx, H5T_NATIVE_INT},
        {"ly", &ly, H5T_NATIVE_INT},
        {"lz", &lz, H5T_NATIVE_INT},
        {"NFreq1", &NFreq1, H5T_NATIVE_INT},
        {"NFreq2", &NFreq2, H5T_NATIVE_INT},
        {"NFeed", &NFeed, H5T_NATIVE_INT},
        {"NPoint", &NPoint, H5T_NATIVE_INT},
        {"Niter", &Niter, H5T_NATIVE_INT},
        {"Ntime", &Ntime, H5T_NATIVE_INT},
        {"Solver_maxiter", &Solver.maxiter, H5T_NATIVE_INT},
        {"Solver_nout", &Solver.nout, H5T_NATIVE_INT},
        {"NGline", &NGline, H5T_NATIVE_INT},
        {"IPlanewave", &IPlanewave, H5T_NATIVE_INT}
    };

    for (int i = 0; i < sizeof(metadata) / sizeof(metadata[0]); i++) {
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, metadata[i].name, metadata[i].type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }

    fprintf(stdout, "meta4.\n");
    // Dtの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    fprintf(stdout, "meta5.\n");
    // Planewaveの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    fprintf(stdout, "meta6.\n");
    // 配列データの書き込み
    struct {
        const char *name;
        double *data;
        size_t size;
    } arrays[] = {
        {"Xn", Xn, Nx + 1},
        {"Yn", Yn, Ny + 1},
        {"Zn", Zn, Nz + 1},
        {"Xc", Xc, Nx},
        {"Yc", Yc, Ny},
        {"Zc", Zc, Nz},
        {"Eiter", Eiter, Niter},
        {"Hiter", Hiter, Niter},
        {"VFeed", VFeed, NFeed * (Solver.maxiter + 1)},
        {"IFeed", IFeed, NFeed * (Solver.maxiter + 1)},
        {"VPoint", VPoint, NPoint * (Solver.maxiter + 1)},
        {"Freq1", Freq1, NFreq1},
        {"Freq2", Freq2, NFreq2},
        //{"Gline", Gline, NGline * 2 * 3}
        {"Gline", reinterpret_cast<double*>(Gline), NGline * 2 * 3}
    };

    for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
        fprintf(stdout, "meta6 1(%d)(%s).\n", i, arrays[i].name);
        hsize_t array_dims[1] = {arrays[i].size};
        dataspace_id = H5Screate_simple(1, array_dims, NULL);
        dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i].data);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    
    fprintf(stdout, "meta7.\n");
    // Surfaceデータの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "NSurface", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NSurface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // surface_t構造体に対応する複合データ型を定義
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
    H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);

    hsize_t surface_dims[1] = {NSurface};
    dataspace_id = H5Screate_simple(1, surface_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Surface", memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(memtype);

    fprintf(stdout, "meta8.\n");
    // Refrection データの書き込み
    //dataspace_id = H5Screate(H5S_SCALAR);
    //dataset_id = H5Dcreate(metadata_group_id, "iEx", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &iEx);
    //H5Dclose(dataset_id);
    //H5Sclose(dataspace_id);

    // surface_t構造体に対応する複合データ型を定義
    /*
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
    H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);

    hsize_t surface_dims[1] = {NSurface};
    dataspace_id = H5Screate_simple(1, surface_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Surface", memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(memtype);
    */

    // 反射強度フィールドデータセットの作成と書き込み
    hsize_t p_dims[2] = {NN, 1};
    dataspace_id = H5Screate_simple(2, p_dims, NULL);
    //hsize_t p_dims[3] = {Nx+3, Ny+3, Nz+3};
    //hsize_t p_dims[4] = {Nx, Ny, Nz, 1};
    //dataspace_id = H5Screate_simple(4, p_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Reflection", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 書き込み用のメモリスペースを修正
    hsize_t mem_dims2[1] = {1};
    memspace_id = H5Screate_simple(1, mem_dims2, NULL);
    for (int nn = 0; nn < NN; nn++)
    {
        //double ref_value[1] = { iEx[nn] };
        double ref_value[1] = { sqrt(Material[iEx[nn]].epsr) };
        //fprintf(stdout, "ref_value = %f.\n", ref_value[0]);

        hsize_t ref_offset[2] = {nn, 0};
        hsize_t ref_count[2] = {1, 1};
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, ref_offset, NULL, ref_count, NULL);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, ref_value);
        if (status < 0) {
            fprintf(stderr, "Error writing Reflection data at nn=%d\n", nn);
        }
    }
    /*
    //Nx
    int datacount = 0;
    //for (int x = iMin + lx; x <= iMax - lx; x++)
    for (id_t x = iMin; x < iMax; x++)
    {
        //Ny
        //for (int y = jMin - ly; y <= jMax + ly; y++)
        for (id_t y = jMin; y < jMax; y++)
        {
            //Ny
            //for (int z = kMin + lz; z <= kMax - lz; z++)
            for (id_t z = kMin; z < kMax; z++)
            {
                //境界線の対応は間引く
                int nn = NA(x, y, z);
                //printf("%d.\", nn)
                //sprintf(str, "%7d", nn);
                //fprintf(stdout, "%s\n", str);
                //屈折率算出
                //double ref_value[1] = { sqrt(Material[iEx[datacount]].epsr) };
                double ref_value[1] = { sqrt(Material[iEx[nn]].epsr) };

                hsize_t ref_offset[4] = {x, y, z, 0};
                hsize_t ref_count[4] = {1, 1, 1, 1};
                //fprintf(stdout, "H5Sselect_hyperslab.\n");
                H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, ref_offset, NULL, ref_count, NULL);

                // 書き込み
                //fprintf(stdout, "H5Dwrite.\n");
                status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, ref_value);
                if (status < 0) {
                    fprintf(stderr, "Error writing Reflection data at nn=%d\n", nn);
                    fprintf(stderr, "Error writing Reflection data at datacount=%d\n", datacount);
                }
                datacount++;
            }
        }
    }
    status = H5Sclose(memspace_id);
    */
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // メタデータグループのクローズ
    H5Gclose(metadata_group_id);

    status = H5Fclose(file_id);

    // free
    memfree2_gpu();

    // copy near3d from device to host
    memcopy3_gpu();

    // free
    memfree3_gpu();
}


// copy from device to host
static void copy_to_host()
{
    if (NFeed) {
        cuda_memcpy(GPU, VFeed, d_VFeed, Feed_size, cudaMemcpyDeviceToHost);
        cuda_memcpy(GPU, IFeed, d_IFeed, Feed_size, cudaMemcpyDeviceToHost);
    }

    if (NPoint) {
        cuda_memcpy(GPU, VPoint, d_VPoint, Point_size, cudaMemcpyDeviceToHost);
    }
}


