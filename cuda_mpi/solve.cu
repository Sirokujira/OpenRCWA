/*
solve.cu (CUDA + MPI)
*/

#include "orcwa.h"
#include "orcwa_cuda.h"
#include "orcwa_prototype.h"

#include "hdf5.h"
#include <mpi.h>
#define FILE_NAME "time_series_data.h5"

static void setup_cuda_mpi();
static void copy_to_host();

void solve(int io, double *tdft, FILE *fp)
{
    // HDF5ファイルの作成
    // 関数から?(fp の入替え?)
    hid_t file_id;
    // local
    hid_t group_id, dataset_id, dataspace_id, memspace_id;
    herr_t status, ret;

    double fmax[] = {0, 0};
    char   str[BUFSIZ];
    int    converged = 0;

    // setup boundary index (MPI)
    setup_mpi();

    // setup host memory
    setup_host();

    // setup (GPU)
    if (GPU) {
        setup_gpu();
        setup_cuda_mpi();
    }

    // initial field
    initfield();

    // ── HDF5 (time_series_data.h5) は 1 プロセス実行のときだけ書く ──────
    // この実装の HDF5 出力は、集団操作 (H5Gcreate / H5Dcreate / H5Gclose) を
    // rank 0 だけが呼ぶ構造で、2 ランク以上では他ランクが H5Fclose / H5Fflush
    // (これも集団) で待ち続けてデッドロックする (実測: n=1 は通り n=2 でハング)。
    // 加えて寸法に rank 0 のローカル NN を使うので多ランクでは中身も正しくない。
    // 直すまでは多ランク時の HDF5 出力を止めて理由を出す (黙ってハングしない)。
    // 計算結果そのもの (.log / .out → *_post) は従来どおり出る。
    // 全ランクで同じ判定になるので、以下の h5on 分岐は集団操作の整合を崩さない。
    const int h5on = (commSize <= 1);
    hid_t plist_id = -1;
    if (h5on) {
        // MPIコミュニケータを使用したファイルアクセスプロパティリストの作成
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

        // HDF5ファイルの作成 (MPI対応)
        file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
    }
    else {
        file_id = -1;
        // 後段の出力関数 (sol/outputZin.c 等) は FILE_NAME を H5Fopen で開き直して
        // 結果を追記する。ファイルが無いので失敗するが、HDF5 の API は無効な id を
        // 受けても落ちずにエラーを返すだけ (結果は .log に出ている)。エラースタックの
        // 大量出力だけを止める (このプロセスでは以後 HDF5 を使わない)。
        H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
        if (commRank == 0) {
            fprintf(stderr, "*** note: %s is not written with %d processes "
                "(the parallel HDF5 output of this binary deadlocks with more than one rank); "
                "run with 1 process, or use the CPU / CUDA build for the time series\n",
                FILE_NAME, commSize);
            fflush(stderr);
        }
    }

    // time step iteration
    int itime;
    double t = 0;
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
            if (Npx > 1) {
                comm_cuda_X(1);
            }
            else {
                pbcx();
            }
        }
        if (PBCy) {
            if (Npy > 1) {
                comm_cuda_Y(1);
            }
            else {
                pbcy();
            }
        }
        if (PBCz) {
            if (Npz > 1) {
                comm_cuda_Z(1);
            }
            else {
                pbcz();
            }
        }

        // share boundary H (MPI)
        if (Npx > 1) {
            comm_cuda_X(0);
        }
        if (Npy > 1) {
            comm_cuda_Y(0);
        }
        if (Npz > 1) {
            comm_cuda_Z(0);
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
        const double t0 = comm_cputime();
        dftNear3d(itime);
        if (GPU) cudaDeviceSynchronize();
        *tdft += comm_cputime() - t0;

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // allreduce average (MPI)
            if (commSize > 1) {
                comm_average(fsum);
            }

            // average (post)
            if (commRank == 0) {
                Eiter[Niter] = fsum[0];
                Hiter[Niter] = fsum[1];
                Niter++;
            }

            // monitor
            if (io) {
                sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
                fprintf(fp,     "%s\n", str);
                fprintf(stdout, "%s\n", str);
                fflush(fp);
                fflush(stdout);
                
                // copy near3d from device to host
                memcopy3_gpu();

            }

            // HDF5 : 1 プロセス実行のときだけ (上の h5on を参照)
            if (io && h5on) {
                // グループの作成前に同期
                //MPI_Barrier(MPI_COMM_WORLD);

                // 各時間ステップごとにグループを作成
                char group_name[32];
                snprintf(group_name, sizeof(group_name), "/data%06d", itime);
                group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                //sprintf(str, "group_name : %s", group_name);
                //fprintf(stdout, "%s\n", str);

                // Eフィールドデータセットの作成と書き込み
                hsize_t e_dims[4] = {1, NFreq2, NN, 6};
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
                        // データ書き込み (MPI対応)
                        plist_id = H5Pcreate(H5P_DATASET_XFER);
                        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // H5FD_MPIO_COLLECTIVE または H5FD_MPIO_INDEPENDENT
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, plist_id, e_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing E data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                        H5Pclose(plist_id);
                    }
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
                            cHx_r[n0 + nn], cHy_r[n0 + nn], cEz_r[n0 + nn],
                            cHx_i[n0 + nn], cHy_i[n0 + nn], cHz_i[n0 + nn]
                        };

                        hsize_t h_offset[4] = {0, ifreq, nn, 0};
                        hsize_t h_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, h_offset, NULL, h_count, NULL);
                        
                        // データ書き込み (MPI対応)
                        plist_id = H5Pcreate(H5P_DATASET_XFER);
                        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, plist_id, h_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                        H5Pclose(plist_id);
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // グループのクローズ
                H5Gclose(group_id);
                
                // グループ作成後の同期
                //MPI_Barrier(MPI_COMM_WORLD);
                // メモリスペース、データセットとデータスペースのクローズ
                status = H5Sclose(memspace_id);
            }

            // check convergence
            fmax[0] = MAX(fmax[0], fsum[0]);
            fmax[1] = MAX(fmax[1], fsum[1]);
            if ((fsum[0] < fmax[0] * Solver.converg) &&
                (fsum[1] < fmax[1] * Solver.converg)) {
                converged = 1;
                break;
            }
        }
    }

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

    // グループの作成前に同期
    //MPI_Barrier(MPI_COMM_WORLD);

    if ((commRank == 0) && h5on) {
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
        // データ書き込み (MPI対応)
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, time_metadata);
        H5Pclose(plist_id);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

        // グリッドに関するメタデータの書き込み
        //double grid_metadata[3] = {Dx, Dy, Dz};
        //hsize_t grid_count[1] = {3};
        //dataspace_id = H5Screate_simple(1, grid_count, NULL);
        //dataset_id = H5Dcreate(metadata_group_id, "grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // データ書き込み (MPI対応)
        //plist_id = H5Pcreate(H5P_DATASET_XFER);
        //status = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        //if (status < 0) {
        //    fprintf(stderr, "Error setting MPI parameters for plist_id");
        //}
        //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, grid_metadata);
        //H5Pclose(plist_id);
        //H5Dclose(dataset_id);
        //H5Sclose(dataspace_id);

        // Title
        fprintf(stdout, "meta2.\n");
        hsize_t title_dims[1] = {256};
        dataspace_id = H5Screate_simple(1, title_dims, NULL);
        dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // データ書き込み (MPI対応)
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
        status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, Title);
        H5Pclose(plist_id);
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
            // データ書き込み (MPI対応)
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, plist_id, metadata[i].value);
            H5Pclose(plist_id);
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }

        fprintf(stdout, "meta4.\n");
        // Dtの書き込み
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, &Dt);
        H5Pclose(plist_id);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

        fprintf(stdout, "meta5.\n");
        // Planewaveの書き込み
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // データ書き込み (MPI対応)
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, &Planewave);
        H5Pclose(plist_id);
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
            {"Gline", reinterpret_cast<double*>(Gline), NGline * 2 * 3}
        };

        for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
            fprintf(stdout, "meta6 1(%d)(%s).\n", i, arrays[i].name);
            hsize_t array_dims[1] = {arrays[i].size};
            dataspace_id = H5Screate_simple(1, array_dims, NULL);
            dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // データ書き込み (MPI対応)
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, arrays[i].data);
            H5Pclose(plist_id);
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }
        
        fprintf(stdout, "meta7.\n");
        // Surfaceデータの書き込み
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, "NSurface", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // データ書き込み (MPI対応)
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
        status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, &NSurface);
        H5Pclose(plist_id);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

        fprintf(stdout, "meta8.\n");
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

        // データ書き込み (MPI対応)
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
        status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, plist_id, Surface);
        H5Pclose(plist_id);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

        fprintf(stdout, "meta9.\n");
        // メタデータグループのクローズ
        H5Gclose(metadata_group_id);
        
        fprintf(stdout, "meta10.\n");
    }

    // グループ作成後の同期
    MPI_Barrier(MPI_COMM_WORLD);

    if (h5on) {
        // キャッシュをフラッシュする
        status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
        if (status < 0) {
            fprintf(stderr, "Error H5Fflush\n");
        }

        //MPI 用に対応しているため
        status = H5Fclose(file_id);
        if (status < 0) {
            fprintf(stderr, "Error H5Fclose\n");
        }
    }

    // free
    memfree2_gpu();

    // copy near3d from device to host
    memcopy3_gpu();

    // free
    memfree3_gpu();

    // MPI : send to root
    if (commSize > 1) {
        // feed waveform
        if (NFeed) {
            comm_feed();
        }

        // point waveform
        if (NPoint) {
            comm_point();
        }

        // near3d
        if (NFreq2) {
            comm_near3d();
        }
    }
}


// setup
static void setup_cuda_mpi()
{
    size_t size;
    //printf("%d %d %d %d %d %d %d %d\n", commSize, commRank, bid.numhy_x, bid.numhz_x, bid.numhz_y, bid.numhx_y, bid.numhx_z, bid.numhy_z); fflush(stdout);

    // X boundary
    size = Bid.numhy_x * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hy_x, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hy_x, size);

    size = Bid.numhz_x * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hz_x, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hz_x, size);

    // Y boundary
    size = Bid.numhz_y * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hz_y, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hz_y, size);

    size = Bid.numhx_y * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hx_y, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hx_y, size);

    // Z boundary
    size = Bid.numhx_z * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hx_z, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hx_z, size);

    size = Bid.numhy_z * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hy_z, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hy_z, size);

    // block
    bufBlock = dim3(16, 16);

    // parameter
    //cudaMemcpyToSymbol(d_bid, &bid, sizeof(bid_t));
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

