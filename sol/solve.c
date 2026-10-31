#include "orcwa.h"
#include "complex.h"
#include "orcwa_prototype.h"
#include "ev.h"

#include "hdf5.h"
#include "h5cells.h"

#define FILE_NAME "time_series_data.h5"

/*
発熱量密度 P_loss の計算は sol/powerloss.c に移した (全ビルド構成で共有するため)。

ここには以前、次の 2 つがあった:
  - calculatePowerLoss(): 導電率を material_id = 0 (= 空気) 固定で引き、磁気損失を
    mu'' = 1e-3 のマジック定数にしていた。セルごとの材料を見ていなかったため、
    無損失問題 (PEC + 空気だけの dipole.ofd) でも最大 2.5e11 W/m^3 を返していた
    (正しくは全セル厳密に 0)。
  - updateTemperature(): 毎ステップ全セルに 7 点ステンシルをかけて温度 T を更新して
    いたが、T はどこにも出力されず free されるだけだった (結果を捨てる計算)。
    熱解析を入れるなら、出力経路とセルごとの材料定数から作り直すこと。
*/

void solve(int io, double *tdft, FILE *fp) {
    // HDF5ファイルの作成
        // 関数から?(fp の入替え?)
    hid_t file_id;
        // local
        hid_t group_id, dataset_id, dataspace_id;
    herr_t status;

    double fmax[] = {0, 0};
    char str[BUFSIZ];
    int converged = 0;

    // initial field
    initfield();

    /* 発熱量密度 [W/m^3]。出力時 (nout ごと) にだけ計算する。
       以前は毎ステップ全セルを計算していたが、書き出すのは nout ごとで、
       値は現在の DFT 配列だけで決まるので結果は変わらない。 */
    double *P_losses = (double *)malloc((size_t)NFreq2 * NN * sizeof(double));
    if (P_losses == NULL) {
        fprintf(stderr, "*** P_loss array malloc error (NFreq2=%d NN=%zu)\n", NFreq2, (size_t)NN);
        exit(1);
    }
    memset(P_losses, 0, (size_t)NFreq2 * NN * sizeof(double));

    /* HDF5 出力用の作業バッファ (1 周波数ぶん、E/H の 6 成分。P / P_loss も使い回す) */
    double *cellbuf = (double *)malloc((size_t)NN * 6 * sizeof(double));
    if (cellbuf == NULL) {
        fprintf(stderr, "*** HDF5 output buffer malloc error (NFreq2=%d NN=%zu)\n", NFreq2, (size_t)NN);
        exit(1);
    }

    // HDF5ファイルの作成
    file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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
            murH(numMurHx, fMurHx, Hx);
            murH(numMurHy, fMurHy, Hy);
            murH(numMurHz, fMurHz, Hz);
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
        const double t0 = cputime();
        dftNear3d(itime);
        *tdft += cputime() - t0;

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // average (post)
            Eiter[Niter] = fsum[0];
            Hiter[Niter] = fsum[1];
            Niter++;

            // monitor
            if (io) {
                sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
                fprintf(fp,     "%s\n", str);
                fprintf(stdout, "%s\n", str);
                fflush(fp);
                fflush(stdout);

                /* 発熱量密度は書き出す直前にだけ計算する
                   (値は現在の DFT 配列だけで決まるので毎ステップ回す必要がない) */
                calculatePowerLoss(P_losses);

                // 各時間ステップごとにグループを作成
                char group_name[32];
                snprintf(group_name, sizeof(group_name), "/data%06d", itime);
                group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // E/H/P/P_loss (sol/h5cells.c。CUDA 版と共通)
                write_cell_fields_hdf5(group_id, P_losses, cellbuf, itime);

                /* 近傍界サーフェス (遠方界計算用) はここには出さない。
                   NSurface / SurfaceEx 等は solve() 後の outputChars() で
                   初めて設定されるため、以前ここで作っていた /dataNNNNNN/Surface は
                   全ビルドで常に全要素 0 だった (しかも次元が NN)。並列 HDF5 では
                   領域が先に確保されるので、ゼロのまま MPI 版のファイルの 4 割を
                   占めていた。遠方界は orcwa.out (writeout) 経由で orcwa_post が
                   計算する。 */

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
            
        }
    }
    // メモリの解放
    free(P_losses);
    free(cellbuf);

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

    // メタデータの作成
    hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 時間に関するメタデータの書き込み(収束条件で終了時の対応)
    //int maxIter = MIN(Solver.maxiter, Ntime);
    //int maxNOut = MIN(Solver.nout, Niter);
    //double time_metadata[1] = {maxIter * Dt};
    double time_metadata[1] = {Solver.maxiter * Dt};
    //dataspace_id = H5Screate_simple(1, count, NULL);
    hsize_t time_count[1] = {1};
    dataspace_id = H5Screate_simple(1, time_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_metadata);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // グリッドに関するメタデータの書き込み
    //double grid_metadata[3] = {Dx, Dy, Dz};
    //hsize_t grid_count[1] = {3};
    //dataspace_id = H5Screate_simple(1, grid_count, NULL);
    //dataset_id = H5Dcreate(metadata_group_id, "grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_metadata);
    //H5Dclose(dataset_id);
    //H5Sclose(dataspace_id);

    // その他のメタデータの書き込み
/*
    // title, dt, source, fPlanewave, z0, Ni, Nj, Nk, N0, NN
    const char *title = Title;
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "title", H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 平面波が伝わるときのインピーダンスの値
    //double metadata_values[8] = {Dt, Planewave.z0, Ni, Nj, Nk, N0, NN};
    double metadata_values[8] = {Dt, 0.0, Ni, Nj, Nk, N0, NN};
    hsize_t metadata_count[1] = {8};
    dataspace_id = H5Screate_simple(1, metadata_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "metadata_values", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata_values);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 配列に関するメタデータの書き込み (Xn, Yn, Zn, Freq1, Freq2)
    hsize_t array_count[1];
    double *arrays[] = {Xn, Yn, Zn, Freq1, Freq2};
    const char *array_names[] = {"Xn", "Yn", "Zn", "Freq1", "Freq2"};
    size_t array_sizes[] = {Nx + 1, Ny + 1, Nz + 1, NFreq1, NFreq2};

    for (int i = 0; i < 5; i++) {
        array_count[0] = array_sizes[i];
        dataspace_id = H5Screate_simple(1, array_count, NULL);
        dataset_id = H5Dcreate(metadata_group_id, array_names[i], H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i]);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
*/
    // Title
    hsize_t title_dims[1] = {256};
    dataspace_id = H5Screate_simple(1, title_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

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
        status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }

    // Dtの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Planewaveの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

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
        {"Gline", (double *)Gline, NGline * 2 * 3}
    };

    for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
        hsize_t array_dims[1] = {arrays[i].size};
        dataspace_id = H5Screate_simple(1, array_dims, NULL);
        dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i].data);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    
    /* metadata/NSurface・metadata/Surface は出さない (per-step の Surface と
       同じ理由で、solve() の時点では NSurface = 0 の空データだった)。 */

    /* 屈折率マップ (metadata/Reflection)。セルごとの材料の sqrt(epsr)。
       CUDA 版だけが出していたため、CUDA ビルドに切り替えたときだけ GUI に
       出せるデータが増える (= 逆に言うと他のビルドでは欠ける) 状態だった。
       全ビルドで同じ内容を出す。 */
    {
        hsize_t ref_dims[2] = {(hsize_t)NN, 1};
        hid_t ref_space = H5Screate_simple(2, ref_dims, NULL);
        hid_t ref_set = H5Dcreate(metadata_group_id, "Reflection", H5T_NATIVE_DOUBLE,
                                  ref_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        double *ref = (double *)malloc((size_t)NN * sizeof(double));
        if (ref == NULL) {
            fprintf(stderr, "*** Reflection array malloc error (NN=%zu)\n", (size_t)NN);
            exit(1);
        }
        for (int64_t nn = 0; nn < NN; nn++) {
            ref[nn] = sqrt(Material[iEx[nn]].epsr);
        }
        status = H5Dwrite(ref_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref);
        if (status < 0) {
            fprintf(stderr, "Error writing Reflection data\n");
        }
        free(ref);
        H5Dclose(ref_set);
        H5Sclose(ref_space);
    }

    // メタデータグループのクローズ
    H5Gclose(metadata_group_id);

    status = H5Fclose(file_id);

    // free
    memfree2();
}

