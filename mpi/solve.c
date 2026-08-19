/*
solve.c (MPI)
*/

#include <stdlib.h>
#include "orcwa.h"
#include "orcwa_prototype.h"
#include "hdf5.h"
#include <mpi.h>
#define FILE_NAME "time_series_data.h5"

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

    // setup (MPI)
    setup_mpi();

    // initial field
    initfield();

    // MPIコミュニケータを使用したファイルアクセスプロパティリストの作成
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    /* メタデータ操作を集団化する。これを入れないと rank ごとにメタデータ
       キャッシュの状態がずれ、H5Fclose の内部で一部の rank が
       H5AC__run_sync_point の MPI_Barrier、残りが H5FD_truncate の
       MPI_File_set_size と別々の集団操作で待ち合ってハングする
       (実測: n<=3 では表面化せず n=4 で発生)。 */
    H5Pset_all_coll_metadata_ops(plist_id, 1);
    H5Pset_coll_metadata_write(plist_id, 1);

    // HDF5ファイルの作成 (MPI対応)
    file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /* ── 並列 HDF5 書き込み用のグローバル添字パラメータ ──────────────
       局所配列は袖 (halo/PML) を含む [iMin-l, iMax+l] を覆い、NA(i,j,k) で
       0..NN-1 に平坦化される (Nk=1 なので k 方向が連続)。HDF5 には直列実行時
       と同じ全体配列 (g_NN 要素) を書きたいので、setupSize と同じ式で全体側の
       Ni/Nj/Nk/N0 を求めておく (グローバル変数は書き換えない)。 */
    const int l_x = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
    const int l_y = l_x;
    const int l_z = l_x;
    const int64_t g_Nk = 1;
    const int64_t g_Nj = (int64_t)(Nz + 2 * l_z + 1);
    const int64_t g_Ni = (int64_t)(Ny + 2 * l_y + 1) * g_Nj;
    const int64_t g_N0 = -(((int64_t)(0 - l_x) * g_Ni)
                         + ((int64_t)(0 - l_y) * g_Nj)
                         + ((int64_t)(0 - l_z) * g_Nk));
    const int64_t g_NN = ((int64_t)(Nx + l_x) * g_Ni)
                       + ((int64_t)(Ny + l_y) * g_Nj)
                       + ((int64_t)(Nz + l_z) * g_Nk) + g_N0 + 1;

    /* 各 rank の担当範囲 (グローバル添字, 両端含む)。
       重複書き込みは MPI-IO では未定義動作なので、内部 rank は自分のコア
       [iMin, iMax-1] のみ、領域端の rank だけ袖まで広げる。こうすると全体を
       過不足なくちょうど 1 回ずつ覆う。 */
    const int w_i0 = (Ipx == 0)       ? (0  - l_x) : iMin;
    const int w_i1 = (Ipx == Npx - 1) ? (Nx + l_x) : (iMax - 1);
    const int w_j0 = (Ipy == 0)       ? (0  - l_y) : jMin;
    const int w_j1 = (Ipy == Npy - 1) ? (Ny + l_y) : (jMax - 1);
    const int w_k0 = (Ipz == 0)       ? (0  - l_z) : kMin;
    const int w_k1 = (Ipz == Npz - 1) ? (Nz + l_z) : (kMax - 1);

    /* 担当範囲のセル数 (集団書き込みのバッファ長に使う) */
    const int64_t w_ni    = (int64_t)(w_i1 - w_i0 + 1);
    const int64_t w_nj    = (int64_t)(w_j1 - w_j0 + 1);
    const int64_t w_nk    = (int64_t)(w_k1 - w_k0 + 1);
    const int64_t w_ncell = w_ni * w_nj * w_nk;
    /* y/z 方向が分割されていない (担当が全範囲) なら、担当セルは平坦添字上で
       連続になる。既定の分割は x のみ (Npy=Npz=1) なので通常こちらを通る。 */
    const int w_full_jk = ((w_j0 == (0 - l_y)) && (w_j1 == (Ny + l_y))
                        && (w_k0 == (0 - l_z)) && (w_k1 == (Nz + l_z)));

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
            if (Npx > 1) {
                comm_X(1);
            }
            else {
                pbcx();
            }
        }
        if (PBCy) {
            if (Npy > 1) {
                comm_Y(1);
            }
            else {
                pbcy();
            }
        }
        if (PBCz) {
            if (Npz > 1) {
                comm_Z(1);
            }
            else {
                pbcz();
            }
        }

        // share boundary H (MPI)
        if (Npx > 1) {
            comm_X(0);
        }
        if (Npy > 1) {
            comm_Y(0);
        }
        if (Npz > 1) {
            comm_Z(0);
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
        const double t0 = comm_cputime();
        dftNear3d(itime);
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
            }

            /* HDF5 出力は全 rank が実行する。並列 HDF5 (H5Pset_fapl_mpio で
               作成したファイル) では H5Gcreate / H5Dcreate は集団操作であり、
               rank 0 だけで呼ぶと他 rank が H5Fclose (これも集団) で待ち続けて
               デッドロックする (実測: n=1 は通り n>=2 でハング)。
               書き込む中身は各 rank が自分の担当範囲だけをグローバル添字で置く。 */
            {
                // グループの作成前に同期
                //MPI_Barrier(MPI_COMM_WORLD);
                // 各時間ステップごとにグループを作成
                char group_name[32];
                snprintf(group_name, sizeof(group_name), "/data%06d", itime);
                group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // Eフィールドデータセットの作成と書き込み
                hsize_t e_dims[4] = {1, NFreq2, g_NN, 6};
                dataspace_id = H5Screate_simple(4, e_dims, NULL);
                dataset_id = H5Dcreate(group_id, "E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                /* 集団書き込み: 全 rank が同じ回数だけ H5Dwrite を呼ぶ。
                   1 セルずつ独立書き込みしていた旧実装は、不均等分割だと
                   rank ごとに呼び出し回数が変わり、HDF5 内部のメタデータ
                   キャッシュ同期がずれて H5Fclose でハングした
                   (実測: 均等分割の n=2,3,5,6 は通り n=4,7 でハング)。
                   自 rank の担当セルは (i,j) ごとに k 方向が連続なので、
                   その run を OR で足し合わせて 1 回で書く。 */
                {
                    const int64_t nsel = w_ncell * 6;
                    double *buf = (double *)malloc((size_t)nsel * sizeof(double));
                    hsize_t mdims[1] = {(hsize_t)nsel};
                    hid_t mspace = H5Screate_simple(1, mdims, NULL);
                    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

                    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                        int64_t n0 = ifreq * NN;

                        /* ファイル側の選択 (自 rank の担当セル)。
                           y/z を分割していない (既定の x 方向のみ分割) 場合、
                           担当セルは平坦添字上で単一の連続ブロックになるので
                           1 個のハイパースラブで表せる。多数のブロックを OR で
                           繋いだ不規則な選択にすると、HDF5 内部の集団 I/O 判定が
                           rank ごとに分岐しうるため、可能な限り単純にする。 */
                        H5Sselect_none(dataspace_id);
                        if (w_full_jk) {
                            const int64_t g_beg = ((int64_t)w_i0 * g_Ni) + ((int64_t)w_j0 * g_Nj)
                                                + ((int64_t)w_k0 * g_Nk) + g_N0;
                            hsize_t st[4] = {0, (hsize_t)ifreq, (hsize_t)g_beg, 0};
                            hsize_t ct[4] = {1, 1, (hsize_t)(w_ni * g_Ni), 6};
                            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, st, NULL, ct, NULL);
                        } else {
                        for (int gi = w_i0; gi <= w_i1; gi++) {
                        for (int gj = w_j0; gj <= w_j1; gj++) {
                            const int64_t g_run = ((int64_t)gi * g_Ni) + ((int64_t)gj * g_Nj)
                                                + ((int64_t)w_k0 * g_Nk) + g_N0;
                            hsize_t st[4] = {0, (hsize_t)ifreq, (hsize_t)g_run, 0};
                            hsize_t ct[4] = {1, 1, (hsize_t)w_nk, 6};
                            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_OR, st, NULL, ct, NULL);
                        }
                        }
                        }

                        /* メモリ側は選択と同じ (i,j,k) 順に詰める */
                        int64_t q = 0;
                        for (int gi = w_i0; gi <= w_i1; gi++) {
                        for (int gj = w_j0; gj <= w_j1; gj++) {
                        for (int gk = w_k0; gk <= w_k1; gk++) {
                            const int64_t nn = ((int64_t)gi * Ni) + ((int64_t)gj * Nj)
                                             + ((int64_t)gk * Nk) + N0;
                            const double v[6] = {
                            cEx_r[n0 + nn], cEy_r[n0 + nn], cEz_r[n0 + nn],
                            cEx_i[n0 + nn], cEy_i[n0 + nn], cEz_i[n0 + nn]
                            };
                            for (int c = 0; c < 6; c++) buf[q++] = v[c];
                        }
                        }
                        }

                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, mspace, dataspace_id, dxpl, buf);
                        if (status < 0) {
                            fprintf(stderr, "Error writing E data at itime=%d, ifreq=%d\n", itime, ifreq);
                        }
                    }
                    H5Pclose(dxpl);
                    H5Sclose(mspace);
                    free(buf);
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Hフィールドデータセットの作成と書き込み
                hsize_t h_dims[4] = {1, NFreq2, g_NN, 6};
                dataspace_id = H5Screate_simple(4, h_dims, NULL);
                dataset_id = H5Dcreate(group_id, "H", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                /* 集団書き込み: 全 rank が同じ回数だけ H5Dwrite を呼ぶ。
                   1 セルずつ独立書き込みしていた旧実装は、不均等分割だと
                   rank ごとに呼び出し回数が変わり、HDF5 内部のメタデータ
                   キャッシュ同期がずれて H5Fclose でハングした
                   (実測: 均等分割の n=2,3,5,6 は通り n=4,7 でハング)。
                   自 rank の担当セルは (i,j) ごとに k 方向が連続なので、
                   その run を OR で足し合わせて 1 回で書く。 */
                {
                    const int64_t nsel = w_ncell * 6;
                    double *buf = (double *)malloc((size_t)nsel * sizeof(double));
                    hsize_t mdims[1] = {(hsize_t)nsel};
                    hid_t mspace = H5Screate_simple(1, mdims, NULL);
                    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

                    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                        int64_t n0 = ifreq * NN;

                        /* ファイル側の選択 (自 rank の担当セル)。
                           y/z を分割していない (既定の x 方向のみ分割) 場合、
                           担当セルは平坦添字上で単一の連続ブロックになるので
                           1 個のハイパースラブで表せる。多数のブロックを OR で
                           繋いだ不規則な選択にすると、HDF5 内部の集団 I/O 判定が
                           rank ごとに分岐しうるため、可能な限り単純にする。 */
                        H5Sselect_none(dataspace_id);
                        if (w_full_jk) {
                            const int64_t g_beg = ((int64_t)w_i0 * g_Ni) + ((int64_t)w_j0 * g_Nj)
                                                + ((int64_t)w_k0 * g_Nk) + g_N0;
                            hsize_t st[4] = {0, (hsize_t)ifreq, (hsize_t)g_beg, 0};
                            hsize_t ct[4] = {1, 1, (hsize_t)(w_ni * g_Ni), 6};
                            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, st, NULL, ct, NULL);
                        } else {
                        for (int gi = w_i0; gi <= w_i1; gi++) {
                        for (int gj = w_j0; gj <= w_j1; gj++) {
                            const int64_t g_run = ((int64_t)gi * g_Ni) + ((int64_t)gj * g_Nj)
                                                + ((int64_t)w_k0 * g_Nk) + g_N0;
                            hsize_t st[4] = {0, (hsize_t)ifreq, (hsize_t)g_run, 0};
                            hsize_t ct[4] = {1, 1, (hsize_t)w_nk, 6};
                            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_OR, st, NULL, ct, NULL);
                        }
                        }
                        }

                        /* メモリ側は選択と同じ (i,j,k) 順に詰める */
                        int64_t q = 0;
                        for (int gi = w_i0; gi <= w_i1; gi++) {
                        for (int gj = w_j0; gj <= w_j1; gj++) {
                        for (int gk = w_k0; gk <= w_k1; gk++) {
                            const int64_t nn = ((int64_t)gi * Ni) + ((int64_t)gj * Nj)
                                             + ((int64_t)gk * Nk) + N0;
                            const double v[6] = {
                            cHx_r[n0 + nn], cHy_r[n0 + nn], cHz_r[n0 + nn],
                            cHx_i[n0 + nn], cHy_i[n0 + nn], cHz_i[n0 + nn]
                            };
                            for (int c = 0; c < 6; c++) buf[q++] = v[c];
                        }
                        }
                        }

                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, mspace, dataspace_id, dxpl, buf);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d\n", itime, ifreq);
                        }
                    }
                    H5Pclose(dxpl);
                    H5Sclose(mspace);
                    free(buf);
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // 複素数用のHDF5データ型を定義
                hid_t complex_datatype = H5Tcreate(H5T_COMPOUND, sizeof(d_complex_t));
                H5Tinsert(complex_datatype, "real", HOFFSET(d_complex_t, r), H5T_NATIVE_DOUBLE);
                H5Tinsert(complex_datatype, "imag", HOFFSET(d_complex_t, i), H5T_NATIVE_DOUBLE);

                // Surfaceフィールドデータセットの作成と書き込み
                /* 第 3 次元は E/H/P と同じ「全体配列の添字空間」= g_NN。
                   ここを rank ローカルの NN にしていたため、不均等分割では
                   rank ごとに違う次元で H5Dcreate (集団操作) を呼ぶことになり、
                   HDF5 内部のメタデータキャッシュ同期がずれて H5Fclose で
                   デッドロックしていた (均等分割の n=2,3,5,6 は NN が全 rank
                   で一致するため顕在化しない。n=4,7 でハング)。 */
                hsize_t surf_dims[4] = {1, NFreq2, g_NN, 6};
                dataspace_id = H5Screate_simple(4, surf_dims, NULL);
                dataset_id = H5Dcreate(group_id, "Surface", complex_datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                /* Surface の書き込み用メモリスペース (1 要素あたり 6 成分)。
                   E/H/P を集団書き込みに作り替えた際に、共有していた
                   memspace_id の生成が消えて H5Sclose が
                   "not a dataspace" で失敗していた。ここで作り直す。 */
                hsize_t surf_mem_dims[1] = {6};
                memspace_id = H5Screate_simple(1, surf_mem_dims, NULL);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    //int64_t surf0 = ifreq * NSurface;
                    /* Surface は NSurface (近傍界サーフェス要素) 基準で、
                       セル分割との対応が未確認のため並列化していない。
                       H5Dcreate は上で全 rank が呼んでおり (集団操作)、
                       ここは書き込みだけを rank 0 に限定する (転送は既定 =
                       INDEPENDENT なのでデッドロックしない)。
                       TODO: SurfaceEx 等の分割を確認して並列化する。 */
                    for (int surf = 0; (commRank == 0) && (surf < NSurface); surf++) {
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
                hsize_t p_dims[4] = {1, NFreq2, g_NN, 3};
                dataspace_id = H5Screate_simple(4, p_dims, NULL);
                dataset_id = H5Dcreate(group_id, "P", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                /* 集団書き込み: 全 rank が同じ回数だけ H5Dwrite を呼ぶ。
                   1 セルずつ独立書き込みしていた旧実装は、不均等分割だと
                   rank ごとに呼び出し回数が変わり、HDF5 内部のメタデータ
                   キャッシュ同期がずれて H5Fclose でハングした
                   (実測: 均等分割の n=2,3,5,6 は通り n=4,7 でハング)。
                   自 rank の担当セルは (i,j) ごとに k 方向が連続なので、
                   その run を OR で足し合わせて 1 回で書く。 */
                {
                    const int64_t nsel = w_ncell * 3;
                    double *buf = (double *)malloc((size_t)nsel * sizeof(double));
                    hsize_t mdims[1] = {(hsize_t)nsel};
                    hid_t mspace = H5Screate_simple(1, mdims, NULL);
                    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

                    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                        int64_t n0 = ifreq * NN;

                        /* ファイル側の選択 (自 rank の担当セル)。
                           y/z を分割していない (既定の x 方向のみ分割) 場合、
                           担当セルは平坦添字上で単一の連続ブロックになるので
                           1 個のハイパースラブで表せる。多数のブロックを OR で
                           繋いだ不規則な選択にすると、HDF5 内部の集団 I/O 判定が
                           rank ごとに分岐しうるため、可能な限り単純にする。 */
                        H5Sselect_none(dataspace_id);
                        if (w_full_jk) {
                            const int64_t g_beg = ((int64_t)w_i0 * g_Ni) + ((int64_t)w_j0 * g_Nj)
                                                + ((int64_t)w_k0 * g_Nk) + g_N0;
                            hsize_t st[4] = {0, (hsize_t)ifreq, (hsize_t)g_beg, 0};
                            hsize_t ct[4] = {1, 1, (hsize_t)(w_ni * g_Ni), 3};
                            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, st, NULL, ct, NULL);
                        } else {
                        for (int gi = w_i0; gi <= w_i1; gi++) {
                        for (int gj = w_j0; gj <= w_j1; gj++) {
                            const int64_t g_run = ((int64_t)gi * g_Ni) + ((int64_t)gj * g_Nj)
                                                + ((int64_t)w_k0 * g_Nk) + g_N0;
                            hsize_t st[4] = {0, (hsize_t)ifreq, (hsize_t)g_run, 0};
                            hsize_t ct[4] = {1, 1, (hsize_t)w_nk, 3};
                            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_OR, st, NULL, ct, NULL);
                        }
                        }
                        }

                        /* メモリ側は選択と同じ (i,j,k) 順に詰める */
                        int64_t q = 0;
                        for (int gi = w_i0; gi <= w_i1; gi++) {
                        for (int gj = w_j0; gj <= w_j1; gj++) {
                        for (int gk = w_k0; gk <= w_k1; gk++) {
                            const int64_t nn = ((int64_t)gi * Ni) + ((int64_t)gj * Nj)
                                             + ((int64_t)gk * Nk) + N0;
                            const double v[3] = {
                            cEx_r[n0 + nn] * cHy_r[n0 + nn] - cEy_r[n0 + nn] * cHx_r[n0 + nn],
                            cEy_r[n0 + nn] * cHz_r[n0 + nn] - cEz_r[n0 + nn] * cHy_r[n0 + nn],
                            cEz_r[n0 + nn] * cHx_r[n0 + nn] - cEx_r[n0 + nn] * cHz_r[n0 + nn]
                            };
                            for (int c = 0; c < 3; c++) buf[q++] = v[c];
                        }
                        }
                        }

                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, mspace, dataspace_id, dxpl, buf);
                        if (status < 0) {
                            fprintf(stderr, "Error writing P data at itime=%d, ifreq=%d\n", itime, ifreq);
                        }
                    }
                    H5Pclose(dxpl);
                    H5Sclose(mspace);
                    free(buf);
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

    // グループの作成前に同期
    MPI_Barrier(MPI_COMM_WORLD);

    /* 集団 H5Dcreate は全 rank でデータ空間の次元が一致していなければ
       ならない。Niter は rank 0 でしか加算されず (average の後段)、
       NSurface も rank ごとに異なりうるため、そのままだと Eiter/Hiter 等の
       次元が食い違って集団呼び出しが噛み合わずハングする。
       メタデータを書く前に rank 0 の値へ揃えておく。 */
#ifdef _MPI
    if (commSize > 1) {
        MPI_Bcast(&Niter,    1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&NSurface, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Ntime,    1, MPI_INT, 0, MPI_COMM_WORLD);

        /* 給電点・観測点の波形は、そのセルを担当する rank にしか溜まらない。
           rank 0 が担当していない分割では IFeed/VFeed/VPoint が空のまま
           出力されるため (実測: n=2 は rank 0 が給電セルを持つので一致、
           n=4/7 では IFeed が直列と 6.6e-3 ずれる)、ここで rank 0 に集める。
           comm_feed/comm_point は用意されていたが呼ばれていなかった。
           内部で MPI_Barrier を使うので全 rank が通ること。 */
        comm_feed();
        comm_point();
    }
#endif

    /* 並列 HDF5 では H5Gcreate / H5Dcreate は集団操作であり、全 rank が
       同じ順序で呼ぶ必要がある。ここを if (commRank == 0) で囲うと rank 0 が
       集団呼び出しに入る一方、他 rank は後段の H5Fclose (これも集団) で待ち
       続けデッドロックする (実測: n=1 は通り n>=2 でハング)。
       メタデータは全 rank で同一なので、作成・クローズは全 rank で行い、
       実データの書き込みだけを INDEPENDENT 転送で rank 0 に限定する。 */
    {
        // メタデータの作成
        hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        //sprintf(str, "group_name : %s", group_name);
        //fprintf(stdout, "meta1.\n");

        // 時間に関するメタデータの書き込み
        double time_metadata[1] = {Solver.maxiter * Dt};
        //dataspace_id = H5Screate_simple(1, count, NULL);
        hsize_t time_count[1] = {1};
        dataspace_id = H5Screate_simple(1, time_count, NULL);
        dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // データ書き込み (MPI対応)
        if (commRank == 0) {   /* 書き込みは rank 0 のみ */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, time_metadata);
            H5Pclose(plist_id);
        }
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
        hsize_t title_dims[1] = {256};
        dataspace_id = H5Screate_simple(1, title_dims, NULL);
        dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // データ書き込み (MPI対応)
        if (commRank == 0) {   /* 書き込みは rank 0 のみ */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, Title);
            H5Pclose(plist_id);
        }
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

        // 各種整数型メタデータの書き込み
        /* metadata[] は void* を取るので、全体側の値を同じ型の実体に置く。
           E/H/P を全体配列として書いている以上、その形状パラメータも
           全体側でなければファイルの自己整合性が壊れる。 */
        const int     m_Ni = (int)g_Ni;
        const int     m_Nj = (int)g_Nj;
        const int     m_Nk = (int)g_Nk;
        const int     m_N0 = (int)g_N0;
        const int64_t m_NN = g_NN;

        struct {
            const char *name;
            void *value;
            hid_t type;
        } metadata[] = {
            {"Nx", &Nx, H5T_NATIVE_INT},
            {"Ny", &Ny, H5T_NATIVE_INT},
            {"Nz", &Nz, H5T_NATIVE_INT},
            /* E/H/P は全体配列 (g_NN 要素) として書いているので、その添字を
               解釈するためのパラメータも全体側の値を書く。局所値 (Ni/Nj/...) を
               書くと rank 0 の部分領域の形状になり、ファイルの自己整合性が
               壊れる (直列実行時の出力と解釈が食い違う)。 */
            {"Ni", &m_Ni, H5T_NATIVE_INT},
            {"Nj", &m_Nj, H5T_NATIVE_INT},
            {"Nk", &m_Nk, H5T_NATIVE_INT},
            {"N0", &m_N0, H5T_NATIVE_INT},
            {"NN", &m_NN, H5T_NATIVE_INT64},
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
            if (commRank == 0) {   /* 書き込みは rank 0 のみ */
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
                status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, plist_id, metadata[i].value);
                H5Pclose(plist_id);
            }
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }

        // Dtの書き込み
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (commRank == 0) {   /* 書き込みは rank 0 のみ */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, &Dt);
            H5Pclose(plist_id);
        }
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

        // Planewaveの書き込み
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // データ書き込み (MPI対応)
        if (commRank == 0) {   /* 書き込みは rank 0 のみ */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, &Planewave);
            H5Pclose(plist_id);
        }
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
            {"Gline", Gline, NGline * 2 * 3}
        };

        /* 集団操作である H5Dcreate は全 rank で次元が一致していなければ
           ならない。NGline は rank 0 でしか計算されず (Gline は
           comm_broadcast の対象外)、そのままだと rank ごとに違う次元で
           H5Dcreate を呼ぶことになり、データセットが 0 要素になったり
           H5Fclose でデッドロックしたりする。書き込むのは rank 0 だけ
           なので、次元だけ rank 0 の値へ揃える。 */
#ifdef _MPI
        if (commSize > 1) {
            const int narray = (int)(sizeof(arrays) / sizeof(arrays[0]));
            int64_t asize[sizeof(arrays) / sizeof(arrays[0])];
            for (int i = 0; i < narray; i++) asize[i] = (int64_t)arrays[i].size;
            MPI_Bcast(asize, narray, MPI_INT64_T, 0, MPI_COMM_WORLD);
            for (int i = 0; i < narray; i++) arrays[i].size = (size_t)asize[i];
        }
#endif

        for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
            hsize_t array_dims[1] = {arrays[i].size};
            dataspace_id = H5Screate_simple(1, array_dims, NULL);
            dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // データ書き込み (MPI対応)
            if (commRank == 0) {   /* 書き込みは rank 0 のみ */
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
                status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, arrays[i].data);
                H5Pclose(plist_id);
            }
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }
        
        // Surfaceデータの書き込み
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, "NSurface", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // データ書き込み (MPI対応)
        if (commRank == 0) {   /* 書き込みは rank 0 のみ */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, &NSurface);
            H5Pclose(plist_id);
        }
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);

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
        if (commRank == 0) {   /* 書き込みは rank 0 のみ */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // または H5FD_MPIO_INDEPENDENT
            status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, plist_id, Surface);
            H5Pclose(plist_id);
        }
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
	    H5Tclose(memtype);

        // メタデータグループのクローズ
        H5Gclose(metadata_group_id);
        

        // キャッシュをフラッシュする
        //status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
        //if (status < 0) {
        //    fprintf(stderr, "Error H5Fflush\n");
        //}
    }

	//MPI 用に対応しているためプロセス毎に対応(並列実行用の対応必要?[そのままだと実行時警告?エラー?出力])
    status = H5Fclose(file_id);
    if (status < 0) {
        fprintf(stderr, "Error H5Fclose\n");
    }

    // free
    memfree2();

    // グループ作成後の同期
    MPI_Barrier(MPI_COMM_WORLD);

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

        // near3d (index is changed)
        if (NFreq2) {
            comm_near3d();
        }
    }
}

