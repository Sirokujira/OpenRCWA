/*
h5cells.c (OpenRCWA)

time_series_data.h5 の /dataNNNNNN に、セルごとの周波数領域の場
E / H / P / P_loss を書く。直列版 (sol/solve.c) と CUDA 版 (cuda/solve.cu)
で共通 (どちらもホスト側の全領域の配列を持つ)。MPI 版は rank ごとの
集団書き込みなので別実装 (mpi/solve.c, cuda_mpi/solve.cu)。

以前は 1 セルずつハイパースラブを選んで H5Dwrite していたため、出力 1 回
あたり NFreq2 * NN * 4 回の書き込みになり、pec_sphere_rcs.ofd (NN=14 万,
7 周波数) では計算時間の 85% を占めていた (96 秒 -> 15 秒)。ここでは
周波数ごとに 1 回で書く。データセットの形・並び・値は以前と同じ。
*/

#include <string.h>

#include "orcwa.h"
#include "h5cells.h"

/* [1, NFreq2, NN, ncomp] の double データセットを作り、周波数ごとに
   fill(ifreq, buf) で [nn][成分] の並びを詰めて 1 回ずつ書く。
   バッファは 1 周波数ぶん (NN * ncomp) だけで済む。 */
typedef void (*cell_fill_t)(int ifreq, double *buf);

static void write_cell_dataset(hid_t group_id, const char *name, int ncomp,
                               cell_fill_t fill, double *buf, int itime)
{
    hsize_t dims[4] = {1, (hsize_t)NFreq2, (hsize_t)NN, (hsize_t)ncomp};
    hid_t space = H5Screate_simple(4, dims, NULL);
    hid_t dset = H5Dcreate(group_id, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t mdims[1] = {(hsize_t)NN * ncomp};
    hid_t mspace = H5Screate_simple(1, mdims, NULL);
    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
        hsize_t offset[4] = {0, (hsize_t)ifreq, 0, 0};
        hsize_t count[4] = {1, 1, (hsize_t)NN, (hsize_t)ncomp};
        H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);
        fill(ifreq, buf);
        if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, space, H5P_DEFAULT, buf) < 0) {
            fprintf(stderr, "Error writing %s data at itime=%d, ifreq=%d\n", name, itime, ifreq);
        }
    }
    H5Sclose(mspace);
    H5Dclose(dset);
    H5Sclose(space);
}

static void fill_E(int ifreq, double *buf)
{
    const int64_t n0 = (int64_t)ifreq * NN;
    for (int64_t nn = 0; nn < NN; nn++) {
        double *e = &buf[nn * 6];
        e[0] = cEx_r[n0 + nn]; e[1] = cEy_r[n0 + nn]; e[2] = cEz_r[n0 + nn];
        e[3] = cEx_i[n0 + nn]; e[4] = cEy_i[n0 + nn]; e[5] = cEz_i[n0 + nn];
    }
}

static void fill_H(int ifreq, double *buf)
{
    const int64_t n0 = (int64_t)ifreq * NN;
    for (int64_t nn = 0; nn < NN; nn++) {
        double *h = &buf[nn * 6];
        h[0] = cHx_r[n0 + nn]; h[1] = cHy_r[n0 + nn]; h[2] = cHz_r[n0 + nn];
        h[3] = cHx_i[n0 + nn]; h[4] = cHy_i[n0 + nn]; h[5] = cHz_i[n0 + nn];
    }
}

// P (仮の例。MPI / CUDA 版と同じ式)
static void fill_P(int ifreq, double *buf)
{
    const int64_t n0 = (int64_t)ifreq * NN;
    for (int64_t nn = 0; nn < NN; nn++) {
        const int64_t m = n0 + nn;
        double *pv = &buf[nn * 3];
        pv[0] = cEx_r[m] * cHy_r[m] - cEy_r[m] * cHx_r[m];
        pv[1] = cEy_r[m] * cHz_r[m] - cEz_r[m] * cHy_r[m];
        pv[2] = cEz_r[m] * cHx_r[m] - cEx_r[m] * cHz_r[m];
    }
}

// 発熱量密度。fill_P_loss の前に P_loss_src を [NFreq2][NN] の配列へ向ける
static const double *P_loss_src;

static void fill_P_loss(int ifreq, double *buf)
{
    memcpy(buf, &P_loss_src[(int64_t)ifreq * NN], (size_t)NN * sizeof(double));
}

void write_cell_fields_hdf5(hid_t group_id, const double *P_losses, double *buf, int itime)
{
    write_cell_dataset(group_id, "E", 6, fill_E, buf, itime);
    write_cell_dataset(group_id, "H", 6, fill_H, buf, itime);
    write_cell_dataset(group_id, "P", 3, fill_P, buf, itime);
    P_loss_src = P_losses;
    write_cell_dataset(group_id, "P_loss", 1, fill_P_loss, buf, itime);
}
