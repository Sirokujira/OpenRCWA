/*
rcwa_hdf5.c

RCWA スペクトルの HDF5 出力。レイアウトは include/rcwa_hdf5.h を参照。
既存の HDF5 出力 (sol/solve.c, sol/outputZin.c 等) と同じ C API の使い方に
そろえてある。
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rcwa_hdf5.h"
#include "hdf5.h"

/* スカラーの int / double を 1 つ書く小物 */
static void write_scalar_int(hid_t gid, const char *name, int value)
{
	hid_t sid = H5Screate(H5S_SCALAR);
	hid_t did = H5Dcreate(gid, name, H5T_NATIVE_INT, sid,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (did >= 0) {
		H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
		H5Dclose(did);
	}
	H5Sclose(sid);
}

static void write_scalar_double(hid_t gid, const char *name, double value)
{
	hid_t sid = H5Screate(H5S_SCALAR);
	hid_t did = H5Dcreate(gid, name, H5T_NATIVE_DOUBLE, sid,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (did >= 0) {
		H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
		H5Dclose(did);
	}
	H5Sclose(sid);
}

/* 可変長でない固定長文字列を 1 つ書く */
static void write_string(hid_t gid, const char *name, const char *value, size_t len)
{
	hid_t tid = H5Tcopy(H5T_C_S1);
	H5Tset_size(tid, len);
	H5Tset_strpad(tid, H5T_STR_NULLTERM);

	char *buf = (char *)calloc(len, 1);
	if (buf == NULL) { H5Tclose(tid); return; }
	strncpy(buf, (value != NULL) ? value : "", len - 1);

	hid_t sid = H5Screate(H5S_SCALAR);
	hid_t did = H5Dcreate(gid, name, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (did >= 0) {
		H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
		H5Dclose(did);
	}
	H5Sclose(sid);
	H5Tclose(tid);
	free(buf);
}

/* 偏波ラベルを固定長文字列の 1 次元配列として書く */
static void write_pol_labels(hid_t gid, int npol, const char *const *labels)
{
	hid_t tid = H5Tcopy(H5T_C_S1);
	H5Tset_size(tid, RCWA_POL_LABEL_LEN);
	H5Tset_strpad(tid, H5T_STR_NULLTERM);

	char *buf = (char *)calloc((size_t)npol, RCWA_POL_LABEL_LEN);
	if (buf == NULL) { H5Tclose(tid); return; }
	for (int i = 0; i < npol; i++) {
		const char *src = (labels[i] != NULL) ? labels[i] : "";
		strncpy(buf + (size_t)i * RCWA_POL_LABEL_LEN, src, RCWA_POL_LABEL_LEN - 1);
	}

	hsize_t dims[1];
	dims[0] = (hsize_t)npol;
	hid_t sid = H5Screate_simple(1, dims, NULL);
	hid_t did = H5Dcreate(gid, "pol_labels", tid, sid,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (did >= 0) {
		H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
		H5Dclose(did);
	}
	H5Sclose(sid);
	H5Tclose(tid);
	free(buf);
}

int rcwa_write_hdf5(const char *filename,
                    const rcwa_meta_t *meta,
                    const rcwa_spectrum_row_t *rows,
                    int npol,
                    int nfreq,
                    const char *const *pol_labels)
{
	if ((filename == NULL) || (meta == NULL) || (rows == NULL) ||
	    (npol <= 0) || (nfreq <= 0) || (pol_labels == NULL)) {
		return 1;
	}

	/* RCWA モードでは FDTD の solve() を通らないため、ここでファイルを作る */
	hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "*** file %s create error (HDF5).\n", filename);
		return 1;
	}

	/* ---- /metadata ---- */
	hid_t meta_gid = H5Gcreate(file_id, "/metadata",
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (meta_gid < 0) {
		H5Fclose(file_id);
		return 1;
	}
	/* FDTD 出力と区別するための目印 (GUI 側はこれの有無で判別する) */
	write_string(meta_gid, "solver_mode", "RCWA", 8);
	write_string(meta_gid, "Title", meta->title, 256);
	write_scalar_int(meta_gid, "NFreq", nfreq);
	write_scalar_int(meta_gid, "npol", npol);
	write_pol_labels(meta_gid, npol, pol_labels);
	write_scalar_int(meta_gid, "rcwa_harmonics", meta->nharmonics);
	write_scalar_double(meta_gid, "rcwa_period", meta->period);
	write_scalar_int(meta_gid, "rcwa_harmonics_y", meta->nharmonics_y);
	write_scalar_double(meta_gid, "rcwa_period_y", meta->period_y);
	write_scalar_int(meta_gid, "rcwa_nlayer", meta->nlayer);
	write_scalar_double(meta_gid, "theta", meta->theta);
	write_scalar_double(meta_gid, "phi", meta->phi);
	H5Gclose(meta_gid);

	/* ---- /rcwa/spectrum ---- */
	hid_t rcwa_gid = H5Gcreate(file_id, "/rcwa",
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (rcwa_gid < 0) {
		H5Fclose(file_id);
		return 1;
	}

	hid_t tid = H5Tcreate(H5T_COMPOUND, sizeof(rcwa_spectrum_row_t));
	H5Tinsert(tid, "frequency", HOFFSET(rcwa_spectrum_row_t, frequency), H5T_NATIVE_DOUBLE);
	H5Tinsert(tid, "lambda",    HOFFSET(rcwa_spectrum_row_t, lambda),    H5T_NATIVE_DOUBLE);
	H5Tinsert(tid, "R",         HOFFSET(rcwa_spectrum_row_t, R),         H5T_NATIVE_DOUBLE);
	H5Tinsert(tid, "T",         HOFFSET(rcwa_spectrum_row_t, T),         H5T_NATIVE_DOUBLE);
	H5Tinsert(tid, "A",         HOFFSET(rcwa_spectrum_row_t, A),         H5T_NATIVE_DOUBLE);

	/* [npol][nfreq]: 偏波ごとに 1 本の曲線として読める */
	hsize_t dims[2];
	dims[0] = (hsize_t)npol;
	dims[1] = (hsize_t)nfreq;
	hid_t sid = H5Screate_simple(2, dims, NULL);
	hid_t did = H5Dcreate(rcwa_gid, "spectrum", tid, sid,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	int ierr = 0;
	if (did >= 0) {
		if (H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rows) < 0) {
			fprintf(stderr, "*** HDF5 write error (spectrum).\n");
			ierr = 1;
		}
		H5Dclose(did);
	} else {
		ierr = 1;
	}

	H5Sclose(sid);
	H5Tclose(tid);
	H5Gclose(rcwa_gid);
	H5Fclose(file_id);

	return ierr;
}

/* [npol][nfreq][nx][ny] の double データセットを 1 つ書く */
static int write_orders(hid_t gid, const char *name, int npol, int nfreq,
                        int nx, int ny, const double *v)
{
	hsize_t dims[4];
	dims[0] = (hsize_t)npol;
	dims[1] = (hsize_t)nfreq;
	dims[2] = (hsize_t)nx;
	dims[3] = (hsize_t)ny;
	hid_t sid = H5Screate_simple(4, dims, NULL);
	hid_t did = H5Dcreate(gid, name, H5T_NATIVE_DOUBLE, sid,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	int ierr = (did < 0);
	if (did >= 0) {
		if (H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v) < 0) ierr = 1;
		H5Dclose(did);
	}
	H5Sclose(sid);
	return ierr;
}

/* 次数の並び (−N..N) を書く */
static int write_order_index(hid_t gid, const char *name, int n)
{
	int *idx = (int *)malloc((size_t)n * sizeof(int));
	if (idx == NULL) return 1;
	for (int i = 0; i < n; i++) idx[i] = i - (n - 1) / 2;
	hsize_t dims[1];
	dims[0] = (hsize_t)n;
	hid_t sid = H5Screate_simple(1, dims, NULL);
	hid_t did = H5Dcreate(gid, name, H5T_NATIVE_INT, sid,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	int ierr = (did < 0);
	if (did >= 0) {
		if (H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, idx) < 0) ierr = 1;
		H5Dclose(did);
	}
	H5Sclose(sid);
	free(idx);
	return ierr;
}

int rcwa_write_hdf5_orders(const char *filename,
                           int npol, int nfreq, int nx, int ny,
                           const double *R, const double *T)
{
	if ((filename == NULL) || (R == NULL) || (T == NULL) ||
	    (npol < 1) || (nfreq < 1) || (nx < 1) || (ny < 1)) {
		return 1;
	}
	hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		return 1;
	}
	hid_t gid = H5Gopen(file_id, "/rcwa", H5P_DEFAULT);
	if (gid < 0) {
		H5Fclose(file_id);
		return 1;
	}
	int ierr = 0;
	ierr |= write_orders(gid, "R_orders", npol, nfreq, nx, ny, R);
	ierr |= write_orders(gid, "T_orders", npol, nfreq, nx, ny, T);
	ierr |= write_order_index(gid, "order_m", nx);
	ierr |= write_order_index(gid, "order_n", ny);
	H5Gclose(gid);
	H5Fclose(file_id);
	return ierr;
}
