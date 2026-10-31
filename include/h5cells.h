/*
h5cells.h (OpenRCWA)

/dataNNNNNN の E / H / P / P_loss の書き出し (sol/h5cells.c)。
*/
#ifndef H5CELLS_H
#define H5CELLS_H

#include "hdf5.h"

#ifdef __cplusplus
extern "C" {
#endif

/* group_id のグループに E [1,NFreq2,NN,6] / H [..,6] / P [..,3] /
   P_loss [..,1] を作って書く。P_losses は [NFreq2][NN] の並び。
   buf は NN * 6 個以上の double の作業領域。 */
void write_cell_fields_hdf5(hid_t group_id, const double *P_losses, double *buf, int itime);

#ifdef __cplusplus
}
#endif

#endif  /* H5CELLS_H */
