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

#define COUPLING_DATASET_NAME "coupling"
#define SPARA_DATASET_NAME "s_parameters"
#define CROSS_SECTION_DATASET_NAME "cross_section"
#define ZIN_DATASET_NAME "input_impedance"
typedef struct {
    double magnitude_dB;
    double phase_deg;
} coupling_data_t;

typedef struct {
    double magnitude_dB;
    double phase_deg;
} spara_data_t;

typedef struct {
    double frequency;
    double backward;
    double forward;
} cross_section_data_t;

typedef struct {
    double frequency;
    double rin;
    double xin;
    double gin;
    double bin;
    double ref;
    double vswr;
} input_impedance_data_t;

void readhdf5() {
    // HDF5ファイルのオープン
    hid_t file_id = H5Fopen(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    herr_t status;
    hid_t group_id;

    // /metadata グループのオープン
    hid_t metadata_group_id = H5Gopen(file_id, "/metadata", H5P_DEFAULT);

    // Titleの読み込み
    hid_t dataset_id = H5Dopen(metadata_group_id, "Title", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
    H5Dclose(dataset_id);

    // 各種整数型メタデータの読み込み
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
        dataset_id = H5Dopen(metadata_group_id, metadata[i].name, H5P_DEFAULT);
        status = H5Dread(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
        H5Dclose(dataset_id);
    }

    // Dtの読み込み
    dataset_id = H5Dopen(metadata_group_id, "Dt", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    H5Dclose(dataset_id);

    // Planewaveの読み込み
    dataset_id = H5Dopen(metadata_group_id, "Planewave", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
    H5Dclose(dataset_id);

    // 配列データの読み込み
    struct {
        const char *name;
        double **data;
        size_t size;
    } arrays[] = {
        {"Xn", &Xn, Nx + 1},
        {"Yn", &Yn, Ny + 1},
        {"Zn", &Zn, Nz + 1},
        {"Xc", &Xc, Nx},
        {"Yc", &Yc, Ny},
        {"Zc", &Zc, Nz},
        {"Eiter", &Eiter, Niter},
        {"Hiter", &Hiter, Niter},
        {"VFeed", &VFeed, NFeed * (Solver.maxiter + 1)},
        {"IFeed", &IFeed, NFeed * (Solver.maxiter + 1)},
        {"VPoint", &VPoint, NPoint * (Solver.maxiter + 1)},
        {"Freq1", &Freq1, NFreq1},
        {"Freq2", &Freq2, NFreq2},
        {"Gline", &Gline, NGline * 2 * 3}
    };

    for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
        *arrays[i].data = (double *)malloc(sizeof(double) * arrays[i].size);
        dataset_id = H5Dopen(metadata_group_id, arrays[i].name, H5P_DEFAULT);
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *arrays[i].data);
        H5Dclose(dataset_id);
    }

    //サイズ確保
    Feed   =         (feed_t *)malloc(sizeof(feed_t)      * NFeed);
    Zin    =    (d_complex_t *)malloc(sizeof(d_complex_t) * NFeed * NFreq1);
    Ref    =         (double *)malloc(sizeof(double)      * NFeed * NFreq1);
    Pin[0] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
    Pin[1] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
    Spara  =    (d_complex_t *)malloc(sizeof(d_complex_t) * NPoint * NFreq1);
    Gline  = (double (*)[2][3])malloc(sizeof(double)      * NGline * 2 * 3);

	//plot1d?(Far?)
	size_t size;
    size = NFreq2 * sizeof(d_complex_t *);
    SurfaceEx = (d_complex_t **)malloc(size);
    SurfaceEy = (d_complex_t **)malloc(size);
    SurfaceEz = (d_complex_t **)malloc(size);
    SurfaceHx = (d_complex_t **)malloc(size);
    SurfaceHy = (d_complex_t **)malloc(size);
    SurfaceHz = (d_complex_t **)malloc(size);

	//plot2d?
    size = NN * NFreq2 * sizeof(float);
    cEx_r = (float *)malloc(size);
    cEx_i = (float *)malloc(size);
    cEy_r = (float *)malloc(size);
    cEy_i = (float *)malloc(size);
    cEz_r = (float *)malloc(size);
    cEz_i = (float *)malloc(size);
    cHx_r = (float *)malloc(size);
    cHx_i = (float *)malloc(size);
    cHy_r = (float *)malloc(size);
    cHy_i = (float *)malloc(size);
    cHz_r = (float *)malloc(size);
    cHz_i = (float *)malloc(size);

    // Surfaceデータの読み込み
    dataset_id = H5Dopen(metadata_group_id, "NSurface", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NSurface);
    H5Dclose(dataset_id);

    // surface_t構造体に対応する複合データ型を再定義
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
    H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);

	//surface_t
    //Surface = (double *)malloc(sizeof(double) * NSurface);
    if (Surface == NULL) {
        Surface = (surface_t *)malloc(NSurface * sizeof(surface_t));
    }
    dataset_id = H5Dopen(metadata_group_id, "Surface", H5P_DEFAULT);
    status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: Surface\n");
    }
    H5Dclose(dataset_id);
    H5Tclose(memtype);

    //size = NSurface * sizeof(surface_t);
    //Surface = (surface_t *)malloc(size);

//
    // coupling データの読み込み
    coupling_data_t coupling_data[NFreq1][NFeed][NPoint];

    // coupling_data_t のデータ型を定義
    hid_t datatype_id;
    datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(coupling_data_t));
    H5Tinsert(datatype_id, "magnitude_dB", HOFFSET(coupling_data_t, magnitude_dB), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "phase_deg", HOFFSET(coupling_data_t, phase_deg), H5T_NATIVE_DOUBLE);

    dataset_id = H5Dopen(metadata_group_id, COUPLING_DATASET_NAME, H5P_DEFAULT);
    status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, coupling_data);
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: %s\n", COUPLING_DATASET_NAME);
    }
    //coupling_data から
    H5Dclose(dataset_id);
    H5Tclose(datatype_id);

    // s-parameters データの読み込み
    spara_data_t spara_data[NFreq1][NPoint];

    // データタイプを作成
    datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(spara_data_t));
    H5Tinsert(datatype_id, "magnitude_dB", HOFFSET(spara_data_t, magnitude_dB), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "phase_deg", HOFFSET(spara_data_t, phase_deg), H5T_NATIVE_DOUBLE);

    dataset_id = H5Dopen(metadata_group_id, SPARA_DATASET_NAME, H5P_DEFAULT);
    status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, spara_data);
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: %s\n", SPARA_DATASET_NAME);
    }
    // spara_data から Spara への加工]
    if (Spara == NULL) {
        Spara = (d_complex_t *)malloc(sizeof(d_complex_t) * NPoint * NFreq1);
    }
    for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
        for (int ipoint = 0; ipoint < NPoint; ipoint++) {
            const int id = (ipoint * NFreq1) + ifreq;
            //Spara[id] = d_polar(pow(10, spara_data[ifreq][ipoint].magnitude_dB / 20), spara_data[ifreq][ipoint].phase_deg * (M_PI / 180.0));
            const double mag = pow(10, spara_data[ifreq][ipoint].magnitude_dB / 20);
            const double ph  = spara_data[ifreq][ipoint].phase_deg * (M_PI / 180.0);
            Spara[id] = d_complex(mag * cos(ph), mag * sin(ph));
        }
    }
    H5Dclose(dataset_id);
    H5Tclose(datatype_id);

    // cross_section データの読み込み
/*
    cross_section_data_t cross_section_data[NFreq2];
    
    // cross_section_data_t 型のデータ型を定義
    datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(cross_section_data_t));
    H5Tinsert(datatype_id, "frequency", HOFFSET(cross_section_data_t, frequency), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "backward", HOFFSET(cross_section_data_t, backward), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "forward", HOFFSET(cross_section_data_t, forward), H5T_NATIVE_DOUBLE);

    dataset_id = H5Dopen(metadata_group_id, CROSS_SECTION_DATASET_NAME, H5P_DEFAULT);
    status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, cross_section_data);
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: %s\n", CROSS_SECTION_DATASET_NAME);
    }
    cross_section_data から 
    H5Dclose(dataset_id);
    H5Tclose(datatype_id);
*/

    // input_impedance データの読み込み
    fprintf(stdout, "input_impedance_data (start)\n");
    input_impedance_data_t input_impedance_data[NFeed][NFreq1];
    
    hid_t memtype_id;

    // メモリタイプを定義
    memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(input_impedance_data_t));
    H5Tinsert(memtype_id, "frequency", HOFFSET(input_impedance_data_t, frequency), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "rin", HOFFSET(input_impedance_data_t, rin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "xin", HOFFSET(input_impedance_data_t, xin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "gin", HOFFSET(input_impedance_data_t, gin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "bin", HOFFSET(input_impedance_data_t, bin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "ref", HOFFSET(input_impedance_data_t, ref), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "vswr", HOFFSET(input_impedance_data_t, vswr), H5T_NATIVE_DOUBLE);

    dataset_id = H5Dopen(metadata_group_id, ZIN_DATASET_NAME, H5P_DEFAULT);
    status = H5Dread(dataset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, input_impedance_data);
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: %s\n", ZIN_DATASET_NAME);
    }
    //input_impedance_data から Zin への加工
    //Zin = (d_complex_t *)malloc(sizeof(d_complex_t) * NFeed * NFreq1);
    //Ref = (double *)malloc(sizeof(double)      * NFeed * NFreq1);
    for (int ifeed = 0; ifeed < NFeed; ifeed++) {
        for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
            const int id = (ifeed * NFreq1) + ifreq;
            Zin[id].r = input_impedance_data[ifeed][ifreq].rin;
            Zin[id].i = input_impedance_data[ifeed][ifreq].xin;
            Ref[id] = input_impedance_data[ifeed][ifreq].ref;
        }
    }
    fprintf(stdout, "input_impedance_data (end)\n");
    H5Dclose(dataset_id);
    H5Tclose(memtype_id);

    // メタデータグループのクローズ
    H5Gclose(metadata_group_id);

    // HDF5ファイルのクローズ
    H5Fclose(file_id);
}

