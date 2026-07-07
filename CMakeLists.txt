cmake_minimum_required(VERSION 3.18.3)

project(openRCWA
  DESCRIPTION "A free and open source RCWA solver"
  LANGUAGES C CXX)

# cmake/Modules ディレクトリをモジュールの検索パスに追加
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(sploosh_utils)
# compile option
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# CUDA を使用するかどうかを設定するオプションを追加
option(WITH_CUDA "Enable CUDA support" OFF)
if(WITH_CUDA)
    enable_language(CUDA)
    #set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++)
    message(STATUS ${CMAKE_CUDA_COMPILER})
    message(STATUS ${CMAKE_CUDA_HOST_COMPILER})
    # CUDAToolkit を探す
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_ARCHITECTURES 60)

    # MAGMA (差分法 GPU ソルバの固有値計算で使用)
    find_path(MAGMA_INCLUDE_DIR magma_v2.h
        PATHS /usr/local/include /usr/include)
    find_library(MAGMA_LIBRARY NAMES magma
        PATHS /usr/local/lib /usr/lib /usr/lib/x86_64-linux-gnu)
    if(MAGMA_INCLUDE_DIR AND MAGMA_LIBRARY)
        message(STATUS "MAGMA が見つかりました: ${MAGMA_LIBRARY}")
        set(MAGMA_LIBRARIES ${MAGMA_LIBRARY} OpenMP::OpenMP_CXX)
    else()
        message(WARNING "MAGMA が見つかりません。GPU の差分法ソルバはビルドされません。")
        set(MAGMA_INCLUDE_DIR "")
        set(MAGMA_LIBRARIES "")
    endif()
endif()

# MPI を使用するかどうかを設定するオプションを追加
option(WITH_MPI "Enable MPI support" OFF)
# オプションに基づいて _MPI を定義
if(WITH_MPI)
    find_package(MPI REQUIRED C)
    add_compile_definitions(_MPI)
    message(STATUS "_MPI マクロが定義されました")
else()
    message(STATUS "_MPI マクロは定義されていません")
endif()

# MKL はオプション扱い: 見つかれば LAPACKE 経由の固有値ソルバを使用し、
# 見つからなければ Eigen::ComplexEigenSolver にフォールバックする
find_package(MKL)
if(MKL_FOUND)
    message(STATUS "MKL libraries: ${MKL_LIBRARIES}")
    message(STATUS "MKL include directory: ${MKL_INCLUDE_DIR}")
    message(STATUS "MKL OpenMP type: ${MKL_OPENMP_TYPE}")
    message(STATUS "MKL OpenMP library: ${MKL_OPENMP_LIBRARY}")
    add_definitions(-DEIGEN_WITH_MKL_ALL -DUSE_MKL)
else()
    message(STATUS "MKL が見つかりません。Eigen の固有値ソルバを使用します。")
    # gdstk/utils.cpp が LAPACK (dgesv_) を直接呼ぶため、システム LAPACK で代替する
    find_package(LAPACK REQUIRED)
    set(MKL_LIBRARIES ${LAPACK_LIBRARIES})
    set(MKL_INCLUDE_DIR "")

    # LAPACKE (zgeev/cgeev) があれば、RCWA の固有値計算で
    # Eigen::ComplexEigenSolver より高速な LAPACKE 経由のソルバを使用する
    find_path(LAPACKE_INCLUDE_DIR lapacke.h)
    find_library(LAPACKE_LIBRARY NAMES lapacke)
    if(LAPACKE_INCLUDE_DIR AND LAPACKE_LIBRARY)
        message(STATUS "LAPACKE が見つかりました: ${LAPACKE_LIBRARY}")
        add_definitions(-DUSE_LAPACKE)
        list(APPEND MKL_LIBRARIES ${LAPACKE_LIBRARY})
        set(MKL_INCLUDE_DIR ${LAPACKE_INCLUDE_DIR})
    else()
        message(STATUS "LAPACKE が見つかりません。Eigen の固有値ソルバを使用します。")
    endif()
endif()

find_package(Eigen3 CONFIG REQUIRED)  # バージョン指定なし: Homebrew の Eigen 5.x も受理する
if (EIGEN3_FOUND)
   message(STATUS "Found Eigen version ${EIGEN3_VERSION_STRING}: ${EIGEN3_INCLUDE_DIRS}")
endif ()

find_package(TBB)
if (TBB_FOUND)
    message(STATUS "Found Intel TBB version ${TBB_VERSION}: ${TBB_INCLUDE_DIRS} ${TBB_LIBRARIES}")
else()
    message(STATUS "TBB が見つかりません。TBB なしでビルドします。")
    set(TBB_LIBRARIES "")
    set(TBB_INCLUDE_DIRS "")
endif ()

include_directories(
    ${MKL_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
)
if (WITH_CUDA)
    include_directories(${CUDAToolkit_INCLUDE_DIRS})
endif()

#対応 geometry数
#add_compile_definitions(_ID64)
add_compile_definitions(_ID32)
#add_compile_definitions(_ID16)

#output folder
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

if(WITH_MPI)
    file(GLOB MPI_SOURCES "mpi/*.c")
else()
    file(GLOB SOURCES "sol/*.c")
    # setupNear.c は旧バージョン由来で未使用 (orcwa.h にない Fnorm を参照する)
    list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/sol/setupNear.c")
endif()

# HDF5パッケージの検索
find_package(HDF5 REQUIRED COMPONENTS C)
if (HDF5_FOUND)
    message(STATUS "HDF5 が見つかりました。HDF5 サポートを有効にします。")
    if (HDF5_IS_PARALLEL)
      if (NOT WITH_MPI)
        message(STATUS "HDF5 is built with parallel support requiring MPI, but MPI is disabled")
        message(STATUS "Disabling HDF5 support")
      endif()
    else()
      message(STATUS "HDF5 is not built with parallel support.")
    endif()
endif()

# ============================================================
# RCWA ソルバ (rcwa/ core/ gdstk/ gpu/)
# OpenFDTD 由来の FDTD ソルバとは独立にビルドする
# ============================================================
option(WITH_RCWA "Build RCWA solver and tests" ON)

set(ROOT_SRC_DIR ${CMAKE_SOURCE_DIR})
set(ROOT_TESTS_DIR ${CMAKE_SOURCE_DIR}/bin/tests)

file(GLOB RCWA_SRC_FILES "${CMAKE_SOURCE_DIR}/rcwa/*.cpp")
file(GLOB GDSTK_SRC_FILES
    "${CMAKE_SOURCE_DIR}/gdstk/*.cpp"
    "${CMAKE_SOURCE_DIR}/gdstk/clipperlib/*.cpp")
file(GLOB CORE_SRC_FILES "${CMAKE_SOURCE_DIR}/core/*.cpp")
if(WITH_CUDA)
    file(GLOB GPU_SRC_FILES
        "${CMAKE_SOURCE_DIR}/gpu/*.cu"
        "${CMAKE_SOURCE_DIR}/gpu/*.cpp")
endif()

# rcwa/, core/, gdstk/, utils/ を "rcwa/xxx.h" 形式で include できるようにする
include_directories(${CMAKE_SOURCE_DIR})

if(WITH_RCWA)
    # tests/CMakeLists.txt は BUILD_GPU を参照する
    set(BUILD_GPU ${WITH_CUDA})
    add_subdirectory(tests)
endif()

# ============================================================
# FDTD ソルバ (OpenFDTD 由来: sol/ mpi/ cuda/ cuda_mpi/)
# ============================================================
# CUDA が見つかった場合のみターゲットをビルド
if (WITH_CUDA)
    if (CUDAToolkit_FOUND)
        #CC : nvcc
        message(STATUS "CUDA が見つかりました。CUDA サポートを有効にします。")

        set(SOURCES2
            ${CMAKE_SOURCE_DIR}/sol/input_data.c
            ${CMAKE_SOURCE_DIR}/sol/input2.c
            ${CMAKE_SOURCE_DIR}/sol/monitor.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree1.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree3.c
            ${CMAKE_SOURCE_DIR}/sol/setup.c
            ${CMAKE_SOURCE_DIR}/sol/setupSize.c
            ${CMAKE_SOURCE_DIR}/sol/setupDispersion.c
            ${CMAKE_SOURCE_DIR}/sol/setupDft.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPml.c
            ${CMAKE_SOURCE_DIR}/sol/vfeed.c
            ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
            ${CMAKE_SOURCE_DIR}/sol/outputChars.c
            ${CMAKE_SOURCE_DIR}/sol/outputZin.c
            ${CMAKE_SOURCE_DIR}/sol/outputSpara.c
            ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
            ${CMAKE_SOURCE_DIR}/sol/outputCross.c
            ${CMAKE_SOURCE_DIR}/sol/setup_vector.c
            ${CMAKE_SOURCE_DIR}/sol/setupId.c
            ${CMAKE_SOURCE_DIR}/sol/ingeometry.c
            ${CMAKE_SOURCE_DIR}/sol/geomlines.c
            ${CMAKE_SOURCE_DIR}/sol/farfield.c
            ${CMAKE_SOURCE_DIR}/sol/writeout.c
            ${CMAKE_SOURCE_DIR}/sol/utils.c
            ${CMAKE_SOURCE_DIR}/sol/cputime.c
            ${CMAKE_SOURCE_DIR}/sol/plot3dGeom.c
            ${CMAKE_SOURCE_DIR}/sol/ev3d.c
        )

        # CUDA の設定を指定
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_60,code=sm_60")

        # Add CUDA executable
        if (WITH_MPI)
            # Add CUDA source files(reject solve.cu)
            set(CUDA_SOURCES2
                ${CMAKE_SOURCE_DIR}/cuda/average.cu
                ${CMAKE_SOURCE_DIR}/cuda/check_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/cuda_memory.cu
                ${CMAKE_SOURCE_DIR}/cuda/dftNear3d.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/efeed.cu
                ${CMAKE_SOURCE_DIR}/cuda/eload.cu
                ${CMAKE_SOURCE_DIR}/cuda/info_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/initfield.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree2_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree3_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/murH.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/setup_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/vpoint.cu
            )

            set(MPI_SOURCES2
                ${CMAKE_SOURCE_DIR}/mpi/comm.c
                ${CMAKE_SOURCE_DIR}/mpi/setup_mpi.c
            )

            file(GLOB CUDA_MPI_SOURCES "cuda_mpi/*.cu")
            add_library(orcwa_cuda_libmpi ${SOURCES2} ${MPI_SOURCES2} ${CUDA_MPI_SOURCES} ${CUDA_SOURCES2})
            target_include_directories(orcwa_cuda_libmpi PUBLIC "${PROJECT_SOURCE_DIR}/include")
            target_link_libraries(orcwa_cuda_libmpi PRIVATE ${CUDA_LIBRARIES})
            target_link_libraries(orcwa_cuda_libmpi PRIVATE OpenMP::OpenMP_C pthread m)
            target_link_libraries(orcwa_cuda_libmpi PRIVATE MPI::MPI_C hdf5::hdf5)

            add_executable(orcwa_cuda_mpi src/cuda_mpi_Main.cu)
            target_include_directories(orcwa_cuda_mpi PUBLIC "${PROJECT_SOURCE_DIR}/include")
            target_link_libraries(orcwa_cuda_mpi PRIVATE orcwa_cuda_libmpi)
            target_link_libraries(orcwa_cuda_mpi PRIVATE ${CUDA_LIBRARIES})
            target_link_libraries(orcwa_cuda_mpi PRIVATE OpenMP::OpenMP_C pthread m)
            target_link_libraries(orcwa_cuda_mpi PRIVATE MPI::MPI_C hdf5::hdf5)
        else ()
            # Add CUDA source files
            set(CUDA_SOURCES
                ${CMAKE_SOURCE_DIR}/cuda/average.cu
                ${CMAKE_SOURCE_DIR}/cuda/check_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/cuda_memory.cu
                ${CMAKE_SOURCE_DIR}/cuda/dftNear3d.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/efeed.cu
                ${CMAKE_SOURCE_DIR}/cuda/eload.cu
                ${CMAKE_SOURCE_DIR}/cuda/info_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/initfield.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree2_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree3_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/murH.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/setup_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/solve.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/vpoint.cu
            )

            add_executable(orcwa_cuda ${SOURCES2} ${CUDA_SOURCES} src/cuda_Main.cu)
            target_include_directories(orcwa_cuda PUBLIC "${PROJECT_SOURCE_DIR}/include")

            # Link CUDA libraries
            target_link_libraries(orcwa_cuda PRIVATE ${CUDA_LIBRARIES})

            # Link OpenMP libraries to the target
            if(OpenMP_C_FOUND)
                message(STATUS "OpenMP FOUND.")
                target_link_libraries(orcwa_cuda PRIVATE OpenMP::OpenMP_C pthread)
                target_link_libraries(orcwa_cuda PRIVATE m hdf5::hdf5)
            else()
                message(STATUS "OpenMP not FOUND.")
                target_link_libraries(orcwa_cuda PRIVATE pthread)
                target_link_libraries(orcwa_cuda PRIVATE m hdf5::hdf5)
            endif()
        endif()
    else()
        message(WARNING "CUDA が見つかりませんでした。CUDA サポートは無効です。")
    endif()
else()
    if(WITH_MPI)
        #solve.c reject
        set(SOURCES3
            ${CMAKE_SOURCE_DIR}/sol/average.c
            ${CMAKE_SOURCE_DIR}/sol/cputime.c
            ${CMAKE_SOURCE_DIR}/sol/dftNear3d.c
            ${CMAKE_SOURCE_DIR}/sol/dispersionEx.c
            ${CMAKE_SOURCE_DIR}/sol/dispersionEy.c
            ${CMAKE_SOURCE_DIR}/sol/dispersionEz.c
            ${CMAKE_SOURCE_DIR}/sol/efeed.c
            ${CMAKE_SOURCE_DIR}/sol/eload.c
            ${CMAKE_SOURCE_DIR}/sol/ev3d.c
            ${CMAKE_SOURCE_DIR}/sol/farfield.c
            ${CMAKE_SOURCE_DIR}/sol/geomlines.c
            ${CMAKE_SOURCE_DIR}/sol/ingeometry.c
            ${CMAKE_SOURCE_DIR}/sol/initfield.c
            ${CMAKE_SOURCE_DIR}/sol/input_data.c
            ${CMAKE_SOURCE_DIR}/sol/input2.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree1.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree2.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree3.c
            ${CMAKE_SOURCE_DIR}/sol/monitor.c
            ${CMAKE_SOURCE_DIR}/sol/murH.c
            ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
            ${CMAKE_SOURCE_DIR}/sol/outputChars.c
            ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
            ${CMAKE_SOURCE_DIR}/sol/outputCross.c
            ${CMAKE_SOURCE_DIR}/sol/outputSpara.c
            ${CMAKE_SOURCE_DIR}/sol/outputZin.c
            ${CMAKE_SOURCE_DIR}/sol/pbcx.c
            ${CMAKE_SOURCE_DIR}/sol/pbcy.c
            ${CMAKE_SOURCE_DIR}/sol/pbcz.c
            ${CMAKE_SOURCE_DIR}/sol/plot3dGeom.c
            ${CMAKE_SOURCE_DIR}/sol/pmlEx.c
            ${CMAKE_SOURCE_DIR}/sol/pmlEy.c
            ${CMAKE_SOURCE_DIR}/sol/pmlEz.c
            ${CMAKE_SOURCE_DIR}/sol/pmlHx.c
            ${CMAKE_SOURCE_DIR}/sol/pmlHy.c
            ${CMAKE_SOURCE_DIR}/sol/pmlHz.c
            ${CMAKE_SOURCE_DIR}/sol/setup.c
            ${CMAKE_SOURCE_DIR}/sol/setup_vector.c
            ${CMAKE_SOURCE_DIR}/sol/setupDft.c
            ${CMAKE_SOURCE_DIR}/sol/setupDispersion.c
            ${CMAKE_SOURCE_DIR}/sol/setupId.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPml.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupSize.c
            ${CMAKE_SOURCE_DIR}/sol/updateEx.c
            ${CMAKE_SOURCE_DIR}/sol/updateEy.c
            ${CMAKE_SOURCE_DIR}/sol/updateEz.c
            ${CMAKE_SOURCE_DIR}/sol/updateHx.c
            ${CMAKE_SOURCE_DIR}/sol/updateHy.c
            ${CMAKE_SOURCE_DIR}/sol/updateHz.c
            ${CMAKE_SOURCE_DIR}/sol/utils.c
            ${CMAKE_SOURCE_DIR}/sol/vfeed.c
            ${CMAKE_SOURCE_DIR}/sol/vpoint.c
            ${CMAKE_SOURCE_DIR}/sol/writeout.c
        )

        add_executable(orcwa_mpi ${SOURCES3} ${MPI_SOURCES} src/mpi_Main.c)

        # Include directories for header files
        target_include_directories(orcwa_mpi PUBLIC "${PROJECT_SOURCE_DIR}/include")

        target_link_libraries(orcwa_mpi PRIVATE OpenMP::OpenMP_C Threads::Threads m)
        target_link_libraries(orcwa_mpi PRIVATE MPI::MPI_C hdf5::hdf5)
    else ()
        add_executable(orcwa ${SOURCES} src/sol_Main.c)

        # Include directories for header files
        target_include_directories(orcwa PUBLIC "${PROJECT_SOURCE_DIR}/include")

        # Link OpenMP libraries to the target
        if(OpenMP_C_FOUND)
            target_link_libraries(orcwa PRIVATE OpenMP::OpenMP_C m hdf5::hdf5)
        else()
            target_link_libraries(orcwa PRIVATE m hdf5::hdf5)
        endif()
    endif()
endif()

option(WITH_POST "Enable Prepost Tool" ON)
if (WITH_POST)
    file(GLOB POST_SOURCES "post/*.c")
    set(POST_SOL_SOURCES
        ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
        ${CMAKE_SOURCE_DIR}/sol/farfield.c
        ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
        ${CMAKE_SOURCE_DIR}/sol/utils.c
        ${CMAKE_SOURCE_DIR}/sol/ev3d.c
    )

    add_executable(orcwa_post ${POST_SOURCES} ${POST_SOL_SOURCES} src/post_Main.c)

    # Include directories for header files
    target_include_directories(orcwa_post PUBLIC "${PROJECT_SOURCE_DIR}/include")

    target_link_libraries(orcwa_post PRIVATE OpenMP::OpenMP_C m hdf5::hdf5)
endif()

option(WITH_DATALIB "Build C data-creation library and samples" ON)
if(WITH_DATALIB)
    add_library(orcwa_datalib STATIC datalib/orcwa_datalib.c)
    target_include_directories(orcwa_datalib PUBLIC "${CMAKE_SOURCE_DIR}/datalib")

    add_executable(sample_rcwa_grating datalib/sample_rcwa_grating.c)
    target_link_libraries(sample_rcwa_grating PRIVATE orcwa_datalib)
endif()
