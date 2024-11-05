cmake_minimum_required(VERSION 3.18.3)

project(openRCWA
  DESCRIPTION "A free and open source RCWA solver"
  LANGUAGES C CXX CUDA)

set (CMAKE_CXX_STANDARD 17)
find_package(Threads REQUIRED)

# CUDA を使用するかどうかを設定するオプションを追加
option(USE_CUDA "Enable CUDA support" ON)
if(USE_CUDA)
    #set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++)
    message(STATUS ${CMAKE_CUDA_COMPILER})
    message(STATUS ${CMAKE_CUDA_HOST_COMPILER})
    # CUDA を探す（CMake 3.18 より前）
    #find_package(CUDA REQUIRED)
    # CUDAToolkit を探す
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_ARCHITECTURES 60)
    #message(STATUS "_MPI マクロが定義されました")
    # CUDA のコンパイラフラグを設定
    #set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler -pthread")
else()
    #message(STATUS "CUDA マクロは定義されていません")
endif()


# MPI を使用するかどうかを設定するオプションを追加
option(USE_MPI "Enable MPI support" OFF)
# オプションに基づいて _MPI を定義
if(USE_MPI)
    #before : -D_MPI
    find_package(MPI REQUIRED C)
    #+cuda
    #find_package(MPI REQUIRED)
    add_compile_definitions(_MPI)
    message(STATUS "_MPI マクロが定義されました")
else()
    message(STATUS "_MPI マクロは定義されていません")
endif()

#対応 geometry数
#add_compile_definitions(_ID64)
add_compile_definitions(_ID32)
#add_compile_definitions(_ID16)

#set(CMAKE_SOURCE_DIR ..)
#include_directories(${CMAKE_SOURCE_DIR}/include)
#output folder
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#find_package(glfw3 REQUIRED)
#find_package(OpenGL REQUIRED)

# Find CUDA package
# use CUDA?
# Find OpenMP package
find_package(OpenMP REQUIRED)
# Find MPI package

if(USE_MPI)
    file(GLOB MPI_SOURCES "mpi/*.c")
else()
    file(GLOB SOURCES "sol/*.c")
endif()

# HDF5パッケージの検索

#conda
#set(HDF5_ROOT "$ENV{CONDA_PREFIX}")
find_package(HDF5 REQUIRED COMPONENTS C)
#find_package(HDF5 COMPONENTS)
#find_package(HDF5 COMPONENTS C)
#find_package(hdf5 REQUIRED CONFIG PATHS )
if (HDF5_FOUND)
    message(STATUS "HDF5 が見つかりました。HDF5 サポートを有効にします。")
    # インクルードディレクトリの追加
    #target_include_directories(TimeSeriesHDF5 PRIVATE ${HDF5_INCLUDE_DIRS})

    # ライブラリのリンク
    #target_link_libraries(TimeSeriesHDF5 PRIVATE ${HDF5_LIBRARIES})
    #https://forum.hdfgroup.org/t/cmake-parallel-hdf5-config-scripts-references-hdf-hdf-shared/4605/3
    #set(OP2_HDF5_DEFINITIONS ${HDF5_DEFINITIONS})
    #set(OP2_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    #set(OP2_HDF5_LIBRARIES ${HDF5_LIBRARIES})
    # If HDF5 is build with parallel support it requires MPI, which needs
    # to be enabled
    if (HDF5_IS_PARALLEL)
      if (USE_MPI)
      #  set(OP2_HDF5_INCLUDE_DIRS ${OP2_HDF5_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
      #  set(OP2_HDF5_LIBRARIES ${OP2_HDF5_LIBRARIES} ${MPI_LIBRARIES})
      #  # MPI is always built with HDF5 if available
      #  set(OP2_MPI_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} ${OP2_MPI_INCLUDE_DIRS})
      #  set(OP2_MPI_LIBRARIES ${HDF5_LIBRARIES} ${OP2_MPI_LIBRARIES})
      else()
        message(STATUS "HDF5 is built with parallel support requiring MPI, but MPI is disabled")
        message(STATUS "Disabling HDF5 support")
        #set(OP2_WITH_HDF5 OFF)
      endif()
    else()
      message(STATUS "HDF5 is not built with parallel support.")
    endif()
endif()

# CUDA が見つかった場合のみターゲットをビルド
if (USE_CUDA)
    #if (CUDA_FOUND)
    if (CUDAToolkit_FOUND)
        #CC : nvcc
        message(STATUS "CUDA が見つかりました。CUDA サポートを有効にします。")
        #message(STATUS ${OBJECT_FILES})
        #TODO: Makefile の記述を元に必要なオブジェクトファイルのみに絞る?

        # オブジェクトファイルを生成するカスタムコマンド
        #add_custom_command(
        #    OUTPUT ${OUTPUT_OBJECT}
        #    COMMAND ${CMAKE_C_COMPILER} -c ${SOURCE_FILE} -o ${OUTPUT_OBJECT} -Wall -O2
        #    DEPENDS ${SOURCE_FILE}
        #    COMMENT "Compiling ${SOURCE_FILE} to ${OUTPUT_OBJECT}"
        #    VERBATIM
        #)

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
        if (USE_MPI)
            # Add CUDA source files(reject solve.cu)
            #file(GLOB CUDA_SOURCES "cuda/*.cu")
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
            add_library(orcwa_cuda_mpi_lib ${SOURCES2} ${MPI_SOURCES2} ${CUDA_MPI_SOURCES} ${CUDA_SOURCES2})
            target_include_directories(orcwa_cuda_mpi_lib PUBLIC "${PROJECT_SOURCE_DIR}/include")
            target_link_libraries(orcwa_cuda_mpi_lib PRIVATE ${CUDA_LIBRARIES})
            target_link_libraries(orcwa_cuda_mpi_lib PRIVATE OpenMP::OpenMP_C pthread m)
            target_link_libraries(orcwa_cuda_mpi_lib PRIVATE MPI::MPI_C hdf5::hdf5)
            #error
            #target_link_libraries(orcwa_cuda_mpi_lib PRIVATE MPI::MPI_C hdf5::hdf5-shared)

            #add_executable(orcwa_cuda_mpi ${SOURCES2} ${MPI_SOURCES2} ${CUDA_MPI_SOURCES} ${CUDA_SOURCES2} src/cuda_mpi_Main.cu)
            add_executable(orcwa_cuda_mpi src/cuda_mpi_Main.cu)

            #
            target_include_directories(orcwa_cuda_mpi PUBLIC "${PROJECT_SOURCE_DIR}/include")

            # Link CUDA libraries
            #orcwa_cuda_mpi_lib
            target_link_libraries(orcwa_cuda_mpi PRIVATE orcwa_cuda_mpi_lib)
            target_link_libraries(orcwa_cuda_mpi PRIVATE ${CUDA_LIBRARIES})
            target_link_libraries(orcwa_cuda_mpi PRIVATE OpenMP::OpenMP_C pthread m)
            target_link_libraries(orcwa_cuda_mpi PRIVATE MPI::MPI_C hdf5::hdf5)
        else ()
            # Add CUDA source files
            file(GLOB CUDA_SOURCES "cuda/*.cu")
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
            #reduction_sum.cu

            #add_library(orcwa_cuda_lib ${SOURCES2} ${CUDA_SOURCES})

            #add_executable(orcwa_cuda ${CUDA_SOURCES} src/cuda_Main.cu)
            #add_executable(orcwa_cuda ${CUDA_SOURCES} ${OBJECT_FILES} src/cuda_Main.cu)
            add_executable(orcwa_cuda ${SOURCES2} ${CUDA_SOURCES} src/cuda_Main.cu)
            #add_executable(orcwa_cuda ${SOURCES2} ${CUDA_SOURCES} ${OBJECT_FILES} src/cuda_Main.cu)

            target_include_directories(orcwa_cuda PUBLIC "${PROJECT_SOURCE_DIR}/include")

            # Link CUDA libraries
            target_link_libraries(orcwa_cuda PRIVATE ${CUDA_LIBRARIES})

            # オプションとして C++ 実行ファイルも追加（CUDA 無効時のデフォルト動作用）
            #set(DEFAULT_TARGET orcwa_cuda)

            # Optional: Specify CUDA architectures
            #set(CUDA_ARCHITECTURES "60;70;75;80")
            #add_compile_options(--expt-relaxed-constexpr)

            # Link OpenMP libraries to the target
            if(OpenMP_C_FOUND)
                message(STATUS "OpenMP FOUND.")
                #target_link_libraries(orcwa_cuda PRIVATE OpenMP::OpenMP_C)
                #target_link_libraries(orcwa_cuda PRIVATE OpenMP::OpenMP_C m)
                #target_link_libraries(orcwa_cuda PRIVATE OpenMP::OpenMP_C Threads::Threads)
                #target_link_libraries(orcwa_cuda PRIVATE OpenMP::OpenMP_C Threads::Threads m)
                target_link_libraries(orcwa_cuda PRIVATE OpenMP::OpenMP_C pthread)
                #target_link_libraries(orcwa PRIVATE MPI::MPI_C OpenMP::OpenMP_C m)
                target_link_libraries(orcwa_cuda PRIVATE m hdf5::hdf5)
            else()
                message(STATUS "OpenMP not FOUND.")
                #target_link_libraries(orcwa_cuda PUBLIC m)
                #target_link_libraries(orcwa_cuda PRIVATE Threads::Threads m)
                #target_link_libraries(orcwa_cuda PRIVATE Threads::Threads)
                target_link_libraries(orcwa_cuda PRIVATE pthread)
                target_link_libraries(orcwa_cuda PRIVATE m hdf5::hdf5)
            endif()

        endif()

        # CUDA のプロパティ設定 (オプション)(x)
        #set_target_properties(orcwa_cuda PROPERTIES
        #    CUDA_SEPARABLE_COMPILATION ON          # 分離コンパイルを有効にする(x?)
        #    #CUDA_ARCHITECTURES "60;70"            # 対応する CUDA アーキテクチャを指定
        #    CUDA_ARCHITECTURES "60"                # 対応する CUDA アーキテクチャを指定
        #    CUDA_STANDARD 14                       # CUDA の標準 (C++14) を指定
        #)
    else()
        message(WARNING "CUDA が見つかりませんでした。CUDA サポートは無効です。")
    endif()
else()
    if(USE_MPI)
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
            ${CMAKE_SOURCE_DIR}/sol/Makefile
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

        #file(GLOB MPI_SOURCES "mpi/*.c")
        #target_compile_definitions(orcwa_mpi PRIVATE _MPI)
        add_executable(orcwa_mpi ${SOURCES3} ${MPI_SOURCES} src/mpi_Main.c)

        # Include directories for header files
        target_include_directories(orcwa_mpi PUBLIC "${PROJECT_SOURCE_DIR}/include")

        #target_link_libraries(orcwa_mpi PRIVATE MPI::MPI_C OpenMP::OpenMP_C m)
        target_link_libraries(orcwa_mpi PRIVATE OpenMP::OpenMP_C Threads::Threads m)
        target_link_libraries(orcwa_mpi PRIVATE MPI::MPI_C hdf5::hdf5)
    else ()
        add_executable(orcwa ${SOURCES} src/sol_Main.c)
        #target_link_libraries(orcwa glfw)
        #target_link_libraries(orcwa OpenGL::GL)

        # Include directories for header files
        target_include_directories(orcwa PUBLIC "${PROJECT_SOURCE_DIR}/include")
        #target_link_libraries(orcwa ${CMAKE_DL_LIBS})
        #target_link_libraries(orcwa ${CMAKE_DL_LIBS})
        # Link OpenMP libraries to the target
        if(OpenMP_C_FOUND)
            #target_link_libraries(orcwa PRIVATE OpenMP::OpenMP_C m)
            target_link_libraries(orcwa PRIVATE OpenMP::OpenMP_C m hdf5::hdf5)
            #target_link_libraries(orcwa PRIVATE OpenMP::OpenMP_CXX m)
        else()
            #target_link_libraries(orcwa PRIVATE m)
            target_link_libraries(orcwa PRIVATE m hdf5::hdf5)
        endif()

        # ターゲットのオブジェクトディレクトリを設定
        set_target_properties(orcwa PROPERTIES
            OBJECT_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/obj"
        )
    endif()
endif()

option(USE_POST "Enable Prepost Tool" ON)
if (USE_POST)
    file(GLOB POST_SOURCES "post/*.c")
    set(SOURCES3
        ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
        ${CMAKE_SOURCE_DIR}/sol/farfield.c
        ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
        ${CMAKE_SOURCE_DIR}/sol/utils.c
        ${CMAKE_SOURCE_DIR}/sol/ev3d.c
    )

    add_executable(orcwa_post ${POST_SOURCES} ${SOURCES3} src/post_Main.c)

    # Include directories for header files
    target_include_directories(orcwa_post PUBLIC "${PROJECT_SOURCE_DIR}/include")

    #target_link_libraries(orcwa_post PRIVATE OpenMP::OpenMP_C m)
    target_link_libraries(orcwa_post PRIVATE OpenMP::OpenMP_C m hdf5::hdf5)

endif()

