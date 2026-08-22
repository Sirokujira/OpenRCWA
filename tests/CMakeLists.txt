if (BUILD_GPU)
  include_directories(
    ${MKL_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${MAGMA_INCLUDE_DIR}
  )
else()
  include_directories(
    ${MKL_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
  )
endif()


# GPU 版 RCWA のテスト (test_gpu / test_metasurface_gpu / test_oe_model_gpu /
# test_gds_gpu) は gdstk/ と core/ に依存する (下の旧世代テストと同じ理由で
# Windows ではリンクできない) うえ、cuSOLVER / cuBLAS / MAGMA を要求する。
# orcwa_cuda (FDTD の CUDA 版) には不要なので、旧世代テストと同じスイッチで
# 切り離す (MSVC では既定 OFF)。
if (BUILD_GPU AND WITH_RCWA_LEGACY_TESTS)
  add_exe(test_gpu
      test_gpu.cpp
      ${RCWA_SRC_FILES}
      ${GDSTK_SRC_FILES}
      ${CORE_SRC_FILES}
      ${GPU_SRC_FILES}
    LINK_LIBRARIES
      ${MKL_LIBRARIES}
      ${TBB_LIBRARIES}
      ${MAGMA_LIBRARIES}
      ${MATH_LIB}
      cusolver
      cublas
    OUT_DIR
      ${ROOT_TESTS_DIR})
  set_property(TARGET test_gpu PROPERTY CXX_STANDARD 17)

  add_exe(test_metasurface_gpu
      test_metasurface_gpu.cpp
      ${RCWA_SRC_FILES}
      ${GDSTK_SRC_FILES}
      ${CORE_SRC_FILES}
      ${GPU_SRC_FILES}
    LINK_LIBRARIES
      ${MKL_LIBRARIES}
      ${TBB_LIBRARIES}
      ${MAGMA_LIBRARIES}
      ${MATH_LIB}
      cusolver
      cublas
    OUT_DIR
      ${ROOT_TESTS_DIR})
  set_property(TARGET test_metasurface_gpu PROPERTY CXX_STANDARD 17)

  add_exe(test_oe_model_gpu
      test_oe_model_gpu.cpp
      ${RCWA_SRC_FILES}
      ${GDSTK_SRC_FILES}
      ${CORE_SRC_FILES}
      ${GPU_SRC_FILES}
    LINK_LIBRARIES
      ${MKL_LIBRARIES}
      ${TBB_LIBRARIES}
      ${MAGMA_LIBRARIES}
      ${MATH_LIB}
      cusolver
      cublas
    OUT_DIR
      ${ROOT_TESTS_DIR})
  set_property(TARGET test_oe_model_gpu PROPERTY CXX_STANDARD 17)

  add_exe(test_gds_gpu
      test_gds_gpu.cpp
      ${RCWA_SRC_FILES}
      ${GDSTK_SRC_FILES}
      ${CORE_SRC_FILES}
      ${GPU_SRC_FILES}
    LINK_LIBRARIES
      ${MKL_LIBRARIES}
      ${TBB_LIBRARIES}
      ${MAGMA_LIBRARIES}
      ${MATH_LIB}
      cusolver
      cublas
    OUT_DIR
      ${ROOT_TESTS_DIR})
  set_property(TARGET test_gds_gpu PROPERTY CXX_STANDARD 17)
endif()

# ── 旧世代テスト (gdstk/ + core/ 依存) ───────────────────────────────
# gdstk/utils.cpp が LAPACK の dgesv_ を直接呼ぶため、標準 LAPACK の無い
# Windows ではリンクできない。WITH_RCWA_LEGACY_TESTS=OFF で丸ごと除外する。
# RCWA 本体と test_rcwa_input はこれらに依存しないので影響を受けない。
if(WITH_RCWA_LEGACY_TESTS)

add_exe(test_gds_cpu
    test_gds_cpu.cpp
    ${RCWA_SRC_FILES}
    ${GDSTK_SRC_FILES}
    ${ROOT_SRC_DIR}/core/LayerSampler.cpp
    ${ROOT_SRC_DIR}/core/VarLayerSampler.cpp
    ${ROOT_SRC_DIR}/core/RedhefferIntegrator.cpp
    ${ROOT_SRC_DIR}/core/RCWAIntegrator.cpp
    ${ROOT_SRC_DIR}/core/DifferentialIntegrator.cpp
    ${ROOT_SRC_DIR}/core/GDSIISampler.cpp
    ${ROOT_SRC_DIR}/core/Material.cpp
  LINK_LIBRARIES
    ${MKL_LIBRARIES}
    ${TBB_LIBRARIES}
    ${MATH_LIB}
  OUT_DIR
    ${ROOT_TESTS_DIR})
set_property(TARGET test_gds_cpu PROPERTY CXX_STANDARD 17)

add_exe(test_cpu
    test_cpu.cpp
    ${RCWA_SRC_FILES}
    ${ROOT_SRC_DIR}/core/LayerSampler.cpp
    ${ROOT_SRC_DIR}/core/VarLayerSampler.cpp
    ${ROOT_SRC_DIR}/core/RedhefferIntegrator.cpp
    ${ROOT_SRC_DIR}/core/RCWAIntegrator.cpp
    ${ROOT_SRC_DIR}/core/DifferentialIntegrator.cpp
  LINK_LIBRARIES
    ${MKL_LIBRARIES}
    ${TBB_LIBRARIES}
    ${MATH_LIB}
  OUT_DIR
    ${ROOT_TESTS_DIR})
set_property(TARGET test_cpu PROPERTY CXX_STANDARD 17)

add_exe(test_metasurface_cpu
    test_metasurface_cpu.cpp
    ${RCWA_SRC_FILES}
    ${ROOT_SRC_DIR}/core/LayerSampler.cpp
    ${ROOT_SRC_DIR}/core/VarLayerSampler.cpp
    ${ROOT_SRC_DIR}/core/RedhefferIntegrator.cpp
    ${ROOT_SRC_DIR}/core/RCWAIntegrator.cpp
    ${ROOT_SRC_DIR}/core/DifferentialIntegrator.cpp
  LINK_LIBRARIES
    ${MKL_LIBRARIES}
    ${TBB_LIBRARIES}
    ${MATH_LIB}
  OUT_DIR
    ${ROOT_TESTS_DIR})
set_property(TARGET test_metasurface_cpu PROPERTY CXX_STANDARD 17)

add_exe(test_cpu_diff
    test_cpu_differential_method.cpp
    ${RCWA_SRC_FILES}
    ${ROOT_SRC_DIR}/core/LayerSampler.cpp
    ${ROOT_SRC_DIR}/core/VarLayerSampler.cpp
    ${ROOT_SRC_DIR}/core/RedhefferIntegrator.cpp
    ${ROOT_SRC_DIR}/core/RCWAIntegrator.cpp
    ${ROOT_SRC_DIR}/core/DifferentialIntegrator.cpp
  LINK_LIBRARIES
    ${MKL_LIBRARIES}
    ${TBB_LIBRARIES}
    ${MATH_LIB}
  OUT_DIR
    ${ROOT_TESTS_DIR})
set_property(TARGET test_cpu_diff PROPERTY CXX_STANDARD 17)

endif()  # WITH_RCWA_LEGACY_TESTS

# RCWAInput パーサ + RCWADriver の単体テスト。
# rcwa/*.cpp のみに依存するため全プラットフォームでビルドできる。
add_exe(test_rcwa_input
    test_rcwa_input.cpp
    ${RCWA_SRC_FILES}
  LINK_LIBRARIES
    ${MKL_LIBRARIES}
    ${TBB_LIBRARIES}
    ${MATH_LIB}
  OUT_DIR
    ${ROOT_TESTS_DIR})
set_property(TARGET test_rcwa_input PROPERTY CXX_STANDARD 17)