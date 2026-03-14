cmake_minimum_required(VERSION 3.16)

project(BasefoldOverGR LANGUAGES C CXX ASM)

set(BASEFOLD_BLAKE3_SOURCES
  third_party/blake3/blake3.c
  third_party/blake3/blake3_dispatch.c
  third_party/blake3/blake3_portable.c
)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|X86_64|amd64|AMD64)$" AND UNIX AND
   (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID MATCHES "Clang"))
  list(APPEND BASEFOLD_BLAKE3_SOURCES
       third_party/blake3/blake3_avx2_x86-64_unix.S
       third_party/blake3/blake3_avx512_x86-64_unix.S
       third_party/blake3/blake3_sse2_x86-64_unix.S
       third_party/blake3/blake3_sse41_x86-64_unix.S)
endif()

add_library(pcs_hash_blake3 STATIC ${BASEFOLD_BLAKE3_SOURCES})
target_include_directories(pcs_hash_blake3 PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/blake3
)
target_compile_features(pcs_hash_blake3 PUBLIC c_std_99)
if(NOT (CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|X86_64|amd64|AMD64)$" AND UNIX AND
        (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID MATCHES "Clang")))
  target_compile_definitions(pcs_hash_blake3 PRIVATE
    BLAKE3_USE_NEON=0
    BLAKE3_NO_SSE2
    BLAKE3_NO_SSE41
    BLAKE3_NO_AVX2
    BLAKE3_NO_AVX512
  )
endif()

function(configure_repo_cxx_target target_name)
  target_include_directories(${target_name} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  )
  target_compile_features(${target_name} PUBLIC cxx_std_17)

  if(MSVC)
    target_compile_options(${target_name} PRIVATE /W4)
  else()
    target_compile_options(${target_name} PRIVATE -Wall -Wextra -pedantic)
  endif()
endfunction()

add_library(galoisring_algebra
  src/GaloisRing/Basis.cpp
  src/GaloisRing/utils.cpp
  src/GaloisRing/Inverse.cpp
  src/GaloisRing/HenselLift.cpp
  src/GaloisRing/PrimitiveElement.cpp
  src/GaloisRing/FrobeniusBasis.cpp
)

add_library(pcs_common
  src/PCS/Common/Hash.cpp
  src/PCS/Common/Merkle.cpp
  src/PCS/Common/Transcript.cpp
  src/PCS/Common/Multilinear.cpp
  src/PCS/Common/Sumcheck.cpp
  src/PCS/Common/Profile.cpp
)

add_library(basefold_core
  src/PCS/BaseFold/FoldableCode.cpp
  src/PCS/BaseFold/IOPP.cpp
  src/PCS/BaseFold/BaseFoldPCSCommon.cpp
  src/PCS/BaseFold/BaseFoldPCSCommit.cpp
  src/PCS/BaseFold/BaseFoldPCSProve.cpp
  src/PCS/BaseFold/BaseFoldPCSVerify.cpp
  src/PCS/BaseFold/BaseFoldPCSExtension.cpp
  src/PCS/BaseFold/ProofDeserialize.cpp
  src/PCS/BaseFold/ProofSerialize.cpp
  src/PCS/BaseFold/ProofSize.cpp
)

add_library(z2k_compiler_core
  src/Compiler/Z2k/PCSBackend.cpp
  src/Compiler/Z2k/BaseFoldBackendAdapter.cpp
  src/Compiler/Z2k/FrobeniusPCS.cpp
  src/Compiler/Z2k/FrobeniusProofSerialize.cpp
  src/Compiler/Z2k/RingSwitchPCS.cpp
  src/Compiler/Z2k/RingSwitchProofSerialize.cpp
)

configure_repo_cxx_target(galoisring_algebra)
configure_repo_cxx_target(pcs_common)
configure_repo_cxx_target(basefold_core)
configure_repo_cxx_target(z2k_compiler_core)

find_library(NTL_LIBRARY NAMES ntl REQUIRED)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)
find_library(M_LIBRARY NAMES m)
find_package(OpenMP)

target_link_libraries(galoisring_algebra PUBLIC
  ${NTL_LIBRARY}
  ${GMP_LIBRARY}
)

if(M_LIBRARY)
  target_link_libraries(galoisring_algebra PUBLIC ${M_LIBRARY})
endif()

target_link_libraries(pcs_common PUBLIC
  galoisring_algebra
  pcs_hash_blake3
)

target_link_libraries(basefold_core PUBLIC
  pcs_common
)

target_link_libraries(z2k_compiler_core PUBLIC
  basefold_core
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(pcs_common PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(pcs_common PUBLIC BASEFOLD_USE_OPENMP=1)
  target_link_libraries(basefold_core PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(basefold_core PUBLIC BASEFOLD_USE_OPENMP=1)
endif()

add_executable(test_galois_ring tests/test_galois_ring_basic.cpp)
target_link_libraries(test_galois_ring PRIVATE galoisring_algebra)
target_include_directories(test_galois_ring PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_foldable_codes tests/test_foldable_codes.cpp)
target_link_libraries(test_foldable_codes PRIVATE basefold_core)
target_include_directories(test_foldable_codes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_iopp tests/test_iopp.cpp)
target_link_libraries(test_iopp PRIVATE basefold_core)
target_include_directories(test_iopp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_pcs tests/test_pcs.cpp)
target_link_libraries(test_pcs PRIVATE basefold_core)
target_include_directories(test_pcs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_z2k_ring_switch_pcs tests/test_z2k_ring_switch_pcs.cpp)
target_link_libraries(test_z2k_ring_switch_pcs PRIVATE z2k_compiler_core)
target_include_directories(test_z2k_ring_switch_pcs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_z2k_frobenius_feasibility
  tests/test_z2k_frobenius_feasibility.cpp)
target_link_libraries(test_z2k_frobenius_feasibility PRIVATE pcs_common)
target_include_directories(test_z2k_frobenius_feasibility PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_z2k_frobenius_pcs
  tests/test_z2k_frobenius_pcs.cpp)
target_link_libraries(test_z2k_frobenius_pcs PRIVATE z2k_compiler_core)
target_include_directories(test_z2k_frobenius_pcs PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_z2k_frobenius_bench_cli
  tests/test_z2k_frobenius_bench_cli.cpp)
configure_repo_cxx_target(test_z2k_frobenius_bench_cli)
target_include_directories(test_z2k_frobenius_bench_cli PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_z2k_ring_switch_bench_cli
  tests/test_z2k_ring_switch_bench_cli.cpp)
configure_repo_cxx_target(test_z2k_ring_switch_bench_cli)
target_include_directories(test_z2k_ring_switch_bench_cli PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_path_hygiene tests/test_path_hygiene.cpp)
configure_repo_cxx_target(test_path_hygiene)
target_include_directories(test_path_hygiene PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(bench_pcs_commit bench/bench_pcs_commit.cpp)
target_link_libraries(bench_pcs_commit PRIVATE basefold_core)

add_executable(bench_pcs_eval bench/bench_pcs_eval.cpp)
target_link_libraries(bench_pcs_eval PRIVATE basefold_core)

add_executable(bench_pcs_prove bench/bench_pcs_prove.cpp)
target_link_libraries(bench_pcs_prove PRIVATE basefold_core)

add_executable(bench_pcs_verify bench/bench_pcs_verify.cpp)
target_link_libraries(bench_pcs_verify PRIVATE basefold_core)

add_executable(dump_pcs_eval_artifact bench/dump_pcs_eval_artifact.cpp)
target_link_libraries(dump_pcs_eval_artifact PRIVATE basefold_core)

add_executable(bench_pcs_verify_artifact bench/bench_pcs_verify_artifact.cpp)
target_link_libraries(bench_pcs_verify_artifact PRIVATE basefold_core)

add_executable(bench_z2k_ring_switch_commit bench/bench_z2k_ring_switch_commit.cpp)
target_link_libraries(bench_z2k_ring_switch_commit PRIVATE z2k_compiler_core)

add_executable(bench_z2k_ring_switch_eval bench/bench_z2k_ring_switch_eval.cpp)
target_link_libraries(bench_z2k_ring_switch_eval PRIVATE z2k_compiler_core)

add_executable(bench_z2k_frobenius_commit bench/bench_z2k_frobenius_commit.cpp)
target_link_libraries(bench_z2k_frobenius_commit PRIVATE z2k_compiler_core)

add_executable(bench_z2k_frobenius_eval bench/bench_z2k_frobenius_eval.cpp)
target_link_libraries(bench_z2k_frobenius_eval PRIVATE z2k_compiler_core)

add_executable(bench_z2k_frobenius_outer_commit bench/bench_z2k_frobenius_outer_commit.cpp)
target_link_libraries(bench_z2k_frobenius_outer_commit PRIVATE z2k_compiler_core)

add_executable(bench_z2k_frobenius_outer_prove bench/bench_z2k_frobenius_outer_prove.cpp)
target_link_libraries(bench_z2k_frobenius_outer_prove PRIVATE z2k_compiler_core)

add_executable(bench_z2k_frobenius_outer_verify bench/bench_z2k_frobenius_outer_verify.cpp)
target_link_libraries(bench_z2k_frobenius_outer_verify PRIVATE z2k_compiler_core)

add_dependencies(test_z2k_frobenius_bench_cli
  bench_z2k_frobenius_commit
  bench_z2k_frobenius_eval
  bench_z2k_frobenius_outer_commit
  bench_z2k_frobenius_outer_prove
  bench_z2k_frobenius_outer_verify)

add_dependencies(test_z2k_ring_switch_bench_cli
  bench_z2k_ring_switch_commit
  bench_z2k_ring_switch_eval
  bench_z2k_ring_switch_outer_commit
  bench_z2k_ring_switch_outer_prove
  bench_z2k_ring_switch_outer_verify)

add_executable(bench_z2k_ring_switch_outer_commit bench/bench_z2k_ring_switch_outer_commit.cpp)
target_link_libraries(bench_z2k_ring_switch_outer_commit PRIVATE z2k_compiler_core)

add_executable(bench_z2k_ring_switch_outer_prove bench/bench_z2k_ring_switch_outer_prove.cpp)
target_link_libraries(bench_z2k_ring_switch_outer_prove PRIVATE z2k_compiler_core)

add_executable(bench_z2k_ring_switch_outer_verify bench/bench_z2k_ring_switch_outer_verify.cpp)
target_link_libraries(bench_z2k_ring_switch_outer_verify PRIVATE z2k_compiler_core)

add_executable(calc_iopp_params bench/calc_iopp_params.cpp)
target_compile_features(calc_iopp_params PRIVATE cxx_std_17)

enable_testing()
add_test(NAME test_galois_ring COMMAND $<TARGET_FILE:test_galois_ring>)
add_test(NAME test_foldable_codes COMMAND $<TARGET_FILE:test_foldable_codes>)
add_test(NAME test_iopp COMMAND $<TARGET_FILE:test_iopp>)
add_test(NAME test_pcs COMMAND $<TARGET_FILE:test_pcs>)
add_test(NAME test_z2k_ring_switch_pcs COMMAND $<TARGET_FILE:test_z2k_ring_switch_pcs>)
add_test(NAME test_z2k_frobenius_feasibility
  COMMAND $<TARGET_FILE:test_z2k_frobenius_feasibility>)
add_test(NAME test_z2k_frobenius_pcs
  COMMAND $<TARGET_FILE:test_z2k_frobenius_pcs>)
add_test(NAME test_z2k_frobenius_bench_cli
  COMMAND $<TARGET_FILE:test_z2k_frobenius_bench_cli>)
add_test(NAME test_z2k_ring_switch_bench_cli
  COMMAND $<TARGET_FILE:test_z2k_ring_switch_bench_cli>)
add_test(NAME test_path_hygiene COMMAND $<TARGET_FILE:test_path_hygiene>)
