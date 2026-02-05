cmake_minimum_required(VERSION 3.16)

project(BasefoldOverGR LANGUAGES CXX)

add_library(galoisring
  src/GaloisRing/utils.cpp
  src/GaloisRing/Inverse.cpp
  src/GaloisRing/HenselLift.cpp
  src/GaloisRing/PrimitiveElement.cpp
  src/BaseFold/FoldableCode.cpp
  src/BaseFold/IOPP.cpp
  src/BaseFold/Multilinear.cpp
  src/BaseFold/Sumcheck.cpp
  src/BaseFold/BaseFoldPCS.cpp
  src/BaseFold/ProofSize.cpp
)

target_include_directories(galoisring PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(galoisring PUBLIC cxx_std_17)

if(MSVC)
  target_compile_options(galoisring PRIVATE /W4)
else()
  target_compile_options(galoisring PRIVATE -Wall -Wextra -pedantic)
endif()

find_library(NTL_LIBRARY NAMES ntl REQUIRED)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)
find_library(M_LIBRARY NAMES m)
find_package(OpenSSL REQUIRED)

target_link_libraries(galoisring PUBLIC
  ${NTL_LIBRARY}
  ${GMP_LIBRARY}
  OpenSSL::Crypto
)
if(M_LIBRARY)
  target_link_libraries(galoisring PUBLIC ${M_LIBRARY})
endif()

add_executable(test_galois_ring tests/test_galois_ring_basic.cpp)
target_link_libraries(test_galois_ring PRIVATE galoisring)
target_include_directories(test_galois_ring PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_foldable_codes tests/test_foldable_codes.cpp)
target_link_libraries(test_foldable_codes PRIVATE galoisring)
target_include_directories(test_foldable_codes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_iopp tests/test_iopp.cpp)
target_link_libraries(test_iopp PRIVATE galoisring)
target_include_directories(test_iopp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_pcs tests/test_pcs.cpp)
target_link_libraries(test_pcs PRIVATE galoisring)
target_include_directories(test_pcs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(bench_pcs_commit bench/bench_pcs_commit.cpp)
target_link_libraries(bench_pcs_commit PRIVATE galoisring)

add_executable(bench_pcs_eval bench/bench_pcs_eval.cpp)
target_link_libraries(bench_pcs_eval PRIVATE galoisring)

add_executable(bench_pcs_proof_size bench/bench_pcs_proof_size.cpp)
target_link_libraries(bench_pcs_proof_size PRIVATE galoisring)

enable_testing()
add_test(NAME test_galois_ring COMMAND $<TARGET_FILE:test_galois_ring>)
add_test(NAME test_foldable_codes COMMAND $<TARGET_FILE:test_foldable_codes>)
add_test(NAME test_iopp COMMAND $<TARGET_FILE:test_iopp>)
add_test(NAME test_pcs COMMAND $<TARGET_FILE:test_pcs>)
