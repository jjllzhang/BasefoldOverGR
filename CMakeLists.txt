cmake_minimum_required(VERSION 3.16)

project(BasefoldOverGR LANGUAGES CXX)

add_library(galoisring
  src/GaloisRing/utils.cpp
  src/GaloisRing/Inverse.cpp
  src/GaloisRing/HenselLift.cpp
  src/GaloisRing/PrimitiveElement.cpp
  src/BaseFold/FoldableCode.cpp
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

target_link_libraries(galoisring PUBLIC
  ${NTL_LIBRARY}
  ${GMP_LIBRARY}
)
if(M_LIBRARY)
  target_link_libraries(galoisring PUBLIC ${M_LIBRARY})
endif()

add_executable(test_galois_ring tests/test_galois_ring_basic.cpp)
target_link_libraries(test_galois_ring PRIVATE galoisring)

add_executable(test_foldable_codes tests/test_foldable_codes.cpp)
target_link_libraries(test_foldable_codes PRIVATE galoisring)

enable_testing()
add_test(NAME test_galois_ring COMMAND $<TARGET_FILE:test_galois_ring>)
add_test(NAME test_foldable_codes COMMAND $<TARGET_FILE:test_foldable_codes>)
