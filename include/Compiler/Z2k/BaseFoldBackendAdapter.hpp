#ifndef BASEFOLD_Z2K_BASEFOLD_BACKEND_ADAPTER_HPP_
#define BASEFOLD_Z2K_BASEFOLD_BACKEND_ADAPTER_HPP_

#include "Compiler/Z2k/PCSBackend.hpp"
#include "PCS/BaseFold/BaseFoldPCS.hpp"

namespace basefold {

Z2kPCSBackendHandle MakeBaseFoldZ2kPCSBackend(
    const FoldableCodeParams &params);

}  // namespace basefold

#endif  // BASEFOLD_Z2K_BASEFOLD_BACKEND_ADAPTER_HPP_
