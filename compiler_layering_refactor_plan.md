# Compiler / PCS 分层重构方案

## 1. 背景与目标

基于论文《Polylogarithmic Proofs for Multilinears over Z2k》的实现语义，当前仓库实际承载了两步：

1. `Z_{2^k}` 上的 compiler，将问题降到 `GR(2^k, r)` 上。
2. `GR(2^k, r)` 上的底层 PCS，实现 `Setup / Commit / Eval / Verify`。

当前代码在功能上已经基本具备这两层，但结构上仍有几个明显问题：

- `compiler` 相关代码仍挂在 `BaseFold/` 目录下，概念上会误导为 “compiler 属于 BaseFold 子模块”。
- `Z2kPCSBackend` 的抽象接口和 `BaseFold` 具体适配实现耦合在一起，不利于未来并列接入其他 `GR PCS backend`。
- `src/BaseFold/BaseFoldPCS.cpp` 过于臃肿，混合了 transcript、commit、prove、verify、extension-challenge 支路和并行配置。
- `CMakeLists.txt` 只有一个 `galoisring` target，没有把层次约束固定在构建图里。

本方案的目标不是立刻大改 API，而是先把层次关系整理成长期可维护的形式：

- 明确区分 `代数层 / 公共组件层 / BaseFold PCS core / Z2k compiler`。
- 为将来接入论文第 3.2 节的 Frobenius compiler 预留稳定位置。
- 在迁移过程中尽量保持现有类型名、命名空间和协议语义稳定。
- 将“目录结构上的分层”进一步落实为“CMake target 分层”。
- 最终仓库不保留 `include/BaseFold/`，也不保留任何转发兼容头。

## 2. 设计原则

### 2.1 依赖方向

目标依赖图应为：

```text
galoisring_algebra
    ^
    |
pcs_common
    ^
    |
basefold_core
    ^
    |
z2k_compiler_core
```

约束如下：

- `BaseFold core` 不反向依赖 `compiler`。
- `compiler` 只能通过抽象 backend 接口调用底层 PCS。
- `tests/bench` 尽量直接链接最窄的 target，而不是默认拉全栈。

### 2.2 稳定性优先级

迁移过程中优先保证以下内容不变：

- `basefold` 命名空间保持不变。
- 证明对象结构和序列化语义保持不变。
- benchmark 口径保持不变。

以下内容可以在后续阶段逐步演进：

- 新目录和新 include 路径。
- 新的 CMake target。
- `BaseFoldPCS.cpp` 的文件内职责拆分。
- 仓库内 include 的批量切换。

最终态要求：

- 不保留 `include/BaseFold/` 目录。
- 不保留任何转发兼容头。
- 不保留仅用于迁移的聚合 target。

## 3. 目标目录树

建议目标目录树如下：

```text
include/
  GaloisRing/
    utils.hpp
    Inverse.hpp
    HenselLift.hpp
    PrimitiveElement.hpp

  PCS/
    Common/
      Hash.hpp
      Merkle.hpp
      Transcript.hpp
      MerkleMultiproofPlanner.hpp
      MerkleMultiproofReplay.hpp
      Multilinear.hpp
      Sumcheck.hpp
      Profile.hpp

    BaseFold/
      FoldableCode.hpp
      IOPP.hpp
      BaseFoldPCS.hpp
      ProofSerialize.hpp
      ProofSize.hpp

  Compiler/
    Z2k/
      PCSBackend.hpp
      BaseFoldBackendAdapter.hpp
      RingSwitchPCS.hpp
      RingSwitchProofSerialize.hpp
      FrobeniusPCS.hpp              # 预留，后续实现

src/
  GaloisRing/
    utils.cpp
    Inverse.cpp
    HenselLift.cpp
    PrimitiveElement.cpp

  PCS/
    Common/
      Hash.cpp
      Merkle.cpp
      Transcript.cpp
      Multilinear.cpp
      Sumcheck.cpp
      Profile.cpp

    BaseFold/
      FoldableCode.cpp
      IOPP.cpp
      BaseFoldPCSCommit.cpp
      BaseFoldPCSProve.cpp
      BaseFoldPCSVerify.cpp
      BaseFoldPCSExtension.cpp
      ProofSerialize.cpp
      ProofSize.cpp

  Compiler/
    Z2k/
      PCSBackend.cpp
      BaseFoldBackendAdapter.cpp
      RingSwitchPCS.cpp
      RingSwitchProofSerialize.cpp
      FrobeniusPCS.cpp              # 预留，后续实现
```

说明：

- `PCS/Common` 用来容纳不应被视为 BaseFold 私有的共享组件。
- `PCS/Common/Hash.hpp` 承载 `Byte / Bytes / Digest` 等字节与摘要基础类型，去掉当前 `Hash.hpp -> IOPP.hpp` 的反向依赖。
- `PCS/Common/Merkle.hpp` 承载从 `IOPP.hpp` 中拆出的通用 Merkle 类型与 API。
- `PCS/Common/Transcript.hpp` 承载通用 transcript 接口；`IOPP` 只保留与 BaseFold IOPP 调度直接相关的派生逻辑。
- `PCS/BaseFold` 承接论文第 4 节和附录 B.2/B.3 的底层 PCS。
- `Compiler/Z2k` 承接论文第 3 节的上层 compiler。
- 最终目录树中不再保留 `include/BaseFold/`。

## 4. 当前文件到目标文件的映射

### 4.1 共享公共层

| 当前路径 | 目标路径 | 说明 |
| --- | --- | --- |
| `include/BaseFold/Hash.hpp` | `include/PCS/Common/Hash.hpp` | 共享 hash 层 |
| `src/BaseFold/Hash.cpp` | `src/PCS/Common/Hash.cpp` | 共享 hash 实现 |
| `include/BaseFold/IOPP.hpp` 中的 `Byte / Bytes / Digest` | `include/PCS/Common/Hash.hpp` | 字节与摘要基础类型回收到 hash 公共层 |
| `include/BaseFold/IOPP.hpp` 中的 `MerkleRoot / MerkleTree / MerkleMultiproof*` | `include/PCS/Common/Merkle.hpp` | 从 IOPP 头中剥离通用 Merkle 抽象 |
| `src/BaseFold/IOPP.cpp` 中的通用 Merkle 实现 | `src/PCS/Common/Merkle.cpp` | Merkle commit/open/verify 实现迁移到公共层 |
| `include/BaseFold/MerkleMultiproofPlanner.hpp` | `include/PCS/Common/MerkleMultiproofPlanner.hpp` | multiproof planner 属于共享 Merkle 辅助 |
| `include/BaseFold/MerkleMultiproofReplay.hpp` | `include/PCS/Common/MerkleMultiproofReplay.hpp` | multiproof replay 属于共享 Merkle 辅助 |
| `include/BaseFold/IOPP.hpp` 中的 `FiatShamirTranscript` | `include/PCS/Common/Transcript.hpp` | 从 IOPP 头中剥离通用 transcript 接口 |
| `include/BaseFold/Multilinear.hpp` | `include/PCS/Common/Multilinear.hpp` | 共享多线性辅助 |
| `src/BaseFold/Multilinear.cpp` | `src/PCS/Common/Multilinear.cpp` | 同上 |
| `include/BaseFold/Sumcheck.hpp` | `include/PCS/Common/Sumcheck.hpp` | 共享 sumcheck；不要被 compiler 命名误导 |
| `src/BaseFold/Sumcheck.cpp` | `src/PCS/Common/Sumcheck.cpp` | 同上 |
| `include/BaseFold/Profile.hpp` | `include/PCS/Common/Profile.hpp` | 共享 profiling |
| `src/BaseFold/IOPP.cpp` 中的 `g_active_profile` 定义点 | `src/PCS/Common/Profile.cpp` | 让公共 profiling 层不再绑在 IOPP 实现文件上 |

### 4.2 BaseFold PCS core

| 当前路径 | 目标路径 | 说明 |
| --- | --- | --- |
| `include/BaseFold/FoldableCode.hpp` | `include/PCS/BaseFold/FoldableCode.hpp` | 折叠码定义 |
| `src/BaseFold/FoldableCode.cpp` | `src/PCS/BaseFold/FoldableCode.cpp` | 折叠码实现 |
| `include/BaseFold/IOPP.hpp` | `include/PCS/BaseFold/IOPP.hpp` | 仅保留 BaseFold IOPP 协议相关类型与派生逻辑 |
| `src/BaseFold/IOPP.cpp` | `src/PCS/BaseFold/IOPP.cpp` | 移除通用 Merkle/Profile 定义点后的 IOPP 实现 |
| `include/BaseFold/BaseFoldPCS.hpp` | `include/PCS/BaseFold/BaseFoldPCS.hpp` | PCS public API |
| `src/BaseFold/BaseFoldPCS.cpp` | 拆分为多个 `src/PCS/BaseFold/BaseFoldPCS*.cpp` | 按职责拆文件 |
| `include/BaseFold/ProofSerialize.hpp` | `include/PCS/BaseFold/ProofSerialize.hpp` | 精确序列化合同 |
| `src/BaseFold/ProofSerialize.cpp` | `src/PCS/BaseFold/ProofSerialize.cpp` | 当前几乎是薄壳 |
| `include/BaseFold/ProofSize.hpp` | `include/PCS/BaseFold/ProofSize.hpp` | proof-size 接口 |
| `src/BaseFold/ProofSize.cpp` | `src/PCS/BaseFold/ProofSize.cpp` | proof-size 薄包装 |

### 4.3 Z2k compiler 层

| 当前路径 | 目标路径 | 说明 |
| --- | --- | --- |
| `include/BaseFold/Z2kPCSBackend.hpp` | `include/Compiler/Z2k/PCSBackend.hpp` | 泛化 backend 接口 |
| `src/BaseFold/Z2kPCSBackend.cpp` | 拆分为 `PCSBackend.cpp` + `BaseFoldBackendAdapter.cpp` | 将抽象层和 BaseFold 适配层分离 |
| `include/BaseFold/Z2kRingSwitchPCS.hpp` | `include/Compiler/Z2k/RingSwitchPCS.hpp` | ring-switch compiler |
| `src/BaseFold/Z2kRingSwitchPCS.cpp` | `src/Compiler/Z2k/RingSwitchPCS.cpp` | compiler 实现 |
| `include/BaseFold/Z2kRingSwitchProofSerialize.hpp` | `include/Compiler/Z2k/RingSwitchProofSerialize.hpp` | compiler proof serializer |
| `src/BaseFold/Z2kRingSwitchProofSerialize.cpp` | `src/Compiler/Z2k/RingSwitchProofSerialize.cpp` | 同上 |

## 5. `BaseFoldPCS.cpp` 的建议拆分

不建议继续维持单个超大实现文件。建议按协议职责拆分，而不是按“helper 大小”随意切：

### 5.1 建议拆分后的职责

- `BaseFoldPCSCommit.cpp`
  - `BaseFoldPCSBuildCommitArtifacts*`
  - `BaseFoldPCSCommit`

- `BaseFoldPCSProve.cpp`
  - base challenge 路径的 prove 逻辑
  - `BaseFoldPCSProveEval*`
  - `BaseFoldPCSProveEvalFromCommittedTopOracle*`

- `BaseFoldPCSVerify.cpp`
  - base challenge 路径的 verify 逻辑
  - multiproof replay 和 query verification

- `BaseFoldPCSExtension.cpp`
  - `BaseFoldPCSChallengeConfig`
  - extension-challenge prove / verify 路径
  - `ZZ_pEX` 相关 helper

- `PCS/Common/Transcript.cpp`
  - `HashTranscript`
  - challenge stream
  - transcript absorb / challenge helper

### 5.2 拆分边界建议

- `HashTranscript` 不应同时在 `BaseFoldPCS.cpp` 和 `RingSwitchPCS.cpp` 各自复制一份。
- ring-switch 如果仍需独立 transcript domain，可以复用通用 transcript 基础设施，只保留自己的 domain tag 和 absorb helper。
- 环境变量解析与并行配置可以留在 `BaseFoldPCSVerify.cpp`，但不应和 transcript、extension 算术混在一起。
- `FiatShamirDeriveChallenges / FiatShamirDeriveQueryPlans` 这类依赖 `FoldableCodeParams` 与 IOPP 调度语义的 helper 继续留在 `PCS/BaseFold/IOPP.hpp|cpp`，不强行下沉到 `PCS/Common`。

## 6. 迁移期兼容策略

## 6.1 最终约束

最终仓库不保留 `include/BaseFold/`，因此兼容层不能作为最终解的一部分。

- 不能依赖长期存在的转发头来“完成”本次重构。
- 不能把 `include/BaseFold/...` 保留为 deprecated 但仍存在的公开入口。
- 迁移完成时，仓库内外公开推荐路径只应剩下新目录结构。

### 6.2 迁移期手段

如果实施过程中确实需要降低单次提交风险，可以在本地迁移分支短暂使用兼容头或兼容 target，但要满足以下约束：

- 兼容头只允许作为短期中间状态存在。
- 合并前必须全部删除。
- 最终一次完整验证必须在“无兼容头、无兼容 target”的状态下通过。

因此，更推荐的执行方式不是“先长期保留转发层”，而是：

1. 新建目标目录。
2. 一次性改写 in-tree include。
3. 直接让 tests / bench 链接新 target。
4. 在没有兼容层的状态下完成编译与测试。

### 6.3 include 重写范围

仓库内需要统一切换的 include 目标如下：

- `BaseFold/Hash.hpp` -> `PCS/Common/Hash.hpp`
- `BaseFold/MerkleMultiproofPlanner.hpp` -> `PCS/Common/MerkleMultiproofPlanner.hpp`
- `BaseFold/MerkleMultiproofReplay.hpp` -> `PCS/Common/MerkleMultiproofReplay.hpp`
- `BaseFold/Multilinear.hpp` -> `PCS/Common/Multilinear.hpp`
- `BaseFold/Sumcheck.hpp` -> `PCS/Common/Sumcheck.hpp`
- `BaseFold/Profile.hpp` -> `PCS/Common/Profile.hpp`
- `BaseFold/FoldableCode.hpp` -> `PCS/BaseFold/FoldableCode.hpp`
- `BaseFold/IOPP.hpp` -> `PCS/BaseFold/IOPP.hpp`
- `BaseFold/BaseFoldPCS.hpp` -> `PCS/BaseFold/BaseFoldPCS.hpp`
- `BaseFold/ProofSerialize.hpp` -> `PCS/BaseFold/ProofSerialize.hpp`
- `BaseFold/ProofSize.hpp` -> `PCS/BaseFold/ProofSize.hpp`
- `BaseFold/Z2kPCSBackend.hpp` -> `Compiler/Z2k/PCSBackend.hpp`
- `BaseFold/Z2kRingSwitchPCS.hpp` -> `Compiler/Z2k/RingSwitchPCS.hpp`
- `BaseFold/Z2kRingSwitchProofSerialize.hpp` -> `Compiler/Z2k/RingSwitchProofSerialize.hpp`

额外约束：

- 不能再依赖 `PCS/BaseFold/IOPP.hpp` 间接提供 `MerkleRoot / MerkleTree / MerkleMultiproof / FiatShamirTranscript / Byte / Bytes / Digest`。
- 任何只使用上述通用类型的文件，都必须显式改为包含 `PCS/Common/Hash.hpp`、`PCS/Common/Merkle.hpp`、`PCS/Common/Transcript.hpp` 中的最小必需集合。

### 6.4 命名空间策略

第一阶段建议保持 `namespace basefold` 不变，不引入命名空间层面的破坏性变更。

原因：

- 当前仓库内测试、bench、序列化、proof-size 路径都大量直接引用 `basefold::*`。
- 目录分层和 target 分层已经能解决 80% 的结构问题。
- 命名空间重写会增加大量低价值迁移噪声。

如果后续确实需要更强的 API 语义分组，可在下一轮引入：

- `basefold::pcs::basefold`
- `basefold::compiler::z2k`

但不建议和本轮目录/CMake 分层绑定一起做。

## 7. CMake target 划分

### 7.1 目标 target 图

建议把当前单一 `galoisring` target 拆成以下 target：

```text
pcs_hash_blake3
galoisring_algebra
pcs_common
basefold_core
z2k_compiler_core
```

### 7.2 每个 target 的职责

- `pcs_hash_blake3`
  - 继续承载 `third_party/blake3/*`

- `galoisring_algebra`
  - `src/GaloisRing/*.cpp`
  - `PUBLIC` 导出仓库根 `include/`
  - `PUBLIC` 链接 `NTL / GMP / m`

- `pcs_common`
  - `src/PCS/Common/Hash.cpp`
  - `src/PCS/Common/Merkle.cpp`
  - `src/PCS/Common/Transcript.cpp`
  - `src/PCS/Common/Multilinear.cpp`
  - `src/PCS/Common/Sumcheck.cpp`
  - `src/PCS/Common/Profile.cpp`
  - `PUBLIC` 导出仓库根 `include/`
  - `PUBLIC` 依赖 `galoisring_algebra` 和 `pcs_hash_blake3`
  - OpenMP compile definitions / link 也在这一层生效，因为 `Merkle.cpp` 直接使用并行 Merkle 构建

- `basefold_core`
  - `src/PCS/BaseFold/*.cpp`
  - `PUBLIC` 导出仓库根 `include/`
  - `PUBLIC` 依赖 `pcs_common`
  - 这里继续挂 OpenMP 相关 compile definitions

- `z2k_compiler_core`
  - `src/Compiler/Z2k/*.cpp`
  - `PUBLIC` 导出仓库根 `include/`
  - `PUBLIC` 依赖 `basefold_core`

### 7.3 tests / bench 的链接建议

- `test_galois_ring`
  - 链接 `galoisring_algebra`

- `test_foldable_codes`
  - 链接 `basefold_core`

- `test_iopp`
  - 链接 `basefold_core`

- `test_pcs`
  - 链接 `basefold_core`

- `test_z2k_ring_switch_pcs`
  - 链接 `z2k_compiler_core`

- `bench_pcs_commit`
  - 链接 `basefold_core`

- `bench_pcs_eval`
  - 链接 `basefold_core`

- `bench_z2k_ring_switch_commit`
  - 链接 `z2k_compiler_core`

- `bench_z2k_ring_switch_eval`
  - 链接 `z2k_compiler_core`

- `calc_iopp_params`
  - 如仅依赖参数公式，可保持独立；若后续需要共享类型，再最小化链接到 `basefold_core`

## 8. `Z2kPCSBackend` 的重构建议

当前 `Z2kPCSBackend` 同时承载了：

- 泛化 backend vtable / handle / validator；
- `BaseFoldPCS` 到该 backend 抽象的具体适配。

建议拆开：

### 8.1 `Compiler/Z2k/PCSBackend.hpp|cpp`

保留通用抽象：

- `Z2kPCSBackendVTable`
- `Z2kPCSBackendHandle`
- `Z2kPCSBackendCommitArtifacts`
- `Z2kPCSBackendEvalProof`
- `Z2kPCSBackendCommit/Prove/Verify/...`

这部分不应直接依赖 `BaseFoldPCS.hpp`，只依赖公共类型声明。

### 8.2 `Compiler/Z2k/BaseFoldBackendAdapter.hpp|cpp`

承载 BaseFold 适配：

- `MakeBaseFoldZ2kPCSBackend(const FoldableCodeParams&)`
- BaseFold 相关 `AsBaseFoldParams / AsBaseFoldCommitArtifacts / AsBaseFoldEvalProof`
- vtable 注册表 `kBaseFoldBackendVTable`

这样层次就变为：

```text
RingSwitchPCS
    -> generic PCSBackend interface
        -> BaseFoldBackendAdapter
            -> BaseFoldPCS
```

这比当前“RingSwitchPCS -> Z2kPCSBackend(内部直接知道 BaseFold)”更符合论文里的“两步结构”。

## 9. 分阶段执行方案

## Phase 0: 方案冻结

- [x] 明确目标层次和依赖方向
- [x] 写出执行文档
- [x] 决定最终目录命名采用 `PCS/` 与 `Compiler/`
- [x] 冻结 `Merkle` 与 `FiatShamirTranscript` 从 `IOPP.hpp` 抽出到 `PCS/Common`
- [x] 冻结 `MerkleMultiproofPlanner.hpp` 与 `MerkleMultiproofReplay.hpp` 进入 `PCS/Common`
- [x] 冻结 `Profile` 继续留在公共层，并把 `g_active_profile` 移到 `src/PCS/Common/Profile.cpp`

交付物：

- 本文档

Phase 0 冻结结论：

- `PCS/Common` 的最终职责明确为：`Hash / Merkle / Transcript / MerkleMultiproof helpers / Multilinear / Sumcheck / Profile`。
- `PCS/BaseFold/IOPP.hpp` 不再承担通用 Merkle 和通用 transcript 接口定义，只保留 BaseFold IOPP 协议本身的数据结构与派生流程。
- `Byte / Bytes / Digest` 统一归到 `PCS/Common/Hash.hpp`，避免公共 hash 层继续反向依赖 IOPP 头。
- `pcs_common` 必须在不依赖 `basefold_core` 的前提下独立成立，因此所有公共层定义点都不能继续挂在 `src/BaseFold/IOPP.cpp`。

## Phase 1: 新目录与 include 重写

- [x] 新建 `include/PCS/Common`
- [x] 新建 `include/PCS/BaseFold`
- [x] 新建 `include/Compiler/Z2k`
- [x] 新建对应 `src/PCS/*` 与 `src/Compiler/Z2k`
- [x] 把现有头文件移动到新路径
- [x] 批量重写仓库内旧 include 路径
- [x] 删除 `include/BaseFold/`

Phase 1 完成说明：

- 目录迁移采用“真实搬迁”而不是兼容层过渡：头文件与 `.cpp` 已直接落到 `include/PCS/*`、`include/Compiler/Z2k/*`、`src/PCS/*`、`src/Compiler/Z2k/*`。
- `Merkle` 通用类型与实现已从 `IOPP.hpp|cpp` 拆出到 `PCS/Common/Merkle.hpp|cpp`。
- `FiatShamirTranscript` 抽象接口已迁到 `PCS/Common/Transcript.hpp`；其上层的具体 hash transcript 实现将在 Phase 4 统一抽取到公共层。
- `Profile` 继续位于公共层，`g_active_profile` 已移到 `src/PCS/Common/Profile.cpp`。

验收条件：

- 仓库内源码全部使用新 include 路径
- 仓库树中不存在 `include/BaseFold/`

## Phase 2: CMake target 拆分

- [x] 引入 `galoisring_algebra`
- [x] 引入 `pcs_common`
- [x] 引入 `basefold_core`
- [x] 引入 `z2k_compiler_core`
- [x] 逐个切换 tests / bench 的链接目标

Phase 2 完成说明：

- 旧的聚合 `galoisring` target 已移除，CMake 现在直接以 `galoisring_algebra -> pcs_common -> basefold_core -> z2k_compiler_core` 表达依赖方向。
- 由于仓库内头文件路径采用 `#include "PCS/..."`、`#include "Compiler/..."`、`#include "GaloisRing/..."`，四个库 target 都 `PUBLIC` 导出仓库根 `include/`，而不是只导出各自子目录。
- `test_galois_ring` 仅链接 `galoisring_algebra`；`test_foldable_codes / test_iopp / test_pcs / bench_pcs_*` 仅链接 `basefold_core`；`test_z2k_ring_switch_pcs / bench_z2k_ring_switch_*` 仅链接 `z2k_compiler_core`。
- `calc_iopp_params` 继续保持独立可执行文件，不强行接入库 target 依赖链。
- OpenMP 相关 compile definitions / link 同时作用于 `pcs_common` 与 `basefold_core`，以覆盖 `Merkle.cpp` 和 BaseFold core 中的并行代码路径。

验收条件：

- `cmake --build build` 通过
- `ctest --test-dir build --output-on-failure` 通过
- CMake 中不存在迁移专用聚合 target

## Phase 3: compiler 层抽象清理

- [x] 拆分 `Z2kPCSBackend` 为 generic backend 和 BaseFold adapter
- [x] 更新 ring-switch include 路径和实现依赖
- [x] 保持 `RingSwitchPCSEvalProof` 和 serializer 语义不变

Phase 3 完成说明：

- `Compiler/Z2k/PCSBackend.hpp|cpp` 现在只承载 generic backend handle / vtable / 校验与调度逻辑，不再直接依赖 `PCS/BaseFold/BaseFoldPCS.hpp`。
- `Compiler/Z2k/BaseFoldBackendAdapter.hpp|cpp` 承载 `MakeBaseFoldZ2kPCSBackend(...)`、BaseFold 特定类型转换和 `kBaseFoldBackendVTable` 注册。
- `test_z2k_ring_switch_pcs` 与 `bench_z2k_ring_switch_{commit,eval}` 已改为显式通过 adapter 入口接入 BaseFold backend；`RingSwitchPCSEvalProof` 结构与 serializer 字节语义保持不变。

验收条件：

- `test_z2k_ring_switch_pcs` 全部通过
- `bench_z2k_ring_switch_commit --help`
- `bench_z2k_ring_switch_eval --help`

## Phase 4: `BaseFoldPCS.cpp` 职责拆分

- [x] 提取 transcript 基础设施
- [x] 拆分 commit / prove / verify / extension
- [x] 保持 public API 不变

Phase 4 完成说明：

- `HashTranscript`、challenge stream、按协议可配置的 domain separator / byte-order 已集中到 `PCS/Common/Transcript.hpp|cpp`；BaseFold 与 ring-switch 现在复用同一套 transcript 基础设施，只各自保留 domain 配置与 public-input absorb helper。
- 原先单文件 `src/PCS/BaseFold/BaseFoldPCS.cpp` 已拆为 `BaseFoldPCSCommon.cpp`、`BaseFoldPCSCommit.cpp`、`BaseFoldPCSProve.cpp`、`BaseFoldPCSVerify.cpp`、`BaseFoldPCSExtension.cpp`，并通过内部头 `BaseFoldPCSInternal.hpp` 共享仅限 BaseFold 内部使用的 helper。
- `include/PCS/BaseFold/BaseFoldPCS.hpp` 的 public API 保持不变；现有 `test_pcs`、`bench_pcs_commit`、`bench_pcs_eval` 调用面无需改动。

验收条件：

- `test_pcs` 全部通过
- `bench_pcs_commit --help`
- `bench_pcs_eval --help`

## Phase 5: 仓库内 include 迁移与清理

- [x] 清理迁移过程残留的旧路径和临时桥接代码
- [x] README 更新为新目录树
- [x] 新代码禁止继续引入旧路径命名习惯
- [x] 确认仓库中不存在 `include/BaseFold/`
- [x] 确认仓库中不存在兼容头或兼容聚合 target

Phase 5 当前进展补充：

- 仓库新增 `test_path_hygiene`（纳入 `ctest`），对 `include/`、`src/`、`tests/`、`bench/`、`scripts/`、`README.md`、`README_zh.md`、`CMakeLists.txt` 执行文本扫描，显式拒绝 `#include "BaseFold/..."`
  / `#include <BaseFold/...>`、`include/BaseFold/`、`src/BaseFold/` 等旧路径命名习惯。

验收条件：

- 仓库内源码默认使用新路径
- 最终树中只保留重构后的目录与 target

## 10. 推荐执行顺序

建议按下面顺序执行，而不是把目录改名、文件拆分、target 拆分混在一个提交里：

1. 先加新目录并一次性改写 include，不保留旧目录。
2. 再拆 CMake target，验证依赖图。
3. 再拆 `Z2kPCSBackend` 的抽象与适配。
4. 最后拆 `BaseFoldPCS.cpp`。

原因：

- 这样每一步都能维持“功能不变、结构更清晰”的小步前进。
- 一旦出现回归，更容易定位是目录问题、构建图问题，还是协议实现问题。
- 不会把“路径迁移噪声”和“协议逻辑变更”混在一起。

## 11. 风险与注意事项

### 11.1 不要误把 `Sumcheck` 当成 compiler 私有组件

虽然注释里提到了 “Compiler I”，但当前 `Sumcheck` 本质是共享协议组件，既服务底层 PCS，也服务 ring-switch 上层逻辑。

### 11.2 不要在第一阶段改 proof object

当前 proof-size、serializer、bench 输出已经形成稳定语义。目录和 target 重构不应顺手改 proof layout。

### 11.3 transcript 基础设施需要谨慎合并

当前 `BaseFoldPCS` 和 `RingSwitchPCS` 都有各自的 transcript domain tag。抽公共逻辑时，应共享基础设施，不应误把 domain separator 合并掉。

### 11.4 不要把兼容层写进最终交付

如果迁移过程中临时引入了兼容头或兼容 target，它们必须在最终交付前清零。否则目录分层和 target 分层都会留下“名义重构、实际兼容旧结构”的尾巴。

## 12. 最小落地版本

如果只想做一轮低风险、立即见效的整理，建议最小版本如下：

注意：这只是中间 checkpoint，不满足“最终仓库不保留 `include/BaseFold/`”的终态要求。

- 只新建 `Compiler/Z2k/`
- 只迁移以下文件：
  - `Z2kPCSBackend.*`
  - `Z2kRingSwitchPCS.*`
  - `Z2kRingSwitchProofSerialize.*`
- 先保留 `BaseFold` 其余文件在原目录
- 在同一轮中把相关 include 全部改到新路径
- 在 CMake 中先新增 `z2k_compiler_core`，其余仍可暂时保持在 `basefold_core`

这样可以先把“compiler 不属于 BaseFold 子目录”这个最大的结构问题解决掉，同时不把旧 `include/BaseFold/` 兼容层带进最终仓库，再继续做更细的公共层抽取与 `BaseFoldPCS.cpp` 拆分。

## 13. 后续执行时的检查命令

建议每个阶段至少执行以下检查：

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

对关键目标再补以下定向验证：

```bash
./build/test_pcs
./build/test_z2k_ring_switch_pcs
./build/bench_pcs_eval --help
./build/bench_z2k_ring_switch_eval --help
```

如果在 release 构建上做合并前确认，再补：

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release -j
ctest --test-dir build-release --output-on-failure
```
