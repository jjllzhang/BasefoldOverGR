# BasefoldOverGR

语言版本：

- 英文：`README.md`
- 中文：`README_zh.md`

基于 [NTL (Number Theory Library)](https://libntl.org/) 的一组 C++ 代码，用于在 **Galois Ring** `GR(p^k, s)`（直观上可理解为“模 `p^k` 的系数环上，再做次数为 `s` 的多项式扩张”）上实现/辅助一些常用计算：求逆、Hensel 提升、插值，以及不同表示之间的转换。

> 说明：NTL 的 `ZZ_p`/`ZZ_pX`/`ZZ_pE` 等类型依赖全局模数上下文（例如 `ZZ_p::init(mod)`、`ZZ_pE::init(F)`）。调用本项目函数前请确保相关上下文已正确初始化。

此外，本仓库还包含一份 `(c, k0, d)`-foldable linear code 的编码过程实现，支持在有限域 `F_{p^s}` 与 Galois ring `GR(p^k,s)` 上运行（见 `include/PCS/BaseFold/FoldableCode.hpp`）。该编码实现可用于 BaseFold 论文中的 IOPP 与 PCS 构造。

同时，仓库中也实现了 BaseFold 论文里的：

- BaseFold IOPP（folding + query，一并支持有限域与 Galois ring，见 `include/PCS/BaseFold/IOPP.hpp`）。
- 一个最小化的、基于 **Merkle + Fiat–Shamir** 的非交互 BaseFold PCS 单点求值证明（支持 `k0=2^κ`，见 `include/PCS/BaseFold/BaseFoldPCS.hpp`；多项式点维度为 `d+κ`）。

## 目录结构

```
.
├── CMakeLists.txt
├── README.md
├── README_zh.md
├── FRI_Ligero-based_results.md
├── LICENSE
├── bench
│   ├── bench_basefold_pcs_commit.cpp
│   ├── bench_basefold_pcs_eval.cpp
│   ├── bench_z2k_frobenius_commit.cpp
│   ├── bench_z2k_frobenius_eval.cpp
│   ├── calc_iopp_params.cpp
│   └── exp_params_release_c4_lambda128.md
├── include
│   ├── GaloisRing
│   │   ├── utils.hpp
│   │   ├── Inverse.hpp
│   │   ├── FrobeniusBasis.hpp
│   │   ├── HenselLift.hpp
│   │   └── PrimitiveElement.hpp
│   ├── PCS
│   │   ├── Common
│   │   │   ├── Hash.hpp
│   │   │   ├── Merkle.hpp
│   │   │   ├── Multilinear.hpp
│   │   │   ├── Profile.hpp
│   │   │   ├── Sumcheck.hpp
│   │   │   └── Transcript.hpp
│   │   └── BaseFold
│   │       ├── FoldableCode.hpp
│   │       ├── IOPP.hpp
│   │       ├── ProofSerialize.hpp
│   │       ├── ProofSize.hpp
│   │       └── BaseFoldPCS.hpp
│   └── Compiler
│       └── Z2k
│           ├── BaseFoldBackendAdapter.hpp
│           ├── FrobeniusPCS.hpp
│           ├── FrobeniusProofSerialize.hpp
│           ├── PCSBackend.hpp
│           ├── RingSwitchPCS.hpp
│           └── RingSwitchProofSerialize.hpp
├── scripts
│   ├── run_release_c4_lambda128.sh
│   └── plot_benchmark_results.py
├── src
│   ├── GaloisRing
│   │   ├── utils.cpp
│   │   ├── Inverse.cpp
│   │   ├── FrobeniusBasis.cpp
│   │   ├── HenselLift.cpp
│   │   └── PrimitiveElement.cpp
│   ├── PCS
│   │   ├── Common
│   │   │   ├── Hash.cpp
│   │   │   ├── Merkle.cpp
│   │   │   ├── Multilinear.cpp
│   │   │   ├── Profile.cpp
│   │   │   ├── Sumcheck.cpp
│   │   │   └── Transcript.cpp
│   │   └── BaseFold
│   │       ├── FoldableCode.cpp
│   │       ├── IOPP.cpp
│   │       ├── ProofSerialize.cpp
│   │       ├── ProofSize.cpp
│   │       ├── BaseFoldPCSCommon.cpp
│   │       ├── BaseFoldPCSCommit.cpp
│   │       ├── BaseFoldPCSProve.cpp
│   │       ├── BaseFoldPCSVerify.cpp
│   │       └── BaseFoldPCSExtension.cpp
│   └── Compiler
│       └── Z2k
│           ├── BaseFoldBackendAdapter.cpp
│           ├── FrobeniusPCS.cpp
│           ├── FrobeniusProofSerialize.cpp
│           ├── PCSBackend.cpp
│           ├── RingSwitchPCS.cpp
│           └── RingSwitchProofSerialize.cpp
├── tests
│   ├── test_common.hpp
│   ├── test_galois_ring_basic.cpp
│   ├── test_foldable_codes.cpp
│   ├── test_iopp.cpp
│   ├── test_z2k_frobenius_feasibility.cpp
│   ├── test_z2k_frobenius_bench_cli.cpp
│   ├── test_z2k_frobenius_pcs.cpp
│   └── test_pcs.cpp
├── result
│   ├── results-*.csv
│   └── plots/
└── results  # 运行 sweep 脚本后生成（默认不纳入仓库）
    ├── release_c4_lambda128_sweep_<RUN_ID>/  # sweep 脚本默认输出目录
    └── single_runs/
```

## 文件说明

### `include/GaloisRing/utils.hpp`

- 通用工具函数的声明与依赖汇总，主要包括：
  - `vector<long>` 与 NTL 多项式/扩张元素之间的相互转换（`ZZ_pX`、`ZZ_pE`、`ZZ_pEX`、`vec_ZZ_pE`）。
  - 向量工具：补零、裁剪、分段拆分等。
  - 数学辅助：判断/计算 2 的幂（`isPowerOfTwo`、`nextPowerOf2` 等）、因子枚举（`FindFactor`）、最近完全平方数等。
  - `FindPrimitivePoly`：在 `Z_p[x]` 中寻找次数为 `n` 的 primitive/不可约多项式（基于 NTL 的 `BuildIrred` 并加额外判定）。
  - `interpolate_for_GR`：在 `GR(p^l,s)` 上做插值；当 `l>1`（不再是域）时会依赖 `Inverse.hpp` 中的求逆逻辑。

### `src/GaloisRing/utils.cpp`

- `utils.hpp` 中工具函数的具体实现：
  - 各类表示之间的转换实现（例如 `FlattenZZpEVectorToLongs`、`LongVecToZZpEX` 等）。
  - `interpolate_for_GR`：当 `l==1` 直接调用 NTL 的 `interpolate`；当 `l>1` 使用自定义 `Inv()` 来处理环上的“可逆元求逆”，完成插值过程。

### `include/GaloisRing/Inverse.hpp`

- 求逆相关接口声明：
  - `ZZ_pE Inv(ZZ_pE a, long s)`：在给定扩张次数/表示维度 `s` 下求 `a` 的逆（不可逆时返回 0）。
  - `ZZ_pE Inv2(ZZ_pE a, ZZ_pX F, ZZ p, long s, long k)`：另一套基于 `p`-进分解/提升的求逆实现（用于模 `p^k` 的场景）。

### `src/GaloisRing/Inverse.cpp`

- `Inv` 与 `Inv2` 的实现：
  - `Inv`：通过构造线性系统（矩阵）解出逆元表示，并对不可逆情况捕获 `NTL::InvModErrorObject` 后返回 0。
  - `Inv2`：将多项式系数按 `p`-进展开，迭代提升到模 `p^k` 的逆（内部使用 `XGCD`、模/除/乘等辅助函数）。

### `include/GaloisRing/HenselLift.hpp`

- Hensel 提升函数声明：
  - `HenselLift(g_, f, g, p, n)`：已知 `g | f (mod p)`，将因子 `g` 提升到更高模数（迭代到 `p^{n+1}` 量级）得到 `g_`。

### `src/GaloisRing/HenselLift.cpp`

- Hensel 提升实现：
  - 内部一步提升：同时更新 `g,h,s,t`，满足新的同余分解关系。
  - 外层循环：逐步将模数从 `p` 提升到 `p^2 ... p^{n+1}`，输出最终提升后的 `g_`。

### `include/GaloisRing/PrimitiveElement.hpp`

- “本原元素/生成元”相关接口声明：
  - `FindPrimitiveElement(ZZ p, long k, long s, const ZZ_pX& F)`：在 `GR(p^k,s)` 上（由扩张多项式 `F` 指定）构造一个候选的本原元素。
  - `FindTeichmullerGenerator(ZZ p, long k, long s, const ZZ_pX& F, long max_trials=1024)`：在 `GR(p^k,s)` 上寻找 Teichmüller 子群（阶 `p^s-1`）的生成元。

### `src/GaloisRing/PrimitiveElement.cpp`

- `FindPrimitiveElement` 的实现：
  - 由调用者传入扩张多项式 `F`（应为次数 `s` 且首项系数为 1 的多项式），再用 `ZZ_pE::init(F)` 初始化扩张。
  - 以 `b = x (mod F)` 为基础返回 `b^(p^{k-1})`（常见于 Teichmüller 代表/单位根相关构造）。
  - 注意：`F` 的选择会直接影响 `GR(p^k,s)` 的结构；实际使用时应传入与你的 `p,s` 匹配的（模 `p` 约化后不可约的）多项式。
- `FindTeichmullerGenerator` 的实现：
  - 先尝试确定性候选 `x^(p^{k-1})`，若阶已是 `p^s-1` 则直接返回。
  - 若不满足，则随机采样单位元 `u`，投影 `u^(p^{k-1})` 到 Teichmüller 子群后做阶检验，直到找到生成元或达到 `max_trials`。

### `include/PCS/BaseFold/FoldableCode.hpp` / `src/PCS/BaseFold/FoldableCode.cpp`

- Foldable code 的编码实现：
  - 参数结构 `basefold::FoldableCodeParams`：`c, k0, d, p, zeta, G0, diag_T`。
    - `G0`：一个 `[n0=ck0, k0]` 的 MDS 线性码生成矩阵（实现中不强制验证 MDS 性质，仅检查维度）。
    - `diag_T[i]`：对角矩阵 `Ti` 的对角元（向量表示），长度必须为 `n_i = c*k0*2^i`。在域上要求每一项非零；在 Galois ring 上要求每一项是单位元（unit）。
    - `zeta`：固定的 `ζ`，在域上要求 `ζ != 0,1`；在 Galois ring 上要求 `ζ` 与 `(1-ζ)` 都是单位元（等价于 `π(ζ) != 1`）。
    - `p`：可选的素数 `p`（用于在 Galois ring `GR(p^s,r)` 场景下做 unit 校验；若不设置则跳过 unit 校验，仅做非零检查）。
  - 编码接口 `basefold::EncodeFoldable(out, msg, params)`：输入 `msg ∈ R^{k_d}`，输出 `out ∈ R^{n_d}`，其中 `k_d = k0*2^d`、`n_d = c*k_d`。
  - 额外提供 `basefold::EncodeFoldableUnchecked(...)`：用于 bench/热路径，**跳过参数与长度校验**（调用者需自行保证参数正确）。
  - NTL 上下文前置条件：
    - 域 `F_{p^r}`：先 `ZZ_p::init(p)`，再用次数为 `r` 的不可约多项式初始化 `ZZ_pE::init(F)`；
    - Galois ring `GR(p^s,r)`：先 `ZZ_p::init(p^s)`，再用“模 p 约化后不可约”的次数为 `r` 多项式初始化 `ZZ_pE::init(F)`。

### `include/PCS/BaseFold/IOPP.hpp` / `src/PCS/BaseFold/IOPP.cpp`

- BaseFold IOPP（Protocol 2/3）实现：
  - folding prover：`ProverCommitAll/ProverCommitRound`
  - query verifier：`VerifyQueryFromOracles`
  - Merkle 承诺与 multiproof：`MerkleCommitOracle/MerkleOpenOracleMany/MerkleVerifyMultiproof`

### `include/PCS/Common/Profile.hpp`

- 轻量、可选的 profiler（bench-oriented）：用于把 verifier 的时间拆分到 Merkle、折叠一致性（`EvalLineAt`）与环算术（如求逆）等模块。
- 核心接口：
  - `basefold::Profile`：累计器（ns + calls）。
  - `basefold::ProfileGuard`：在当前线程启用/关闭 profiling（基于 `thread_local`）。
  - `basefold::ScopedTimer`：RAII 计时器（仅在启用 profiling 时计时，默认开销极低）。
  - `basefold::PrintProfile`：格式化输出 profile 结果。

### `include/PCS/BaseFold/BaseFoldPCS.hpp` / `src/PCS/BaseFold/BaseFoldPCS*.cpp`

- 一个最小化的非交互 BaseFold PCS 单点求值证明（Protocol 4 + Merkle + Fiat–Shamir）：
  - `BaseFoldPCSCommit/ProveEval/VerifyEval`
  - 新增配置化入口（保持旧接口不变）：`BaseFoldPCSProveEvalWithChallengeConfig` / `BaseFoldPCSProveEvalWithChallengeConfigUnchecked` / `BaseFoldPCSVerifyEvalWithChallengeConfig`
    - 当 `use_extension_challenges=false` 时，行为与旧接口一致；
    - 当 `use_extension_challenges=true` 时，Fiat–Shamir challenge 在 `ZZ_pE` 的外层扩环上采样，并把扩环参数绑定到 transcript；sumcheck / folding 在扩环内执行，并对扩环中间层 `π_0..π_{d-1}` 做可验证 Merkle 承诺，同时顶层 commitment（`π_d` 的 Merkle root）仍保持在原环上。
    - 扩环 challenge 路径的 proof payload 已收敛为 multiproof-only 且做了去冗余压缩：不再重复携带 base 侧 `h_i / pi0_codeword`，base 与 extension query payload 都使用共享的 Merkle multiproof；`extension.r_by_level` 允许省略（由 transcript 重采样恢复）。
    - `challenge_extension_modulus` 现在会校验代数条件：在域模式要求不可约；在环模式要求模 `p` 约化后不可约（basic irreducible）。
  - 支持 `k0=2^κ` 的情况（BaseFold 论文 Remark 3）：IOPP depth 为 `d`，多项式点维度为 `d+κ`；`κ=0` 时退化为 BaseFold 论文中的 Protocol 4（`k0==1`）。
  - 实现现已按职责拆分：
    - `BaseFoldPCSCommit.cpp`：顶层 commit artifacts 与 `BaseFoldPCSCommit`
    - `BaseFoldPCSProve.cpp`：base-challenge prove 路径
    - `BaseFoldPCSVerify.cpp`：base-challenge verify 路径与 verifier query 并行配置
    - `BaseFoldPCSExtension.cpp`：extension-challenge prove / verify 路径
    - `BaseFoldPCSCommon.cpp`：仅供 BaseFold 内部复用的 shared helper

### `include/PCS/BaseFold/ProofSerialize.hpp` / `src/PCS/BaseFold/ProofSerialize.cpp`

- fixed-width proof serializer 的公共契约：
  - `basefold::FixedProofEncodingOptions`
  - `basefold::FixedProofEncodingContext`
  - `basefold::CountingSink`
  - `basefold::CountSerializedBaseFoldPCSEvalProofFixedBytes(...)`

### `include/PCS/BaseFold/ProofSize.hpp` / `src/PCS/BaseFold/ProofSize.cpp`

- 通过 fixed-width proof serializer 做 BaseFold PCS proof size 的精确计数：
  - `basefold::BaseFoldPCSEvalProofSizeBytes(proof)`
  - `basefold::BaseFoldPCSEvalProofSizeKB(proof)`（KiB, 1024 bytes）
- 计数时会省略 verifier 可从 transcript 重建的 query indices / 显式扩环
  challenge（例如 `extension.r_by_level`）。
- 面向用户的 proof size 输出来自 `bench_basefold_pcs_eval`：同一次 prove/eval bench
  会直接打印 `proof_size_bytes` / `proof_size_kb`。

### `include/Compiler/Z2k/FrobeniusPCS.hpp` / `src/Compiler/Z2k/FrobeniusPCS.cpp`

- 基于 Frobenius map 的 `Z_{2^k} -> GR(2^k, r)` compiler（论文 `Protocol 2`），
  构建在现有通用 backend 边界之上。
- 对外接口沿用 ring-switch 那套分层：
  - setup/packing/commit helper，
  - 可复用 commit artifacts，
  - outer-only proof（`FrobeniusPCSOuterEvalProof`），
  - composed proof（`FrobeniusPCSEvalProof`），
  - staged prover/verifier wrapper。
- setup 现在支持两种模式：
  - 在 setup 阶段自动搜索 normal/dual basis，
  - 通过 `FrobeniusPCSProvidedBasisInput` 显式传入论文语义下的
    `beta/alpha` basis 数据。
- 显式 basis 路径仍然会做强校验。若调用方不提供 Teichmüller generator，
  setup 会先自行构造一个，再检查给定的 `beta/alpha` 是否按真实
  Frobenius orbit 顺序排列。
- 当前范围仍然是 correctness-first 的单点评估证明；不包含 multipoint
  opening、Blaze-Orion 风格 proof composition、或 WHIR。

使用显式 basis 的最小 setup 示例：

```cpp
basefold::FrobeniusPCSSetupInput input;
input.ell = ell;
input.kappa = kappa;
input.base_modulus = base_modulus;
input.extension_modulus = extension_modulus;
input.use_provided_basis = true;
input.provided_basis.normal_basis = normal_basis;  // beta + alpha
input.provided_basis.has_teichmuller_generator = false;
input.backend = backend;

const basefold::FrobeniusPCSParams params = basefold::FrobeniusPCSSetup(input);
```

如果你已经有可信的 Teichmüller generator，可以把
`input.provided_basis.has_teichmuller_generator = true`，并填入
`input.provided_basis.teichmuller_generator`。

### `include/Compiler/Z2k/FrobeniusProofSerialize.hpp` / `src/Compiler/Z2k/FrobeniusProofSerialize.cpp`

- 为公开的 Frobenius proof object 提供 fixed-width serializer 和 proof-size helper。
- serializer 契约：
  - outer-only size 只统计 `s_by_i`、`h_by_level`、`t_star`，
  - composed size 在 outer proof 后追加一个显式 8-byte backend proof 长度前缀和 backend proof bytes，
  - commitment 与 public inputs 不计入 proof-size reporting。
- proof size 的“精确”含义是：相对于这套 fixed-width serializer 契约精确；size helper 和 Frobenius bench 共用同一条计数路径。

### `bench/bench_basefold_pcs_commit.cpp`

- 顶层 commit 基准：
  - `encode-only mean`：顶层原始编码时间（`EncodeFoldableUnchecked`，不含校验）
  - `top-commit mean`：对顶层 oracle 做 `MerkleTree::Build` 并提取 root 的时间
  - `commit mean`：headline 的顶层 commit 阶段 = encode + 顶层 Merkle commit
- 有限域与 Galois ring 上下文都可由命令行分别配置。
- 可选 `--k0 <int>`（默认 `1`）用于测试一般 `k0` 的编码性能。
- 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

### `bench/bench_z2k_frobenius_commit.cpp`

- Frobenius compiler 的 commit benchmark。
- 在进入 timed region 之前，bench setup 会先做一次有界、确定性的
  normal-basis 候选搜索，并把这次运行固定到选中的 preferred basis；
  搜索结果会打印到输出里，但不计入 headline timing。
- `packing mean` 统计 `t -> t'` packing 加上 Boolean-table 到 monomial
  conversion 的时间，具体走 `FrobeniusPCSBuildOuterCommitArtifacts(...)`。
- `commit mean` 统计完整的 `FrobeniusPCSCommit(...)` 路径。

### `bench/bench_basefold_pcs_eval.cpp`

- PCS eval proof 性能基准：测量 prove-phase time 与 verifier time；其中 prove-phase 从预先构造好的顶层 commitment artifacts 开始，不包含顶层 encode + commit。
- 默认 prover 走 `BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked`；如需把校验也算进 prove-phase time，可加 `--checked`。
- 同一次运行会通过 fixed-width `CountingSink` 路径输出该 proof 的精确 payload 大小
  （`proof_size_bytes` / `proof_size_kb`，省略 verifier 可自行从 transcript 重建的
  indices/challenges）。
- 可加 `--use-extension-challenges` 切到扩域/扩环 challenge 路径（`BaseFoldPCSChallengeConfig`），此时 sumcheck/folding 在扩域/扩环上运行，顶层 commitment 仍在原域/环。
- 可选 `--field-challenge-ext` / `--ring-challenge-ext` 指定 challenge 扩域模多项式 `E(U)`，格式为 `a0;a1;...;ad`（每个 `ai` 是一个 `ZZ_pE` 元素，写作 `c0,c1,...`）；若不指定，默认 `E(U)=zeta + U + U^2`。
- 可加 `--profile` 输出 prover/verifier 内部耗时拆分（profile 会在 `reps` 次迭代上累加，不包含 `warmup`；建议 `--warmup 0 --reps 1` 方便阅读）。
- 可选 `--k0 <int>`（默认 `1`，要求 2 的幂）；此时多项式点维度为 `d+log2(k0)`，消息长度为 `k_d = k0*2^d`。
- 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

### `bench/bench_z2k_frobenius_eval.cpp`

- Frobenius compiler 的 eval benchmark：测量 outer protocol 加 backend
  单点评估证明的 prove/verify 开销。
- 和其他 Frobenius bench 二进制一样，它会在 timed region 外先固定到
  这次预搜索选出的 preferred normal basis。
- 计时 prove 从预先构造好的 `FrobeniusPCSCommitArtifacts` 开始，和当前
  `bench_basefold_pcs_eval` 一样不把顶层 commit 阶段算进 prove-phase。
- 输出内容包括：
  - 总 prover/verifier 时间，
  - 扣掉 backend 子调用后的 outer prover/verifier 时间，
  - 对公开 `FrobeniusPCSEvalProof` 做 serializer-backed 计数得到的
    `outer proof size` 与 `proof size`。

### `scripts/plot_benchmark_results.py`

- 从一个或多个 benchmark CSV 画出随 `d` 变化的曲线图（commit/prover/verifier/proof size）。
- 默认输出目录 `result/plots/`，默认输出前缀 `benchmark`。
- 默认读取 `proof_size_kb` 列作为 proof-size 曲线；可用 `--proof-size-column` 改成其他列。
- 可用 `--metrics` 选择只画部分指标（`commit / prover / verifier / proof_size / all`）。

## 依赖

- NTL（以及其底层依赖 GMP）。编译/链接方式因系统环境而异。

## 测试

项目提供如下测试：

- `tests/test_galois_ring_basic.cpp`：覆盖主要工具函数、求逆、插值，以及 Hensel 提升、`FindPrimitiveElement`、`FindTeichmullerGenerator` 的 smoke test。
- `tests/test_foldable_codes.cpp`：覆盖 foldable code 编码的正确性测试（递归编码结果与显式构造的 `G_d` 乘法结果一致）。
- `tests/test_iopp.cpp`：覆盖 BaseFold IOPP（有限域与 GR）commit/query 与 Merkle multiproof。
- `tests/test_pcs.cpp`：覆盖 BaseFold PCS（有限域与 GR）生成 proof 并验证通过（以及篡改后应失败）。
- `tests/test_z2k_frobenius_feasibility.cpp`：覆盖小参数 Galois ring 上的 Frobenius 代数可行性回归。
- `tests/test_z2k_frobenius_bench_cli.cpp`：覆盖 Frobenius commit/eval bench CLI 的 `--help`、smoke execution，以及稳定 proof-size 输出字段。
- `tests/test_z2k_frobenius_pcs.cpp`：覆盖 Frobenius PCS 的 setup/packing/prove/verify、provided-basis setup 校验，以及 serializer-backed proof-size 检查。

在安装好 NTL/GMP 后，使用 CMake（推荐 out-of-source 构建）：

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Bench

如需更稳定/可对比的性能数据，建议使用 Release 构建（Debug/未优化构建用于 correctness 测试即可）：

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

哈希固定使用 `BLAKE3`：

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
```

构建后可运行：

```bash
./build-release/bench_basefold_pcs_commit --help
./build-release/bench_basefold_pcs_eval --help
./build-release/bench_z2k_frobenius_commit --help
./build-release/bench_z2k_frobenius_eval --help
./build-release/dump_basefold_pcs_eval_artifact --help
./build-release/bench_basefold_pcs_verify_artifact --help
./build-release/calc_iopp_params --help
```

### 两种 verifier-only 模式

现在有两条 verifier-only benchmark 路径：

- `bench_basefold_pcs_verify`：自包含 verifier benchmark。命令本身会在进程内按
  确定性逻辑生成 proof，然后对该 proof 重复做 verify 并计时。
- `bench_basefold_pcs_verify_artifact`：artifact 驱动的 verifier benchmark。命令先从
  磁盘加载一个已落盘的 proof case，在 timed region 外完成反序列化，再只对
  已加载的 proof 重复做 verify 并计时。

如果你想一次命令里复现完整 logical case，用 `bench_basefold_pcs_verify`；如果你想把
文件 IO 和 proof/public-input 反序列化排除在 verifier mean 之外，用
`bench_basefold_pcs_verify_artifact`。

### Artifact 工作流

artifact benchmarking 是增量能力，不会替代现有 `bench_basefold_pcs_*` 二进制。

`dump_basefold_pcs_eval_artifact` 每次调用只落一个 BaseFold eval-proof case：

```bash
./build-release/dump_basefold_pcs_eval_artifact \
  --artifact-root /tmp/basefold-artifacts \
  --artifact-id field_d10_q2 \
  --mode field --d 10 --queries 2 --seed 5
```

artifact root 的布局为：

```text
/tmp/basefold-artifacts/
  manifest.jsonl
  objects/
    <artifact_id>/
      meta.json
      public_inputs.bin
      proof.bin
```

artifact 语义要点：

- 每次调用只处理一个 case，`--mode both` 会被拒绝。
- 已存在的 `artifact_id` 默认不会被覆盖。
- `manifest.jsonl` 是主索引面，object 目录名保持简短。
- `meta.json` 保存 verifier 需要的关键元数据，包括该 case 实际使用的
  `zeta_coeffs`。
- `public_inputs.bin` 只保存 verifier public inputs：`commitment_root`、
  点 `z`、claimed value `y`。
- artifact 不会保存完整顶层 commit artifacts，也不会保存 `MerkleTree`
  内部结构。

然后用 artifact 路径只测 pure verify：

```bash
./build-release/bench_basefold_pcs_verify_artifact \
  --artifact-root /tmp/basefold-artifacts \
  --artifact-id field_d10_q2 \
  --warmup 0 --reps 3
```

artifact-driven verifier timing 的语义是：

- `artifact load wall time` 与 `artifact deserialize wall time` 会作为诊断信息打印，
  但不计入 headline `verifier mean`。
- `input proof size` 是从磁盘载入的 proof 的精确 fixed-width 大小，使用和
  `bench_basefold_pcs_eval` 相同的 serializer contract。
- 只有 `--verifier-query-*` 会影响 timed verify loop。
- `--merkle-*` 在这里被明确禁止，因为 artifact verify 路径不会在计时段内
  构建 Merkle tree。

如需和自包含 verifier benchmark 对比同一 logical case，可运行：

```bash
./build-release/bench_basefold_pcs_verify \
  --mode field --d 10 --queries 2 --warmup 0 --reps 3 --seed 5
```

### 一键复现实验（`c=4`, `lambda=128`）

仓库内置 sweep 脚本 `scripts/run_release_c4_lambda128.sh`，默认会在 14 组上下文上遍历 `d=3..29`，并输出到 `results/release_c4_lambda128_sweep_<RUN_ID>/`（`RUN_ID` 默认形如 `<timestamp>_pid<pid>`）。

```bash
scripts/run_release_c4_lambda128.sh

# 只跑单个 context + 较小维度区间
CONTEXTS=field-prime128-ext D_MIN=10 D_MAX=20 scripts/run_release_c4_lambda128.sh
```

参数与上下文说明见 `bench/exp_params_release_c4_lambda128.md`。

默认输出目录里会包含：

- `results.csv`：每个 `(context, d)` 的结构化结果汇总。
- `RESULTS.md`：便于快速浏览的 markdown 表格。
- `logs/*.log`：`calc_iopp_params`、`bench_basefold_pcs_commit`、`bench_basefold_pcs_eval` 的原始日志。

`results.csv` 当前表头为：

```text
context_id,context_label,mode,d,poly_dim,c,k0,lambda,gamma,queries,commit_mean_ms,prove_phase_mean_ms,verifier_mean_ms,proof_size_kb,proof_size_bytes,status,error
```

其中常用列含义：

- `poly_dim`：消息/真值表长度 `k_d = k0*2^d`（默认 `k0=1` 时就是 `2^d`）。
- `gamma`：`calc_iopp_params --auto-gamma` 选出的 slack 参数。
- `queries`：推荐查询次数（即 `l_min_for_PCS`）。
- `proof_size_kb/proof_size_bytes`：来自 `bench_basefold_pcs_eval` 的精确 payload 大小
  （fixed-width counting 路径，省略 verifier 可从 transcript 重建的
  indices/challenges，KiB / bytes）。
- 当前 release sweep 不会再从单独的 formula-only benchmark 二进制填这两列。
- `status,error`：该 `(context,d)` 点是否成功及失败原因（若有）。

### 结果文件与绘图

- 单次 profile 拆分（如 `--profile --reps 1`）可放在 `results/single_runs/*_profile_breakdown.md`（例如 `results/single_runs/<context>_d<d>_profile_breakdown.md`）。
- 当前跟踪的对比 CSV 位于 `results-legacy/results-*.csv`（例如 `results-legacy/results-F_2^256.csv`）。
- 使用 `scripts/plot_benchmark_results.py` 从 CSV 生成图：

```bash
# 单个 CSV 出图
python3 scripts/plot_benchmark_results.py results-legacy/results-F_2^256.csv --prefix f2_256

# 多个 CSV 叠加对比
python3 scripts/plot_benchmark_results.py \
  results-legacy/results-F_2^256.csv \
  'results-legacy/results-GR(2^16;128).csv' \
  --prefix compare_f2_256_vs_gr2p16_128

# 直接使用某次 sweep 的 results.csv
python3 scripts/plot_benchmark_results.py \
  results/release_c4_lambda128_sweep_<RUN_ID>/results.csv \
  --prefix sweep_run
```

### 参数选取工具（码距 + 推荐 query 次数）

`bench/calc_iopp_params.cpp` 输入 `d/c/lambda` 与 `q`（或 `p,r,m`）后，自动计算：

- 距离下界 `Delta_Cd >= 1 - t_d / n_d`（按 Theorem 1 / Corollary 1 的 `t_i` 递推）；
- `delta < J_gamma(J_gamma(Delta_Cd))`（严格不等式）；
- 实现策略：固定 `gamma` 时，在满足全部约束（`delta < J_gamma(J_gamma(Delta_Cd))`、`0 < 1-delta+gamma*d < 1`、`3*delta-gamma*d < Delta_Cd`）的可行区间内，取尽量接近上界的 `delta`，以减小查询数；
- 计算并区分：
  - `l_min_iopp_only`：满足 `2d/(gamma^3 q) + (1-delta+gamma*d)^l <= 2^-lambda`；
  - `l_min_for_PCS`（推荐）：满足 `2d/q + 2d/(gamma^3 q) + (1-delta+gamma*d)^l <= 2^-lambda`。
- 支持 `--auto-gamma`：在给定 `c,d,k0,lambda,q` 时自动搜索 `gamma`，目标最小化 `l_min_for_PCS`。
- 搜索实现：`Delta_Cd` 先预计算一次，再对 `gamma` 进行粗到细自适应搜索。

示例：

```bash
# 直接给 q（示例：q = 2^192）
./build-release/calc_iopp_params --d 20 --c 16 --lambda 128 --q 6277101735386680763835789423207666416102355444464034512896 --gamma 0.005

# 按 q = p^(r*m) 给（示例：p=2, r=64, m=3 => q=2^192）
./build-release/calc_iopp_params --d 20 --c 16 --lambda 128 --p 2 --r 64 --m 3 --gamma 0.005

# 查看每一层递推细节（ell_i, t_i）
./build-release/calc_iopp_params --d 16 --c 16 --lambda 128 --q 6277101735386680763835789423207666416102355444464034512896 --gamma 0.00625 --show-levels

# 自动选 gamma（目标：l_min_for_PCS 最小）
./build-release/calc_iopp_params --d 20 --c 16 --k0 1 --lambda 128 --q 6277101735386680763835789423207666416102355444464034512896 --auto-gamma

# 可选：限制搜索区间和搜索预算
./build-release/calc_iopp_params --d 20 --c 16 --k0 1 --lambda 128 --q 6277101735386680763835789423207666416102355444464034512896 --auto-gamma --gamma-min 1e-4 --gamma-max 0.05 --gamma-steps 8000
```

示例（128-bit 可复现参数，一套常量覆盖 release 面向用户的 bench）：

```bash
# ---------- 固定参数（可选；仅在你想用变量写法时需要） ----------
# 128-bit 素数（field）
FIELD_MOD_128=326594724262804054738278293730872375507
# 64-bit 素数 p，以及 p^2（128-bit，ring）
RING_P_64=18446744073709551557
RING_MOD_128=340282366920938461286658806734041124249   # = RING_P_64^2
#
# field: F(x)=x^2+1, zeta=x
FIELD_F=1,0,1
FIELD_ZETA=0,1
#
# ring:  F(x)=x^2+x+1, zeta=x
RING_F=1,1,1
RING_ZETA=0,1

# ---------- Field profile (GF(p^2), p 为 128-bit；无变量版本，直接可运行) ----------
./build-release/bench_basefold_pcs_commit --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --d 10 --reps 3 --warmup 1
./build-release/bench_basefold_pcs_eval --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --d 10 --queries 2 --reps 2 --warmup 1

# ---------- Ring profile (GR(p^2,2), p 为 64-bit, p^2 为 128-bit；无变量版本，直接可运行) ----------
./build-release/bench_basefold_pcs_commit --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --d 10 --reps 3 --warmup 1
./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --d 10 --queries 2 --reps 2 --warmup 1

# ---------- extension-challenge 路径 ----------
# 注意：challenge 多项式参数含 ';'，请使用引号
# E(U) = (0 + 3*x) + U + U^2  => '0,3;1;1'
./build-release/bench_basefold_pcs_eval --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --use-extension-challenges --field-challenge-ext '0,3;1;1' --d 10 --queries 2 --reps 1 --warmup 0
./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --use-extension-challenges --ring-challenge-ext '0,3;1;1' --d 10 --queries 2 --reps 1 --warmup 0
```

### Merkle Build 并行阈值调优

`MerkleTree::Build` 支持通过环境变量或 `bench_basefold_pcs_eval` CLI 调整并行策略：

```bash
# 环境变量（对所有 bench 生效）
export BASEFOLD_MERKLE_LEAVES_PER_THREAD=32768
export BASEFOLD_MERKLE_PARALLEL_LEVEL_THRESHOLD=4096
export BASEFOLD_MERKLE_MAX_THREADS=8

# CLI（仅 bench_basefold_pcs_eval，优先级高于环境变量）
./build-release/bench_basefold_pcs_eval ... --merkle-leaves-per-thread 32768 --merkle-level-threshold 4096 --merkle-max-threads 8
```

自动 sweep（示例）：

```bash
csv=/tmp/merkle_threshold_sweep.csv
echo "threshold,prove_phase_mean_ms,merkle_build_total_ms" > "$csv"
for t in 256 512 1024 2048 4096 8192 16384 32768 65536; do
  out=$(./build-release/bench_basefold_pcs_eval --mode field \
    --field-mod 326594724262804054738278293730872375507 \
    --field-F 1,0,1 --field-zeta 0,1 \
    --d 14 --queries 4 --warmup 1 --reps 2 --profile \
    --merkle-level-threshold "$t")
  prove_phase=$(printf '%s\n' "$out" | awk '/prove-phase mean/{print $3; exit}')
  merkle=$(printf '%s\n' "$out" | awk '/MerkleTree::Build:/{print $2; exit}')
  echo "$t,$prove_phase,$merkle" | tee -a "$csv"
done
cat "$csv"
```

### Verifier Query 并行线程调优

`BaseFoldPCSVerifyEval` 的 query-level 并行支持环境变量或 `bench_basefold_pcs_eval` CLI：

```bash
# 环境变量（对 bench_basefold_pcs_eval 生效）
export BASEFOLD_VERIFY_QUERY_QUERIES_PER_THREAD=1
export BASEFOLD_VERIFY_QUERY_PARALLEL_THRESHOLD=2
export BASEFOLD_VERIFY_QUERY_MAX_THREADS=8

# CLI（优先级高于环境变量）
./build-release/bench_basefold_pcs_eval ... \
  --verifier-query-per-thread 1 \
  --verifier-query-threshold 2 \
  --verifier-query-max-threads 8
```

自动 sweep 最优线程数（示例）：

```bash
csv=/tmp/verifier_query_threads_sweep.csv
echo "max_threads,verifier_mean_ms,prove_phase_mean_ms" > "$csv"
for t in 1 2 4 8 12 16 24 32; do
  out=$(./build-release/bench_basefold_pcs_eval --mode field \
    --field-mod 326594724262804054738278293730872375507 \
    --field-F 1,0,0,0,1 --field-zeta 3 \
    --d 14 --queries 64 --warmup 1 --reps 3 \
    --verifier-query-per-thread 1 \
    --verifier-query-threshold 2 \
    --verifier-query-max-threads "$t")
  verifier=$(printf '%s\n' "$out" | awk '/verifier mean/{print $3; exit}')
  prove_phase=$(printf '%s\n' "$out" | awk '/prove-phase mean/{print $3; exit}')
  echo "$t,$verifier,$prove_phase" | tee -a "$csv"
done
cat "$csv"
```

### Profiling（`--profile`）

`bench_basefold_pcs_eval --profile` 会打印 prover/verifier 的 breakdown（数值用 `...` 省略）：

```text
[Ring] ...  queries=4  warmup=0 reps=1
  prove-phase mean ... ms
  verifier mean ... ms
  [profile-prover]
    BaseFoldPCSProveEval total:  ... ms  (calls 1)
    MerkleTree::Build:           ... ms  (calls ...)
    ...
  [profile-verifier]
    BaseFoldPCSVerifyEval total: ... ms  (calls 1)
    MerkleVerifyMultiproof:      ... ms  (calls ...)
    EvalLineAt:                  ... ms  (calls ...)
    ...
```

- `BaseFoldPCSProveEval total`：总 prover 时间（bench 里的 prover 段）。
- `BaseFoldPCSVerifyEval total`：总 verifier 时间。
- `MerkleTree::Build`（prover）与 `MerkleVerifyMultiproof`（verifier）等条目用于定位 Merkle 相关开销；`EvalLineAt`/`TryInvertUnit` 等条目用于定位折叠一致性与环算术开销。

#### 验证完整性说明（不会减少 proof 校验步骤）

- `--profile` 只是在关键函数周围打点计时，不改变任何验证逻辑。
- verifier 近期的主要加速来自工程优化：对 `FoldableCodeParams` 的合法性检查做缓存，避免在热路径上重复扫描 `diag_T`。这不会减少 Merkle opening 校验次数、folding 一致性检查次数，也不会减少 sumcheck 关系检查等协议级验证步骤。
- 该缓存假设 `FoldableCodeParams` 在构造后是不可变配置（尤其 `diag_T/zeta/G0`）。如果调用方原地修改 `params`，缓存可能导致后续不再重复触发参数错误检查；但这属于 API 误用，与 proof 验证完整性无关。

## Proof size（精确值）

本仓库通过 fixed-width proof serializer 提供 BaseFold PCS proof size 的精确计数：

- `basefold::BaseFoldPCSEvalProofSizeBytes(proof)`（bytes）
- `basefold::BaseFoldPCSEvalProofSizeKB(proof)`（KiB）

对终端用户而言，推荐直接看 `bench_basefold_pcs_eval`：它在生成真实 proof 后，会在同一次运行里输出
`proof_size_bytes` / `proof_size_kb`。

如果走 verifier-only artifact benchmark，则 `bench_basefold_pcs_verify_artifact`
会对磁盘中载入的 proof 输出同样精确的 `input proof size`，同时仍把文件 IO 与
反序列化排除在 headline verifier timing 之外。

## License

本项目采用 MIT License，详见 `LICENSE`。
