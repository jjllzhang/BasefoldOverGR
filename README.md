# BasefoldOverGR

基于 [NTL (Number Theory Library)](https://libntl.org/) 的一组 C++ 代码，用于在 **Galois Ring** `GR(p^k, s)`（直观上可理解为“模 `p^k` 的系数环上，再做次数为 `s` 的多项式扩张”）上实现/辅助一些常用计算：求逆、Hensel 提升、插值，以及不同表示之间的转换。

> 说明：NTL 的 `ZZ_p`/`ZZ_pX`/`ZZ_pE` 等类型依赖全局模数上下文（例如 `ZZ_p::init(mod)`、`ZZ_pE::init(F)`）。调用本项目函数前请确保相关上下文已正确初始化。

此外，本仓库还包含一份 `(c, k0, d)`-foldable linear code 的编码过程实现，支持在有限域 `F_{p^s}` 与 Galois ring `GR(p^k,s)` 上运行（见 `include/BaseFold/FoldableCode.hpp`）。该编码实现可用于 BaseFold 论文中的 IOPP 与 PCS 构造。

同时，仓库中也实现了 BaseFold 论文里的：

- BaseFold IOPP（folding + query，一并支持有限域与 Galois ring，见 `include/BaseFold/IOPP.hpp`）。
- 一个最小化的、基于 **Merkle + Fiat–Shamir** 的非交互 BaseFold PCS 单点求值证明（支持 `k0=2^κ`，见 `include/BaseFold/BaseFoldPCS.hpp`；多项式点维度为 `d+κ`）。

## 目录结构

```
.
├── include
│   └── GaloisRing
│       ├── utils.hpp
│       ├── Inverse.hpp
│       ├── HenselLift.hpp
│       └── PrimitiveElement.hpp
│   └── BaseFold
│       ├── FoldableCode.hpp
│       ├── IOPP.hpp
│       ├── Multilinear.hpp
│       ├── Sumcheck.hpp
│       ├── Profile.hpp
│       ├── ProofSize.hpp
│       └── BaseFoldPCS.hpp
├── src
│   └── GaloisRing
│       ├── utils.cpp
│       ├── Inverse.cpp
│       ├── HenselLift.cpp
│       └── PrimitiveElement.cpp
│   └── BaseFold
│       ├── FoldableCode.cpp
│       ├── IOPP.cpp
│       ├── Multilinear.cpp
│       ├── Sumcheck.cpp
│       ├── ProofSize.cpp
│       └── BaseFoldPCS.cpp
├── bench
│   ├── bench_pcs_commit.cpp
│   ├── bench_pcs_eval.cpp
│   └── bench_pcs_proof_size.cpp
└── tests
    ├── test_common.hpp
    ├── test_galois_ring_basic.cpp
    ├── test_foldable_codes.cpp
    ├── test_iopp.cpp
    └── test_pcs.cpp
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
  - 各类表示之间的转换实现（例如 `VeczzpE2Veclong`、`Long2ZZpEX` 等）。
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

### `include/BaseFold/FoldableCode.hpp` / `src/BaseFold/FoldableCode.cpp`

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

### `include/BaseFold/IOPP.hpp` / `src/BaseFold/IOPP.cpp`

- BaseFold IOPP（Protocol 2/3）实现：
  - folding prover：`ProverCommitAll/ProverCommitRound`
  - query verifier：`VerifyQueryFromOracles/VerifyQueryFromMerkleOpenings`
  - Merkle 承诺与 opening：`MerkleCommitOracle/MerkleOpenOracle/MerkleVerifyOpening`

### `include/BaseFold/Profile.hpp`

- 轻量、可选的 profiler（bench-oriented）：用于把 verifier 的时间拆分到 Merkle、折叠一致性（`EvalLineAt`）与环算术（如求逆）等模块。
- 核心接口：
  - `basefold::Profile`：累计器（ns + calls）。
  - `basefold::ProfileGuard`：在当前线程启用/关闭 profiling（基于 `thread_local`）。
  - `basefold::ScopedTimer`：RAII 计时器（仅在启用 profiling 时计时，默认开销极低）。
  - `basefold::PrintProfile`：格式化输出 profile 结果。

### `include/BaseFold/BaseFoldPCS.hpp` / `src/BaseFold/BaseFoldPCS.cpp`

- 一个最小化的非交互 BaseFold PCS 单点求值证明（Protocol 4 + Merkle + Fiat–Shamir）：
  - `BaseFoldPCSCommit/ProveEval/VerifyEval`
  - 新增配置化入口（保持旧接口不变）：`BaseFoldPCSProveEvalWithChallengeConfig` / `BaseFoldPCSProveEvalWithChallengeConfigUnchecked` / `BaseFoldPCSVerifyEvalWithChallengeConfig`
    - 当 `use_extension_challenges=false` 时，行为与旧接口一致；
    - 当 `use_extension_challenges=true` 时，Fiat–Shamir challenge 在 `ZZ_pE` 的外层扩环上采样，并把扩环参数绑定到 transcript；sumcheck / folding 在扩环内执行，同时顶层 commitment（`π_d` 的 Merkle root）仍保持在原环上。
  - 支持 `k0=2^κ` 的情况（BaseFold 论文 Remark 3）：IOPP depth 为 `d`，多项式点维度为 `d+κ`；`κ=0` 时退化为 `Basefold_over_GR.pdf` 的 Protocol 4（`k0==1`）。

### `include/BaseFold/ProofSize.hpp` / `src/BaseFold/ProofSize.cpp`

- BaseFold PCS proof size 的估算函数（bench-oriented）：
  - `basefold::BaseFoldPCSEvalProofSizeBytes(proof)`
  - `basefold::BaseFoldPCSEvalProofSizeKB(proof)`（KiB, 1024 bytes）

### `bench/bench_pcs_commit.cpp`

- 编码性能基准：测量 **去掉校验后的纯编码时间**（`EncodeFoldableUnchecked`），支持由命令行分别指定有限域与 Galois ring 的上下文参数。
- 可选 `--k0 <int>`（默认 `1`）用于测试一般 `k0` 的编码性能。
- 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

### `bench/bench_pcs_eval.cpp`

- PCS eval proof 性能基准：测量 prover time 与 verifier time（Merkle + Fiat–Shamir）。
- 默认 prover 走 `BaseFoldPCSProveEvalUnchecked`；如需把校验也算进 prover time，可加 `--checked`。
- 可加 `--profile` 输出 prover/verifier 内部耗时拆分（profile 会在 `reps` 次迭代上累加，不包含 `warmup`；建议 `--warmup 0 --reps 1` 方便阅读）。
- 可选 `--k0 <int>`（默认 `1`，要求 2 的幂）；此时多项式点维度为 `d+log2(k0)`，消息长度为 `k_d = k0*2^d`。
- 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

### `bench/bench_pcs_proof_size.cpp`

- PCS eval proof size 估算：
  - 不带 `--formula`：生成一次 **真实 proof** 并输出估算的 proof size（KB）；默认使用 `BaseFoldPCSProveEval`（包含参数/长度检查与 `claimed_y==f(z)` 校验）。
  - 带 `--formula`：完全不运行 prover，仅根据输入参数用公式给出近似/上界估算。
  - 可选 `--k0 <int>`（默认 `1`，要求 2 的幂）；此时多项式点维度为 `d+log2(k0)`，消息长度为 `k_d = k0*2^d`。
  - 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

## 依赖

- NTL（以及其底层依赖 GMP）。编译/链接方式因系统环境而异。

## 测试

项目提供如下测试：

- `tests/test_galois_ring_basic.cpp`：覆盖主要工具函数、求逆、插值，以及 Hensel 提升、`FindPrimitiveElement`、`FindTeichmullerGenerator` 的 smoke test。
- `tests/test_foldable_codes.cpp`：覆盖 foldable code 编码的正确性测试（递归编码结果与显式构造的 `G_d` 乘法结果一致）。
- `tests/test_iopp.cpp`：覆盖 BaseFold IOPP（有限域与 GR）commit/query 与 Merkle openings。
- `tests/test_pcs.cpp`：覆盖 BaseFold PCS（有限域与 GR）生成 proof 并验证通过（以及篡改后应失败）。

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

构建后可运行：

```bash
./build/bench_pcs_commit --help
./build/bench_pcs_eval --help
./build/bench_pcs_proof_size --help
```

示例（4 类用例在 commit/eval/proof_size 的对照）：

```bash
# ===== 1) 大素数域 F_p (p>1000), p=1223, zeta = -1 =====
./build/bench_pcs_commit --mode field --field-mod 1223 --field-F 0,1 --field-zeta 1222 --d 16 --reps 5 --warmup 1
./build/bench_pcs_eval --mode field --field-mod 1223 --field-F 0,1 --field-zeta 1222 --d 16 --queries 4 --reps 3 --warmup 1
./build/bench_pcs_proof_size --mode field --field-mod 1223 --field-F 0,1 --field-zeta 1222 --d 16 --queries 4

# ===== 1b) 大素数域 F_p (p>1000), p=1223, zeta != -1（auto-zeta 版本） =====
./build/bench_pcs_commit --mode field --field-mod 1223 --field-F 0,1 --auto-zeta teich --d 16 --reps 5 --warmup 1
./build/bench_pcs_eval --mode field --field-mod 1223 --field-F 0,1 --auto-zeta teich --d 16 --queries 4 --reps 3 --warmup 1
./build/bench_pcs_proof_size --mode field --field-mod 1223 --field-F 0,1 --auto-zeta teich --d 16 --queries 4

# ===== 2) 2 的扩域 F_(2^r), r>10（r=11）, zeta auto =====
# F(x)=x^11+x^2+1 -> coeffs: 1,0,1,0,0,0,0,0,0,0,0,1
./build/bench_pcs_commit --mode field --field-mod 2 --field-F 1,0,1,0,0,0,0,0,0,0,0,1 --auto-zeta teich --d 16 --reps 5 --warmup 1
./build/bench_pcs_eval --mode field --field-mod 2 --field-F 1,0,1,0,0,0,0,0,0,0,0,1 --auto-zeta teich --d 16 --queries 4 --reps 3 --warmup 1
./build/bench_pcs_proof_size --mode field --field-mod 2 --field-F 1,0,1,0,0,0,0,0,0,0,0,1 --auto-zeta teich --d 16 --queries 4

# ===== 3) GR(2^r,s), s>10（2^r=4, s=11）, zeta auto =====
./build/bench_pcs_commit --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,0,1,0,0,0,0,0,0,0,0,1 --auto-zeta teich --d 16 --reps 5 --warmup 1
./build/bench_pcs_eval --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,0,1,0,0,0,0,0,0,0,0,1 --auto-zeta teich --d 16 --queries 4 --reps 3 --warmup 1
./build/bench_pcs_proof_size --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,0,1,0,0,0,0,0,0,0,0,1 --auto-zeta teich --d 16 --queries 4

# ===== 4) 一般 GR(p^r,s), p!=2 且 p^s>1000（p=43, r=2, s=2）, zeta = -1 =====
./build/bench_pcs_commit --mode ring --ring-mod 1849 --ring-p 43 --ring-F 1,0,1 --ring-zeta 1848 --d 16 --reps 5 --warmup 1
./build/bench_pcs_eval --mode ring --ring-mod 1849 --ring-p 43 --ring-F 1,0,1 --ring-zeta 1848 --d 16 --queries 4 --reps 3 --warmup 1
./build/bench_pcs_proof_size --mode ring --ring-mod 1849 --ring-p 43 --ring-F 1,0,1 --ring-zeta 1848 --d 16 --queries 4

# ===== 4b) 一般 GR(p^r,s), p!=2 且 p^s>1000（p=43, r=2, s=2）, zeta != -1（auto-zeta 版本） =====
./build/bench_pcs_commit --mode ring --ring-mod 1849 --ring-p 43 --ring-F 1,0,1 --auto-zeta teich --d 16 --reps 5 --warmup 1
./build/bench_pcs_eval --mode ring --ring-mod 1849 --ring-p 43 --ring-F 1,0,1 --auto-zeta teich --d 16 --queries 4 --reps 3 --warmup 1
./build/bench_pcs_proof_size --mode ring --ring-mod 1849 --ring-p 43 --ring-F 1,0,1 --auto-zeta teich --d 16 --queries 4
```

### Profiling（`--profile`）

`bench_pcs_eval --profile` 会打印 prover/verifier 的 breakdown（数值用 `...` 省略）：

```text
[Ring] ...  queries=4  warmup=0 reps=1
  prover   mean ... ms
  verifier mean ... ms
  [profile-prover]
    BaseFoldPCSProveEval total:  ... ms  (calls 1)
    EncodeFoldableUnchecked:     ... ms  (calls 1)
    MerkleTree::Build:           ... ms  (calls ...)
    ...
  [profile-verifier]
    BaseFoldPCSVerifyEval total: ... ms  (calls 1)
    MerkleVerifyOpening:         ... ms  (calls ...)
    EvalLineAt:                  ... ms  (calls ...)
    ...
```

- `BaseFoldPCSProveEval total`：总 prover 时间（bench 里的 prover 段）。
- `BaseFoldPCSVerifyEval total`：总 verifier 时间。
- `MerkleTree::Build`（prover）与 `MerkleVerifyOpening`（verifier）等条目用于定位 Merkle 相关开销；`EvalLineAt`/`TryInvertUnit` 等条目用于定位折叠一致性与环算术开销。

#### 验证完整性说明（不会减少 proof 校验步骤）

- `--profile` 只是在关键函数周围打点计时，不改变任何验证逻辑。
- verifier 近期的主要加速来自工程优化：对 `FoldableCodeParams` 的合法性检查做缓存，避免在热路径上重复扫描 `diag_T`。这不会减少 Merkle opening 校验次数、folding 一致性检查次数，也不会减少 sumcheck 关系检查等协议级验证步骤。
- 该缓存假设 `FoldableCodeParams` 在构造后是不可变配置（尤其 `diag_T/zeta/G0`）。如果调用方原地修改 `params`，缓存可能导致后续不再重复触发参数错误检查；但这属于 API 误用，与 proof 验证完整性无关。

## Proof size（估算）

本仓库提供了一个 bench-oriented 的 proof size 估算函数（单位 KiB）：

- `basefold::BaseFoldPCSEvalProofSizeKB(proof)`（见 `include/BaseFold/ProofSize.hpp`）
  也可以直接使用 `bench/bench_pcs_proof_size.cpp`：
  - 不带 `--formula`：生成真实 proof 并输出估算大小；
  - 带 `--formula`：仅参数估计（不运行 prover）。
