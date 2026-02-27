# BasefoldOverGR

语言版本：

- 英文：`README.md`
- 中文：`README_zh.md`

基于 [NTL (Number Theory Library)](https://libntl.org/) 的一组 C++ 代码，用于在 **Galois Ring** `GR(p^k, s)`（直观上可理解为“模 `p^k` 的系数环上，再做次数为 `s` 的多项式扩张”）上实现/辅助一些常用计算：求逆、Hensel 提升、插值，以及不同表示之间的转换。

> 说明：NTL 的 `ZZ_p`/`ZZ_pX`/`ZZ_pE` 等类型依赖全局模数上下文（例如 `ZZ_p::init(mod)`、`ZZ_pE::init(F)`）。调用本项目函数前请确保相关上下文已正确初始化。

此外，本仓库还包含一份 `(c, k0, d)`-foldable linear code 的编码过程实现，支持在有限域 `F_{p^s}` 与 Galois ring `GR(p^k,s)` 上运行（见 `include/BaseFold/FoldableCode.hpp`）。该编码实现可用于 BaseFold 论文中的 IOPP 与 PCS 构造。

同时，仓库中也实现了 BaseFold 论文里的：

- BaseFold IOPP（folding + query，一并支持有限域与 Galois ring，见 `include/BaseFold/IOPP.hpp`）。
- 一个最小化的、基于 **Merkle + Fiat–Shamir** 的非交互 BaseFold PCS 单点求值证明（支持 `k0=2^κ`，见 `include/BaseFold/BaseFoldPCS.hpp`；多项式点维度为 `d+κ`）。

## 目录结构

```
.
├── CMakeLists.txt
├── README.md
├── README_zh.md
├── FRI_Ligero-based_results.md
├── LICENSE
├── bench
│   ├── bench_pcs_commit.cpp
│   ├── bench_pcs_eval.cpp
│   ├── bench_pcs_proof_size.cpp
│   ├── bench_pcs_communication.cpp
│   ├── calc_iopp_params.cpp
│   └── exp_params_release_c4_lambda128.md
├── include
│   ├── GaloisRing
│   │   ├── utils.hpp
│   │   ├── Inverse.hpp
│   │   ├── HenselLift.hpp
│   │   └── PrimitiveElement.hpp
│   └── BaseFold
│       ├── FoldableCode.hpp
│       ├── IOPP.hpp
│       ├── Multilinear.hpp
│       ├── Sumcheck.hpp
│       ├── Profile.hpp
│       ├── ProofSize.hpp
│       └── BaseFoldPCS.hpp
├── scripts
│   ├── run_release_c4_lambda128.sh
│   └── plot_benchmark_results.py
├── src
│   ├── GaloisRing
│   │   ├── utils.cpp
│   │   ├── Inverse.cpp
│   │   ├── HenselLift.cpp
│   │   └── PrimitiveElement.cpp
│   └── BaseFold
│       ├── FoldableCode.cpp
│       ├── IOPP.cpp
│       ├── Multilinear.cpp
│       ├── Sumcheck.cpp
│       ├── ProofSize.cpp
│       └── BaseFoldPCS.cpp
├── tests
│   ├── test_common.hpp
│   ├── test_galois_ring_basic.cpp
│   ├── test_foldable_codes.cpp
│   ├── test_iopp.cpp
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
    - 当 `use_extension_challenges=true` 时，Fiat–Shamir challenge 在 `ZZ_pE` 的外层扩环上采样，并把扩环参数绑定到 transcript；sumcheck / folding 在扩环内执行，并对扩环中间层 `π_0..π_{d-1}` 做可验证 Merkle 承诺，同时顶层 commitment（`π_d` 的 Merkle root）仍保持在原环上。
    - 扩环 challenge 路径的 proof payload 已做去冗余压缩：不再重复携带 base 侧 `h_i / pi0_full`，`query_proofs` 只保留顶层 base opening；`extension.r_by_level` 允许省略（由 transcript 重采样恢复）。
    - `challenge_extension_modulus` 现在会校验代数条件：在域模式要求不可约；在环模式要求模 `p` 约化后不可约（basic irreducible）。
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
- 可加 `--use-extension-challenges` 切到扩域/扩环 challenge 路径（`BaseFoldPCSChallengeConfig`），此时 sumcheck/folding 在扩域/扩环上运行，顶层 commitment 仍在原域/环。
- 可选 `--field-challenge-ext` / `--ring-challenge-ext` 指定 challenge 扩域模多项式 `E(U)`，格式为 `a0;a1;...;ad`（每个 `ai` 是一个 `ZZ_pE` 元素，写作 `c0,c1,...`）；若不指定，默认 `E(U)=zeta + U + U^2`。
- 可加 `--profile` 输出 prover/verifier 内部耗时拆分（profile 会在 `reps` 次迭代上累加，不包含 `warmup`；建议 `--warmup 0 --reps 1` 方便阅读）。
- 可选 `--k0 <int>`（默认 `1`，要求 2 的幂）；此时多项式点维度为 `d+log2(k0)`，消息长度为 `k_d = k0*2^d`。
- 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

### `bench/bench_pcs_proof_size.cpp`

- PCS eval proof size 估算：
  - 不带 `--formula`：生成一次 **真实 proof** 并输出估算的 proof size（KB）；默认使用 `BaseFoldPCSProveEval`（包含参数/长度检查与 `claimed_y==f(z)` 校验）。
  - 带 `--formula`：完全不运行 prover，仅根据输入参数用公式给出近似/上界估算。
  - 可加 `--use-extension-challenges` 切到扩域/扩环 challenge 路径：
    - 非公式模式会调用 `BaseFoldPCSProveEvalWithChallengeConfig`，统计真实扩域/扩环 proof payload；
    - 公式模式会按当前紧凑化后的 `BaseFoldPCSEvalProof` 字段布局计入 payload（base 顶层 commitment/root 与顶层 openings + extension `roots/h/msg0/pi0/query openings`；`r_i` 默认不单独携带）。
  - 可选 `--field-challenge-ext` / `--ring-challenge-ext` 指定 challenge 扩域模多项式 `E(U)`，格式为 `a0;a1;...;ad`（每个 `ai` 是一个 `ZZ_pE` 元素，写作 `c0,c1,...`）；若不指定，默认 `E(U)=zeta + U + U^2`。
  - 可选 `--field-challenge-degree <m>` / `--ring-challenge-degree <m>` 仅指定扩张次数，自动构造默认多项式 `E(U)=zeta + U + U^m`（`m=1` 时为 `E(U)=zeta+U`）。
  - `--*-challenge-ext` 与 `--*-challenge-degree` 互斥。
  - 可选 `--k0 <int>`（默认 `1`，要求 2 的幂）；此时多项式点维度为 `d+log2(k0)`，消息长度为 `k_d = k0*2^d`。
  - 可选 `--auto-zeta teich` 自动从 `(p,k,F)` 推导 Teichmüller 子群生成元作为 `zeta`（启用后忽略 `--field-zeta/--ring-zeta`）。

### `bench/bench_pcs_communication.cpp`

- PCS prover/verifier 通信量估算（仅公式，不运行 prover）：
  - Base challenge 路径：`P -> V` 按当前 `BaseFoldPCSEvalProof` payload 拆分（roots、sumcheck、`pi0_full`、query openings）。
  - 扩域/扩环 challenge 路径（`--use-extension-challenges`）：
    - `P -> V` 估算按当前紧凑 payload：base 顶层 root + 顶层 base openings + extension `roots/h/msg0/pi0/query openings`；
    - 默认不单独计入 extension `r_i`（按 transcript 重算）。
  - `V -> P`：给出交互式等价口径（`d` 个 challenge `r_i` + `queries` 个索引 `mu`）。
  - 同时输出当前 Fiat–Shamir 非交互路径总通信量（`V -> P = 0`）。
  - 可选 `--field-challenge-ext` / `--ring-challenge-ext` 指定 challenge 扩域模 `E(U)`（省略时默认 `E(U)=zeta + U + U^2`）。
  - 可选 `--field-challenge-degree <m>` / `--ring-challenge-degree <m>` 仅指定扩张次数，自动构造默认多项式 `E(U)=zeta + U + U^m`（`m=1` 时为 `E(U)=zeta+U`）。
  - `--*-challenge-ext` 与 `--*-challenge-degree` 互斥。
  - 可选 `--k0 <int>`（默认 `1`，要求 2 的幂）。

### `scripts/plot_benchmark_results.py`

- 从一个或多个 benchmark CSV 画出随 `d` 变化的曲线图（commit/prover/verifier/proof size）。
- 默认输出目录 `result/plots/`，默认输出前缀 `benchmark`。
- 默认读取 `proof_size_kb` 列作为 proof-size 曲线；可用 `--proof-size-column`（或别名 `--communication-column`）改成其他列（例如通信量列）。
- 可用 `--metrics` 选择只画部分指标（`commit / prover / verifier / proof_size / all`）。

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
./build-release/bench_pcs_commit --help
./build-release/bench_pcs_eval --help
./build-release/bench_pcs_proof_size --help
./build-release/bench_pcs_communication --help
./build-release/calc_iopp_params --help
```

### 一键复现实验（`c=4`, `lambda=128`）

仓库内置 sweep 脚本 `scripts/run_release_c4_lambda128.sh`，默认会在 12 组上下文上遍历 `d=3..29`，并输出到 `results/release_c4_lambda128_sweep_<RUN_ID>/`（`RUN_ID` 默认形如 `<timestamp>_pid<pid>`）。

```bash
scripts/run_release_c4_lambda128.sh

# 只跑单个 context + 较小维度区间
CONTEXTS=field-prime128-ext D_MIN=10 D_MAX=20 scripts/run_release_c4_lambda128.sh
```

参数与上下文说明见 `bench/exp_params_release_c4_lambda128.md`。

默认输出目录里会包含：

- `results.csv`：每个 `(context, d)` 的结构化结果汇总。
- `RESULTS.md`：便于快速浏览的 markdown 表格。
- `logs/*.log`：`calc_iopp_params` / `bench_*` 的原始日志。

`results.csv` 当前表头为：

```text
context_id,context_label,mode,d,poly_dim,c,k0,lambda,gamma,queries,commit_mean_ms,prover_mean_ms,verifier_mean_ms,proof_size_kb,proof_size_bytes,status,error
```

其中常用列含义：

- `gamma`：`calc_iopp_params --auto-gamma` 选出的 slack 参数。
- `queries`：推荐查询次数（即 `l_min_for_PCS`）。
- `proof_size_kb/proof_size_bytes`：来自 `bench_pcs_proof_size --formula` 的估算结果（KiB / bytes）。
- `status,error`：该 `(context,d)` 点是否成功及失败原因（若有）。

### 结果文件与绘图

- 单次 profile 拆分（如 `--profile --reps 1`）可放在 `results/single_runs/*_profile_breakdown.md`（例如 `results/single_runs/ring-gr-2p16-64-ext_d15_profile_breakdown.md`）。
- 手工整理/对比用 CSV 可放在 `result/results-*.csv`（例如 `result/results-F_2^256.csv`）。
- 使用 `scripts/plot_benchmark_results.py` 从 CSV 生成图：

```bash
# 单个 CSV 出图
python3 scripts/plot_benchmark_results.py result/results-F_2^256.csv --prefix f2_256

# 多个 CSV 叠加对比
python3 scripts/plot_benchmark_results.py \
  result/results-F_2^256.csv \
  'result/results-GR(2^16;128).csv' \
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

示例（128-bit 可复现参数，一套常量覆盖四个 bench）：

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
./build-release/bench_pcs_commit --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --d 10 --reps 3 --warmup 1
./build-release/bench_pcs_eval --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --d 10 --queries 2 --reps 2 --warmup 1
./build-release/bench_pcs_proof_size --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --d 10 --queries 2 --formula
./build-release/bench_pcs_communication --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --d 10 --queries 2

# ---------- Ring profile (GR(p^2,2), p 为 64-bit, p^2 为 128-bit；无变量版本，直接可运行) ----------
./build-release/bench_pcs_commit --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --d 10 --reps 3 --warmup 1
./build-release/bench_pcs_eval --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --d 10 --queries 2 --reps 2 --warmup 1
./build-release/bench_pcs_proof_size --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --d 10 --queries 2 --formula
./build-release/bench_pcs_communication --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --d 10 --queries 2

# ---------- extension-challenge 路径 ----------
# 注意：challenge 多项式参数含 ';'，请使用引号
# E(U) = (0 + 3*x) + U + U^2  => '0,3;1;1'
./build-release/bench_pcs_eval --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --use-extension-challenges --field-challenge-ext '0,3;1;1' --d 10 --queries 2 --reps 1 --warmup 0
./build-release/bench_pcs_eval --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --use-extension-challenges --ring-challenge-ext '0,3;1;1' --d 10 --queries 2 --reps 1 --warmup 0
./build-release/bench_pcs_proof_size --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --field-zeta 0,1 --use-extension-challenges --field-challenge-ext '0,3;1;1' --d 10 --queries 2 --formula
./build-release/bench_pcs_proof_size --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --ring-zeta 0,1 --use-extension-challenges --ring-challenge-ext '0,3;1;1' --d 10 --queries 2 --formula
```

### Merkle Build 并行阈值调优

`MerkleTree::Build` 支持通过环境变量或 `bench_pcs_eval` CLI 调整并行策略：

```bash
# 环境变量（对所有 bench 生效）
export BASEFOLD_MERKLE_LEAFS_PER_THREAD=32768
export BASEFOLD_MERKLE_PARALLEL_LEVEL_THRESHOLD=4096
export BASEFOLD_MERKLE_MAX_THREADS=8

# CLI（仅 bench_pcs_eval，优先级高于环境变量）
./build-release/bench_pcs_eval ... --merkle-leafs-per-thread 32768 --merkle-level-threshold 4096 --merkle-max-threads 8
```

自动 sweep（示例）：

```bash
csv=/tmp/merkle_threshold_sweep.csv
echo "threshold,prover_mean_ms,merkle_build_total_ms" > "$csv"
for t in 256 512 1024 2048 4096 8192 16384 32768 65536; do
  out=$(./build-release/bench_pcs_eval --mode field \
    --field-mod 326594724262804054738278293730872375507 \
    --field-F 1,0,1 --field-zeta 0,1 \
    --d 14 --queries 4 --warmup 1 --reps 2 --profile \
    --merkle-level-threshold "$t")
  prover=$(printf '%s\n' "$out" | awk '/prover   mean/{print $3; exit}')
  merkle=$(printf '%s\n' "$out" | awk '/MerkleTree::Build:/{print $2; exit}')
  echo "$t,$prover,$merkle" | tee -a "$csv"
done
cat "$csv"
```

### Verifier Query 并行线程调优

`BaseFoldPCSVerifyEval` 的 query-level 并行支持环境变量或 `bench_pcs_eval` CLI：

```bash
# 环境变量（对 bench_pcs_eval 生效）
export BASEFOLD_VERIFY_QUERY_QUERIES_PER_THREAD=1
export BASEFOLD_VERIFY_QUERY_PARALLEL_THRESHOLD=2
export BASEFOLD_VERIFY_QUERY_MAX_THREADS=8

# CLI（优先级高于环境变量）
./build-release/bench_pcs_eval ... \
  --verifier-query-per-thread 1 \
  --verifier-query-threshold 2 \
  --verifier-query-max-threads 8
```

自动 sweep 最优线程数（示例）：

```bash
csv=/tmp/verifier_query_threads_sweep.csv
echo "max_threads,verifier_mean_ms,prover_mean_ms" > "$csv"
for t in 1 2 4 8 12 16 24 32; do
  out=$(./build-release/bench_pcs_eval --mode field \
    --field-mod 326594724262804054738278293730872375507 \
    --field-F 1,0,0,0,1 --field-zeta 3 \
    --d 14 --queries 64 --warmup 1 --reps 3 \
    --verifier-query-per-thread 1 \
    --verifier-query-threshold 2 \
    --verifier-query-max-threads "$t")
  verifier=$(printf '%s\n' "$out" | awk '/verifier mean/{print $3; exit}')
  prover=$(printf '%s\n' "$out" | awk '/prover   mean/{print $3; exit}')
  echo "$t,$verifier,$prover" | tee -a "$csv"
done
cat "$csv"
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

## Communication（估算）

可使用 `bench/bench_pcs_communication.cpp` 仅根据输入参数估算 PCS 双向通信量：

- Base challenge 路径：
  - `P -> V`：proof payload（roots、sumcheck、`pi0_full`、Merkle openings）。
- 扩域/扩环 challenge 路径（`--use-extension-challenges`）：
  - `P -> V`：按紧凑 payload 估算（base 顶层 root/openings + extension `roots/h/msg0/pi0/query openings`）。
  - 默认不单独计入 extension `r_i`（按 transcript 重算）。
- challenge 扩域可用两种方式指定：
  - `--field/ring-challenge-ext`：显式给 `E(U)` 多项式系数。
  - `--field/ring-challenge-degree m`：仅给扩张次数，自动构造默认 `E(U)=zeta + U + U^m`（`m=1` 时为 `zeta+U`）。
- `V -> P`：交互式等价挑战（`r_i` 与 `mu`）。
- 同时给出当前 Fiat–Shamir 路径通信量（`V -> P = 0`）。

示例：

```bash
./build-release/bench_pcs_communication --mode field --field-mod 326594724262804054738278293730872375507 --field-F 1,0,1 --d 10 --queries 2
./build-release/bench_pcs_communication --mode ring --ring-mod 340282366920938461286658806734041124249 --ring-p 18446744073709551557 --ring-F 1,1,1 --d 10 --queries 2
./build-release/bench_pcs_communication --mode field --field-mod 2 --field-F 1,1,1 --use-extension-challenges --field-challenge-ext '0,1;1;1' --d 10 --queries 2
./build-release/bench_pcs_communication --mode field --field-mod 2 --field-F 1,1,1 --use-extension-challenges --field-challenge-degree 4 --d 10 --queries 2
```

## License

本项目采用 MIT License，详见 `LICENSE`。
