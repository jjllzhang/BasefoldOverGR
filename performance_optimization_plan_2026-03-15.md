# BasefoldOverGR 性能优化计划（2026-03-15）

## 目标

在不改变当前论文语义、Fiat-Shamir 调度、Merkle 承诺定义、proof 序列化契约的前提下，优先做原型实现里“高收益、低风险、易验证”的性能优化。

这个计划只覆盖现在值得做的优化，不包含生产级工程化、复杂底层重写、或协议边界附近的高风险改动。

## 测量与约束

- 统一使用 `build-release/` 做性能判断；`build/` 只用于 correctness 调试。
- 区分 `bench_basefold_pcs_commit` 与 `bench_basefold_pcs_eval`：
  - `commit` 看 top encode + top Merkle。
  - `eval` 的 `prove-phase` 不包含 top commit。
- 单点定位优先用：

```bash
./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --queries 4 --profile --warmup 0 --reps 1
./build-release/bench_basefold_pcs_commit --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --warmup 0 --reps 1
./build-release/bench_z2k_frobenius_eval --c 4 --ell 10 --kappa 1 --queries 4 --warmup 0 --reps 1 --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1
```

- 每一项优化都必须满足：
  - tests 通过；
  - bench 输出语义不变；
  - proof bytes 不因“偷删协议负载”而变化，除非明确同步修改 verifier/serializer 契约并单独审计。

## 不做的事

- 不改协议流程，不改 challenge label / 吸收顺序，不改 query 派生规则。
- 不改 Merkle root 定义，不改当前 proof 格式契约。
- 不做生产级 allocator、SIMD 手写、复杂 task graph 并行。
- 不继续深挖 `Inv()` 微优化，除非新的 release profile 证明它重新进入主热点。
- 不优先做 transcript 字符串拼接、小对象拷贝、零散 glue code 抛光类优化。

## 优先级 A

### A1. Frobenius / RingSwitch outer prover 增加 unchecked 热路径

状态：`completed`（2026-03-16）

目的：
- 去掉 prover 热路径里仅用于 honest-witness 自检的重复重算。
- 保留现有 checked API；新增 unchecked/bench-friendly 路径，而不是删除检查。

当前热点依据：
- Frobenius outer prover 已明显大于 backend prover。
- 现有 outer prove 内部包含额外直接求值、partial reconstruction、一致性自检。

主要范围：
- `src/Compiler/Z2k/FrobeniusPCS.cpp`
- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `include/Compiler/Z2k/FrobeniusPCS.hpp`
- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- 对应 benches/tests

拟做内容：
- 为 outer prove / full prove 增加 unchecked 入口，语义参考 BaseFold 现有 checked/unchecked 分层。
- 将下列 prover 侧自检保留在 checked 路径：
  - `claimed_s == t(z)`
  - `recovered_partials == direct_partials`
  - honest Equality Check 1/3
- bench 默认优先走 unchecked 热路径，和 BaseFold eval 的当前语义保持一致。

验收：
- `tests/test_z2k_frobenius_pcs.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`
- `tests/test_z2k_frobenius_bench_cli.cpp`
- `tests/test_z2k_ring_switch_bench_cli.cpp`
- 前后对比 `bench_z2k_frobenius_eval` 的 outer prover mean。

已完成内容：
- 新增 Frobenius / RingSwitch 的 checked / unchecked prove API 分层：
  - `ProveOuterEvalUnchecked(...)`
  - `ProveOuterEvalFromCommitArtifactsUnchecked(...)`
  - `ProveEvalUnchecked(...)`
  - `ProveEvalFromCommitArtifactsUnchecked(...)`
- checked 路径保留 prover 侧 honest-witness 自检；unchecked 路径跳过这些仅用于本地防呆的重算。
- `Prover` / `OuterProver` class 默认仍走 checked API，不静默改变普通调用语义。
- `bench_z2k_frobenius_eval`
- `bench_z2k_ring_switch_eval`
- `bench_z2k_frobenius_outer_prove`
- `bench_z2k_ring_switch_outer_prove`
  默认改为走 unchecked；新增 `--checked` 显式切回 checked 路径。
- 补充 checked/unchecked 一致性回归测试与 CLI `--checked` 覆盖。

本轮验证：
- `cmake --build build-release -j 4 --target bench_z2k_frobenius_eval bench_z2k_ring_switch_eval bench_z2k_frobenius_outer_prove bench_z2k_ring_switch_outer_prove test_z2k_frobenius_bench_cli test_z2k_ring_switch_bench_cli test_z2k_frobenius_pcs test_z2k_ring_switch_pcs`
- `ctest --test-dir build-release --output-on-failure -R 'test_z2k_(frobenius|ring_switch)_(bench_cli|pcs)'`

### A2. 共享并重写 equality table 构造为迭代 DP

状态：`completed`（2026-03-16）

目的：
- 消除 `BuildEqualityTable(point)` 中“每个布尔点都重新构造 `BooleanPointFromIndex` 再调用 `EqPolynomial`”的低效写法。
- 用一个共享 helper 覆盖 Frobenius 与 RingSwitch。

主要范围：
- `src/Compiler/Z2k/FrobeniusPCS.cpp`
- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- 必要时抽到 `src/PCS/Common/` 或 `include/PCS/Common/`

拟做内容：
- 提供迭代式 equality table 构造：
  - 从长度 `1` 的表开始；
  - 每加入一个变量，把旧表扩成两半；
  - 避免 `2^d` 次 `EqPolynomial(...)` 调用。
- 替换下列使用点：
  - prefix equality table
  - Frobenius suffix orbit 上的 equality tables
  - RingSwitch prefix equality table

验收：
- Frobenius / RingSwitch 相关 tests 全绿。
- `bench_z2k_frobenius_eval` outer prover mean 下降。
- 如有必要，再补一个 equality-table helper 单测。

已完成内容：
- 在共享层新增 `EqualityTableFromPoint(...)`，放到 `PCS/Common/Multilinear`。
- 用迭代 DP 直接构造 `eq_point` 在 Boolean hypercube 上的 evaluation table，不再对每个布尔点单独构造向量并调用 `EqPolynomial(...)`。
- Frobenius suffix orbit 上的 `sigma_eq_tables_by_i` 改为走共享 helper。
- Frobenius prefix equality table `lambda_by_i` 改为走共享 helper。
- RingSwitch prefix equality table 改为走共享 helper。
- RingSwitch prove / verify 热路径中，同一轮 prefix equality table 的重复构造收敛为一次。
- 新增共享 helper 单测，保持与 `EqPolynomial` 定义逐点一致。

本轮验证：
- `cmake --build build-release -j 4 --target test_pcs test_z2k_frobenius_pcs test_z2k_ring_switch_pcs bench_z2k_frobenius_eval bench_z2k_ring_switch_eval bench_z2k_frobenius_outer_prove bench_z2k_ring_switch_outer_prove`
- `ctest --test-dir build-release --output-on-failure -R 'test_(pcs|z2k_frobenius_pcs|z2k_ring_switch_pcs)'`
- 对照基准：
  - A1 baseline: commit `d7453cc`
  - A2 current workspace
  - 统一参数：`OMP_NUM_THREADS=1`, `c=4`, `ell=10`, `kappa=1`, `queries=4`, `warmup=1`, `reps=5`, `GR(4,2)`

本轮测得收益：
- `bench_z2k_frobenius_eval` outer prover mean：`7.438 ms -> 1.718 ms`
- `bench_z2k_frobenius_outer_prove`：`7.165 ms -> 1.284 ms`
- `bench_z2k_ring_switch_eval` outer prover mean：`4.597 ms -> 4.472 ms`
- `bench_z2k_ring_switch_outer_prove`：`4.341 ms -> 4.470 ms`

备注：
- Frobenius 收益很大，说明此前热点确实主要在 equality table 的逐点重算。
- RingSwitch 收益较小，但 `eval` 仍有轻微下降；`outer_prove` 在这组小参数下基本视为持平，后续应在更大 `kappa` / `ell_prime` 上复测再决定是否继续深挖 RingSwitch 外层。

### A3. verifier multiproof 查找表一次性构建

状态：`completed`（2026-03-16）

目的：
- 避免 verifier 在 query 内层循环中反复对同一 `queried_indices` 做 `lower_bound`。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSVerify.cpp`
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`
- `include/PCS/Common/MerkleMultiproofReplay.hpp`

拟做内容：
- 为每棵树的 multiproof 构造一次 `index -> position` 辅助结构。
- query 验证阶段直接 O(1) 或近似 O(1) 取值，不再重复二分。

验收：
- `tests/test_pcs.cpp`
- extension-challenge 相关 tests
- 高 query 数下复测 verifier：

```bash
./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --queries 2048 --warmup 1 --reps 3
```

已完成内容：
- 在 `MerkleMultiproofReplay` 共享层新增 verifier 侧 `ValuePositionCache`。
- 为每棵树的 multiproof 一次性构造 `index -> position` 查找表，替代 query 内层循环中的重复 `lower_bound`。
- base verifier 与 extension verifier 都接到这套 shared helper：
  - `src/PCS/BaseFold/BaseFoldPCSVerify.cpp`
  - `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`
- proof layout、serializer、Merkle multiproof 校验语义均未改变。
- 新增 shared helper 单测，覆盖命中、miss、长度不匹配、重复索引拒绝。

本轮验证：
- `cmake --build build-release -j 4 --target test_pcs bench_basefold_pcs_eval`
- `ctest --test-dir build-release --output-on-failure -R '^test_pcs$'`
- 对照基准：
  - A2 baseline: commit `03a808f`
  - A3 current workspace
  - 统一参数：`BASEFOLD_VERIFY_QUERY_MAX_THREADS=1`, `OMP_NUM_THREADS=1`
  - 命令：

```bash
./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --queries 2048 --profile --warmup 1 --reps 3
```

本轮测得收益：
- verifier mean：`42.974 ms -> 40.632 ms`，下降 `2.342 ms`，约 `5.4%`
- `VerifyQueryMultiproofs` 总时间：`115.782 ms -> 108.707 ms`，下降 `7.075 ms`，约 `6.1%`
- `Inside queries other`：`33.146 ms -> 28.014 ms`，下降 `5.132 ms`

备注：
- 这次收益集中在 verifier query 路径，符合 A3 目标。
- `MerkleVerifyMultiproof` 本身基本不变，说明提升主要来自值回放查找，而不是底层 Merkle replay。
- extension verifier 也已切到相同 cache helper；这轮只对 base bench 做了高-query 定量，extension 路径由 `test_pcs` 覆盖 correctness。

## 优先级 B

### B1. BaseFold sumcheck 初始化缓存到 committed witness

状态：`completed`

前提：
- 只在“同一 commitment 反复 opening / prove”是现实工作流时做。
- 如果后续使用场景主要是一条 witness 只 prove 一次，则这项降级。

目的：
- 避免每次 prove 都重复做 monomial-to-boolean transform。

主要范围：
- `src/PCS/Common/Sumcheck.cpp`
- `src/PCS/BaseFold/BaseFoldPCSProve.cpp`
- `include/PCS/BaseFold/BaseFoldPCS.hpp`

拟做内容：
- 研究把 `SumcheckProver` 初始化所需的预处理工件缓存到 prover-local committed witness。
- 只缓存 prover 本地中间态，不进入 proof，不改变 verifier。

本轮已完成：
- 在 `Sumcheck.cpp` 中拆出 `SumcheckMonomialPrecomputation`，只缓存 witness 相关的 base boolean-eval table。
- 将该 precompute 挂到 `BaseFoldPCSCommitArtifacts`，由 `BaseFoldPCSBuildCommitArtifacts*` 构建并由 base prove 路径复用。
- 保持 `z` 相关的 prefix equality table 仍在 prove 时按 opening 点现算，不把 point-dependent 状态混入 commit artifacts。
- 增加 regression test：同一组 `commit_artifacts` 在开启 / 关闭该 cache 时生成的 proof 固定序列化字节一致。

风险点：
- 内存占用会上升。
- 需要严格限定为 prover-local cache，不能把“只是为了加速”的状态混进协议对象。

验收：
- BaseFold tests 全绿。
- 在“同一 commit，多次 prove”微基准下看到收益。

代表性结果：
- `OMP_NUM_THREADS=1`, `GR(4,2)`, `d=12`, `queries=4`, `warmup=1`, `reps=3`
- `bench_basefold_pcs_prove`:
  - prove-phase `18.426 ms -> 17.813 ms`
  - `SumcheckProver init` `6.961 ms -> 1.762 ms`（总计 3 reps）
- `bench_basefold_pcs_eval`:
  - prove-phase `16.924 ms -> 15.538 ms`
  - `SumcheckProver init` `7.103 ms -> 2.777 ms`（总计 3 reps）
- 结论：这轮确实把 base path 里“重复 monomial-to-boolean transform”的大头挪出了 prove 热路径；当前低-query prove-phase 的下一个主要桶已转到 folding commit rounds。

extension-challenge 后续方向：
- 当前 `bench_basefold_pcs_prove --use-extension-challenges` 的主要桶是 `ExtensionCommitRound`，其次是 `ExtensionSumcheck` 的 init / receive。
- 下一步如果要优化 extension-challenge，优先顺序应是：
  - 复用这轮缓存的 base boolean-eval table，并在 prove 时 lift 到 `E(U)`，避免重新跑 `BooleanEvalTableFromMonomialCoeffsExtension`
  - 再看 extension 版 commit round 的 scratch / denom 复用
  - 最后再考虑 extension prefix-equality 的 point-dependent预处理

### B2. folding commit round 的 scratch buffer 复用

状态：`pending`

目的：
- 减少 `denoms` / `inv_denoms` / prefix 中间量的重复分配。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSCommon.cpp`
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

拟做内容：
- 用 thread-local 或调用方传入 scratch 的方式复用临时向量。
- 不引入跨 NTL context 的不安全缓存对象。

验收：
- correctness 不变。
- 对 BaseFold prove-phase 有小幅稳定收益即可，不要求大幅提升。

### B3. 继续针对编码器做“常见参数族”专门化

状态：`completed`

目的：
- 既然 `bench_basefold_pcs_commit` 里 encode 仍明显大于 top Merkle，就继续只做编码器里高信号的专门化。

主要范围：
- `src/PCS/BaseFold/FoldableCode.cpp`
- `bench/bench_basefold_pcs_commit.cpp`

拟做内容：
- 保持当前 iterative encoder 主体不变。
- 只补常见参数分支上的专门化，例如：
  - `diag_T` 常量或全 1 的更强 fast path；
  - 常见 `k0` / `c` 组合下减少 NTL 临时对象；
  - 避免无收益的数据搬运。

本轮已完成：
- 在 `FoldableCode.cpp` 的 `k0==1` level-0 路径上增加 `G0 = [1, zeta, 1, ...]` 探测。
- 对 `c==1`、`c==2` 和一般 `c` 分别走更轻量的赋值/单次乘法路径，避免每列都做通用 `m * G0[0][j]`。
- 在 `tests/test_foldable_codes.cpp` 增加 `[1, zeta, 1, ...]` 行模式的 checked / unchecked 定向回归。

不做：
- 不引入大规模模板膨胀。
- 不做平台相关 SIMD 手写。

验收：
- `bench_basefold_pcs_commit` 的 encode-only mean 下降。
- `tests/test_foldable_codes.cpp` 通过。

代表性结果：
- `OMP_NUM_THREADS=1`, `GR(4,2)`, `d=12`, `warmup=1`, `reps=5`
- `c=2, k0=1`: encode-only `29.507 ms -> 28.909 ms`，commit `32.713 ms -> 32.111 ms`
- `c=4, k0=1`: encode-only `57.094 ms -> 56.857 ms`，commit `63.488 ms -> 63.406 ms`
- 结论：这类 level-0 常见参数专门化有小幅稳定收益，但整体 commit 仍主要由更深层 folding 编码主导。

## Extension-Challenge 专项计划

状态：`pending`

目标：
- 单独针对 `BaseFoldPCSChallengeConfig{use_extension_challenges=true}` 模式做一条可跟踪的优化线。
- 只做 prover-local 或 benchmark-hot-path 优化，不改变 proof 语义、transcript 顺序或 verifier 约束。

当前基线：
- 基准：`OMP_NUM_THREADS=1 BASEFOLD_VERIFY_QUERY_MAX_THREADS=1 ./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --queries 4 --profile --warmup 1 --reps 3 --use-extension-challenges`
- 当前 `prove-phase mean = 65.287 ms`
- 主要桶：
  - `ExtensionCommitRound = 82.923 ms / 3 reps`
  - `ExtensionSumcheck total = 68.239 ms / 3 reps`
    - `init = 14.952 ms`
    - `CurrentPolynomial = 36.578 ms`
    - `ReceiveChallenge = 16.709 ms`
  - `ExtensionMerkleTree::Build = 18.775 ms / 3 reps`
  - `Other = 25.093 ms / 3 reps`
- 当前 `verifier mean = 0.983 ms`，不是主战场。

### EC1. Extension commit round 优化

状态：`completed`

目的：
- 先打当前最大的单桶 `ExtensionCommitRound`。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

拟做内容：
- 复用 `denoms` / `inv_denoms` scratch buffer。
- 补上和 base path 对齐的 `diag_T` all-equal fast path。
- 评估是否把 `inv((zeta-1) * t_j)` 作为 level-local 参数预处理缓存下来；它只依赖 `params.zeta` 和 `diag_T[level]`，不依赖 witness 或 challenge。

不做：
- 不改 extension proof 格式。
- 不把 challenge-dependent 数据放入 commit artifacts。

验收：
- `ExtensionCommitRound` 总时间下降。
- `bench_basefold_pcs_prove --use-extension-challenges --profile` 的 prove-phase 有稳定下降。

结果：
- 已实现按 level 预计算 `inv((zeta-1) * t_j)` 并复用到 extension folding commit round。
- 已补上和 base path 对齐的 `diag_T` all-equal fast path。
- 当前 release 对照（`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=3, ext_deg=2`）：
  - `bench_basefold_pcs_prove --use-extension-challenges`：
    `prove-phase mean 54.687 ms -> 43.802 ms`
    `ExtensionCommitRound 82.695 ms / 3 reps -> 52.114 ms / 3 reps`
  - `bench_basefold_pcs_eval --use-extension-challenges`：
    `prove-phase mean 65.287 ms -> 54.035 ms`
    `ExtensionCommitRound 82.923 ms / 3 reps -> 52.018 ms / 3 reps`

### EC2. 复用 base boolean-eval table 到 extension sumcheck init

状态：`completed`

目的：
- 避免 extension 路径重复跑 `BooleanEvalTableFromMonomialCoeffsExtension(...)`。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`
- `include/PCS/BaseFold/BaseFoldPCS.hpp`

拟做内容：
- 复用 B1 已经缓存到 `BaseFoldPCSCommitArtifacts` 的 base boolean-eval table。
- prove 时把该表 lift 到 `E(U)`，替代“先 lift monomial coeffs，再做 extension zeta transform”的现状。

验收：
- `ExtensionSumcheck init` 明显下降。
- 相关 BaseFold tests 全绿。

结果：
- 已在 extension prove 路径优先复用 `BaseFoldPCSCommitArtifacts::base_sumcheck_precomputation`。
- 当前实现改为“直接 lift 已缓存的 base boolean-eval table”，替代“先 lift monomial coeffs，再跑 `BooleanEvalTableFromMonomialCoeffsExtension(...)`”。
- 当前 release 对照（基于 EC1 之后的同一组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=3, ext_deg=2`）：
  - `bench_basefold_pcs_prove --use-extension-challenges`：
    `prove-phase mean 43.802 ms -> 39.668 ms`
    `ExtensionSumcheck total 37.984 ms / 3 reps -> 26.931 ms / 3 reps`
    `ExtensionSumcheck init 13.888 ms / 3 reps -> 3.152 ms / 3 reps`
  - `bench_basefold_pcs_eval --use-extension-challenges`：
    `prove-phase mean 54.035 ms -> 51.919 ms`
    `ExtensionSumcheck total 67.153 ms / 3 reps -> 58.957 ms / 3 reps`
    `ExtensionSumcheck init 14.808 ms / 3 reps -> 4.404 ms / 3 reps`

### EC3. extension prefix-equality 改为 base-ring 原生表示

状态：`completed`

目的：
- 降低 `ExtensionSumcheck::CurrentPolynomial` 的扩域运算密度。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

拟做内容：
- 把 `z` 和 `prefix_eq_by_vars_` 保持在 base ring，而不是整张表都存 `ZZ_pEX`。
- 保留真正需要在扩域里的 `suffix_eq_prod_` / `f_eval_table_`。
- 去掉 `CurrentPolynomial()` 中反复的 `ExtractBaseConstantCoefficient(...)` 风格绕行。

验收：
- `ExtensionSumcheck current` 下降。
- 不改变 verifier 结果和 proof payload。

结果：
- 已把 `ExtensionSumcheckProver` 内部的 `z` 和 `prefix_eq_by_vars_` 改成 base-ring 表示。
- 已去掉 `CurrentPolynomial()` 中对 `z` / `prefix_eq_by_vars_` 的 `ExtractBaseConstantCoefficient(...)` 绕行；`suffix_eq_prod_` / `f_eval_table_` 仍保留在扩域。
- 当前 release 对照（基于 EC2 之后的同一组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=3, ext_deg=2`）：
  - `bench_basefold_pcs_prove --use-extension-challenges`：
    `prove-phase mean 39.668 ms -> 39.575 ms`，小幅下降
    `ExtensionSumcheck current 8.112 ms / 36 calls -> 7.998 ms / 36 calls`
  - `bench_basefold_pcs_eval --use-extension-challenges`：
    `prove-phase mean 51.919 ms -> 51.165 ms`
    `ExtensionSumcheck total 58.957 ms / 3 reps -> 57.225 ms / 3 reps`
    `ExtensionSumcheck current 37.461 ms / 36 calls -> 35.989 ms / 36 calls`
- 结论：这一步是温和收益，不是 EC1/EC2 那种大幅下降；主要改善体现在 `eval` 路径的 `CurrentPolynomial` 桶，`prove-only` 路径只有小幅改善。

### EC4. unchecked 热路径移出 prover-local 自检

状态：`completed`

目的：
- 去掉 extension unchecked 路径里不属于协议消息本身的本地重算。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

拟做内容：
- 把 `expected_pi0` / `EncodeC0Extension(...)` 这类一致性自检留在 checked 路径。
- unchecked 路径只保留真正生成 proof 需要的计算。

验收：
- unchecked prove/eval 的结果不变。
- `Other` 或 headline `prove-phase mean` 有可见下降。

结果：
- 已把 `expected_pi0` / `EncodeC0Extension(...)` 的一致性自检移出 unchecked 路径，只保留在 checked 包装层。
- 已补 checked / unchecked 一致性回归，固定宽度 proof 序列化保持一致。
- 当前 release 对照（基于 EC3 之后的同一组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=3, ext_deg=2`）：
  - `bench_basefold_pcs_prove --use-extension-challenges`：
    `prove-phase mean 39.575 ms -> 40.041 ms`
  - `bench_basefold_pcs_eval --use-extension-challenges`：
    `prove-phase mean 51.165 ms -> 51.850 ms`
- 结论：这一步主要是把 unchecked/checked 的语义边界理顺，当前 benchmark 尺度下的 end-to-end 收益太小，3 reps 的前后对比基本落在噪声带内，不适合宣称为明显加速。

### EC5. lifted top oracle 与 suffix 工件缓存

状态：`completed`

目的：
- 面向“同一 commitment 多次 extension opening / prove”的工作流，减少重复 lifting 和 suffix 收尾开销。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`
- `include/PCS/BaseFold/BaseFoldPCS.hpp`

本轮完成：
- 在 `BaseFoldPCSCommitArtifacts` 里缓存 lifted top oracle，供 extension prove repeated-prove 热路径直接复用。
- unchecked extension prove 在 cache 缺失时仍保留本地 fallback，不改 proof / verifier 语义。
- 明确不缓存 `Msg0CoeffsAtSuffixChallenges(...)`；它依赖本次 transcript 采样出的 `r_by_level`，不能安全放进通用 commit artifacts。

验收：
- repeated-prove 微基准下降。主验收面使用 `bench_basefold_pcs_prove`，因为它会在 reps 之间复用同一个 `commit_artifacts`；`bench_basefold_pcs_eval` 每轮都会重建 artifacts，不适合作为 EC5 主指标。
- 不把 point-dependent结果错误缓存到 commit artifacts。

当前 release 对照（基于 `57376dd` 的同组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=10, ext_deg=2`）：
- `bench_basefold_pcs_prove --use-extension-challenges`：
  `prove-phase mean 41.822 ms -> 37.359 ms`
- 结论：这一步对 repeated extension prove 是有意义的，下降约 `10.7%`。收益主要来自把 `LiftOracleToExtension(commit_artifacts.pi_d)` 从计时内移到 artifact 预处理，并避免在 prove 热路径里重复构造/复制 top extension oracle。

不建议现在优先做：
- extension verifier 优化；当前 verifier 还不是问题。
- extension Merkle 小修小补；`ExtensionMerkleTree::Build` 不是第一大头。
- transcript / serializer 微优化。

### extension-challenge 后续执行计划（EC6+）

当前热点基线（`build-release`, `OMP_NUM_THREADS=1`, `GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=3, ext_deg=2`）：
- `bench_basefold_pcs_eval --use-extension-challenges`
  - `prove-phase mean 51.759 ms`
  - `ExtensionSumcheck total 57.932 ms / 3 reps`
  - `CurrentPolynomial 37.179 ms / 36 calls`
  - `ReceiveChallenge 17.352 ms / 36 calls`
  - `ExtensionCommitRound 53.482 ms / 3 reps`
- `bench_basefold_pcs_prove --use-extension-challenges`
  - `prove-phase mean 37.774 ms`
  - `ExtensionCommitRound 51.274 ms / 3 reps`
  - `ExtensionSumcheck total 26.819 ms / 3 reps`
  - `ExtensionMerkleTree::Build 17.282 ms / 3 reps`

结论：
- 如果目标是完整 extension `eval` prove-phase，下一优先项应先打 `ExtensionSumcheck`。
- 如果目标是 repeated-prove，`ExtensionCommitRound` 的扩域算术内核也值得继续打。
- `B2` 原始定义里的 scratch-buffer 复用已经不是 extension 主热点；当前快路径大多直接走 precomputed inverse，不再频繁落到 `denoms` / `inv_denoms` 临时分配。

### EC6. `ExtensionSumcheck::CurrentPolynomial()` 因式重排

状态：`completed`

目的：
- 降低 `CurrentPolynomial()` 内层循环里的通用扩域乘法数量。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

本轮完成：
- 已把 `CurrentPolynomial()` 改成先累积
  `S0 = Σ prefix * f0` 与 `S1 = Σ prefix * (f1-f0)`，
  再在循环外组合成
  `a0 = suffix * factor0 * S0`
  `a1 = suffix * (delta * S0 + factor0 * S1)`
  `a2 = suffix * delta * S1`
- 内层循环现在只保留“base 常数乘扩域元素 + 加法”；通用 `MulExtension(...)` 从每个 `mask` 的 4 次压缩到循环外 2 次。

验收：
- `bench_basefold_pcs_eval --use-extension-challenges` 的 `ExtensionSumcheck current` 明显下降。
- proof / verifier 结果不变。

当前 release 对照（基于 EC5 之后的同组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=5, ext_deg=2`）：
- `bench_basefold_pcs_eval --use-extension-challenges`
  - `prove-phase mean 51.367 ms -> 39.519 ms`
  - `ExtensionSumcheck total 95.726 ms / 5 reps -> 50.627 ms / 5 reps`
  - `CurrentPolynomial 61.128 ms / 60 calls -> 17.282 ms / 60 calls`
- `bench_basefold_pcs_prove --use-extension-challenges`
  - `prove-phase mean 37.691 ms -> 37.931 ms`
  - `CurrentPolynomial 13.384 ms / 60 calls -> 11.641 ms / 60 calls`
- 结论：EC6 对完整 extension `eval` 路径是大收益，直接压掉了 `CurrentPolynomial()` 这个最大单桶；但 repeated-prove headline 基本持平，因为剩余主瓶颈已经变成 `ExtensionCommitRound` 与 `ReceiveChallenge()`。

### EC7. `ExtensionSumcheck::ReceiveChallenge()` in-place 与 modulus hoist

状态：`completed`

目的：
- 降低 `ReceiveChallenge()` 的扩域临时对象 churn 和 helper 调度开销。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

本轮完成：
- 在 `ReceiveChallenge()` 里引入 reduced-input 专用乘法 helper，避免每个 fold 都重新走通用 `AddExtension` / `SubExtension` / `MulExtension` 包装。
- 把 `mod_ctx` hoist 到 round 局部，并把 `EqFactorExtensionFromBase` 改成不依赖通用 helper 的 reduced 版本。
- 复用少量局部 scratch `ZZ_pEX`，减少 `delta_f` / folded 值的临时对象 churn。

验收：
- `bench_basefold_pcs_eval --use-extension-challenges` 的 `ExtensionSumcheck receive` 稳定下降。
- correctness 不变，`test_pcs` 继续通过。

当前 release 对照（基于 EC6 之后的同组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=5, ext_deg=2`）：
- `bench_basefold_pcs_eval --use-extension-challenges`
  - `prove-phase mean 39.519 ms -> 39.504 ms`
  - `ExtensionSumcheck total 50.627 ms / 5 reps -> 48.093 ms / 5 reps`
  - `ReceiveChallenge 27.771 ms / 60 calls -> 24.490 ms / 60 calls`
- `bench_basefold_pcs_prove --use-extension-challenges`
  - `prove-phase mean 37.931 ms -> 36.668 ms`
  - `ReceiveChallenge 26.906 ms / 60 calls -> 22.820 ms / 60 calls`
- 结论：EC7 的收益主要体现在 `ReceiveChallenge()` 本身，eval headline 基本持平但 repeated-prove 已有小幅下降；当前 extension 路径剩余最大桶已经更明确地落在 `ExtensionCommitRound`，所以下一步优先做 EC8。

### EC8. extension commit round 的 fused arithmetic kernel

状态：`completed`

目的：
- 继续优化 `EvalLineAtExtensionWithInvDenom(...)` 一带的扩域算术常数，而不是继续追 scratch-buffer 分配。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

本轮完成：
- 已在 extension commit round 内部引入 reduced-input 的 fused helper，直接把
  `delta_y`、乘 base 常数、`delta_x` 和最终加回 `y1` 串成一条 in-place 路径。
- 已把 `mod_ctx` hoist 到 `ProverCommitRoundExtensionNoValidate()` 的 round 局部，让四条热循环共享同一份 modulus context。
- 本轮没有继续做 `lambda_j` 预计算；先吃掉低风险的临时对象和 helper 调度开销。

验收：
- `bench_basefold_pcs_prove --use-extension-challenges` 的 `ExtensionCommitRound` 稳定下降。
- `bench_basefold_pcs_eval --use-extension-challenges` headline `prove-phase mean` 至少有可见改善。

当前 release 对照（基于 EC7 之后的同组参数：`GR(4,2), c=2, k0=1, d=12, queries=4, warmup=1, reps=5, ext_deg=2`）：
- `bench_basefold_pcs_eval --use-extension-challenges`
  - `prove-phase mean 39.633 ms -> 36.365 ms`
  - `ExtensionCommitRound 88.465 ms / 60 calls -> 73.199 ms / 60 calls`
- `bench_basefold_pcs_prove --use-extension-challenges`
  - `prove-phase mean 37.177 ms -> 34.182 ms`
  - `ExtensionCommitRound 86.955 ms / 60 calls -> 72.676 ms / 60 calls`
- 结论：EC8 是 extension 线里又一波大收益，直接把当前最大剩余桶压下去了。到这一步为止，`EC6 + EC7 + EC8` 已经把 extension sumcheck 与 commit round 的主要低风险优化都吃掉了；是否继续做 `EC9`，应取决于你是否真的长期固定在 `ext_deg=2` 工作负载上。

### EC9. `ext_deg=2` 热路径专门化

状态：`pending`

目的：
- 当主要 workload 固定在 `ext_deg=2` 时，进一步降低通用 `ZZ_pEX` helper 的常数开销。

主要范围：
- `src/PCS/BaseFold/BaseFoldPCSExtension.cpp`

拟做内容：
- 评估在 prover 热路径内部引入 degree-2 专门化表示，只在边界保留 `ZZ_pEX`。
- 只在前面的 EC6~EC8 做完后再决定是否值得推进，避免过早把原型代码复杂化。

验收：
- `bench_basefold_pcs_eval --use-extension-challenges` 与 `bench_basefold_pcs_prove --use-extension-challenges` 均有稳定收益。
- 不破坏现有通用扩域路径。

## 执行顺序

1. `A1` Frobenius / RingSwitch unchecked 热路径
2. `A2` equality table 共享迭代 helper
3. `A3` verifier multiproof 查找表
4. 复测 release benches
5. 若 outer prover 仍是主要瓶颈，做 `B1`
6. 若 commit 仍明显重于 Merkle，做 `B3`
7. `B2` 作为补充整理项穿插进行
8. 若开始专门优化 extension-challenge，先完成 `EC1 -> EC2 -> EC3 -> EC4 -> EC5`
9. extension 后续优化按 `EC6 -> EC7 -> EC8` 顺序推进
10. `B2` 原始 scope 暂时后置；只有在 fallback/no-cache 路径重新变热时再做
11. `EC9` 仅在 `ext_deg=2` 仍是主 workload 且 EC6~EC8 收益已吃完时再评估

## 每轮变更后的固定检查

```bash
cmake --build build-release
ctest --test-dir build-release --output-on-failure
./build-release/bench_basefold_pcs_eval --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --queries 4 --profile --warmup 0 --reps 1
./build-release/bench_basefold_pcs_commit --mode ring --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 12 --warmup 0 --reps 1
./build-release/bench_z2k_frobenius_eval --c 4 --ell 10 --kappa 1 --queries 4 --warmup 0 --reps 1 --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1
```

## 跟踪表

| ID | 任务 | 优先级 | 状态 | 主要文件 | 验收基准 |
|---|---|---|---|---|---|
| A1 | Frobenius / RingSwitch unchecked 热路径 | A | completed | `src/Compiler/Z2k/*PCS.cpp` | outer prover mean |
| A2 | equality table 迭代 helper | A | completed | `src/Compiler/Z2k/*PCS.cpp`, `src/PCS/Common/*` | outer prover mean |
| A3 | multiproof 查找表 | A | completed | `BaseFoldPCSVerify.cpp`, `BaseFoldPCSExtension.cpp` | verifier mean |
| B1 | BaseFold sumcheck init cache | B | completed | `Sumcheck.cpp`, `BaseFoldPCSProve.cpp` | repeated-prove 微基准 |
| B2 | scratch buffer 复用 | B | pending | `BaseFoldPCSCommon.cpp`, `BaseFoldPCSExtension.cpp` | prove-phase 小幅下降 |
| B3 | 编码器常见参数专门化 | B | completed | `FoldableCode.cpp` | encode-only mean |
| EC1 | Extension commit round 优化 | Ext | completed | `BaseFoldPCSExtension.cpp`, `BaseFoldPCSCommit.cpp`, `BaseFoldPCS.hpp` | `ExtensionCommitRound` |
| EC2 | 复用 base eval table 到 extension sumcheck init | Ext | completed | `BaseFoldPCSExtension.cpp` | `ExtensionSumcheck init` |
| EC3 | extension prefix-equality base-ring 化 | Ext | completed | `BaseFoldPCSExtension.cpp` | `ExtensionSumcheck current` |
| EC4 | extension unchecked 自检移出热路径 | Ext | completed | `BaseFoldPCSExtension.cpp` | ext prove-phase mean |
| EC5 | lifted top oracle / suffix 工件缓存 | Ext | completed | `BaseFoldPCSExtension.cpp`, `BaseFoldPCS.hpp` | repeated-prove 微基准 |
| EC6 | `CurrentPolynomial()` 因式重排 | Ext | completed | `BaseFoldPCSExtension.cpp` | `ExtensionSumcheck current` |
| EC7 | `ReceiveChallenge()` in-place / modulus hoist | Ext | completed | `BaseFoldPCSExtension.cpp` | `ExtensionSumcheck receive` |
| EC8 | commit round fused arithmetic kernel | Ext | completed | `BaseFoldPCSExtension.cpp` | `ExtensionCommitRound` |
| EC9 | `ext_deg=2` 热路径专门化 | Ext | pending | `BaseFoldPCSExtension.cpp` | ext prove-phase mean |
