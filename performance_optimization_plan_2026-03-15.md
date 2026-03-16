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

状态：`pending`

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

风险点：
- 内存占用会上升。
- 需要严格限定为 prover-local cache，不能把“只是为了加速”的状态混进协议对象。

验收：
- BaseFold tests 全绿。
- 在“同一 commit，多次 prove”微基准下看到收益。

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

## 执行顺序

1. `A1` Frobenius / RingSwitch unchecked 热路径
2. `A2` equality table 共享迭代 helper
3. `A3` verifier multiproof 查找表
4. 复测 release benches
5. 若 outer prover 仍是主要瓶颈，做 `B1`
6. 若 commit 仍明显重于 Merkle，做 `B3`
7. `B2` 作为补充整理项穿插进行

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
| B1 | BaseFold sumcheck init cache | B | pending | `Sumcheck.cpp`, `BaseFoldPCSProve.cpp` | repeated-prove 微基准 |
| B2 | scratch buffer 复用 | B | pending | `BaseFoldPCSCommon.cpp`, `BaseFoldPCSExtension.cpp` | prove-phase 小幅下降 |
| B3 | 编码器常见参数专门化 | B | completed | `FoldableCode.cpp` | encode-only mean |
