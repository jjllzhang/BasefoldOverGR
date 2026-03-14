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

状态：`pending`

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

### A2. 共享并重写 equality table 构造为迭代 DP

状态：`pending`

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

### A3. verifier multiproof 查找表一次性构建

状态：`pending`

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

状态：`pending`

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

不做：
- 不引入大规模模板膨胀。
- 不做平台相关 SIMD 手写。

验收：
- `bench_basefold_pcs_commit` 的 encode-only mean 下降。
- `tests/test_foldable_codes.cpp` 通过。

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
| A1 | Frobenius / RingSwitch unchecked 热路径 | A | pending | `src/Compiler/Z2k/*PCS.cpp` | outer prover mean |
| A2 | equality table 迭代 helper | A | pending | `src/Compiler/Z2k/*PCS.cpp`, `src/PCS/Common/*` | outer prover mean |
| A3 | multiproof 查找表 | A | pending | `BaseFoldPCSVerify.cpp`, `BaseFoldPCSExtension.cpp` | verifier mean |
| B1 | BaseFold sumcheck init cache | B | pending | `Sumcheck.cpp`, `BaseFoldPCSProve.cpp` | repeated-prove 微基准 |
| B2 | scratch buffer 复用 | B | pending | `BaseFoldPCSCommon.cpp`, `BaseFoldPCSExtension.cpp` | prove-phase 小幅下降 |
| B3 | 编码器常见参数专门化 | B | pending | `FoldableCode.cpp` | encode-only mean |
