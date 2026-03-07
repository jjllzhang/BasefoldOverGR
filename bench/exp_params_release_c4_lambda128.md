# Release 实验参数（`c=4`, `lambda=128`, 维度 sweep）

本文档对应脚本 `scripts/run_release_c4_lambda128.sh`。

## 1) 目标场景

- 码率：`1/4`，即 `c=4`
- 安全性：`lambda=128`
- 多项式维度：`2^3` 到 `2^29`
  - 在 bench 参数中对应 `d=3..29`（`k0=1` 时，`poly_dim = k_d = 2^d`）
- 构建：`Release`
- 指标：
  - `commit time`：`bench_pcs_commit` 的 `encode-only mean`
  - `prover time`：`bench_pcs_eval` 的 `prover mean`
  - `verifier time`：`bench_pcs_eval` 的 `verifier mean`
  - `proof size`：`bench_pcs_eval` 对真实 proof 走 fixed-width `CountingSink`
    得到的精确序列化大小（`proof_size_bytes` / `proof_size_kb`）

## 2) 实验上下文（默认 `all` 共 12 组）

- `Field-255`
  - `--mode field`
  - `--field-mod = 2^255 - 19`
  - `--field-F = 1,1`
  - `--field-zeta = 0,1`
  - `calc_iopp_params` 使用 `q = p^(r*m), p=2^255-19, r=1, m=1`

- `GR(2^16,162)`
  - `--mode ring`
  - `--ring-mod = 2^16 = 65536`
  - `--ring-p = 2`
  - `--ring-F = RING_F_162`（脚本内置的 162 次模 2 不可约多项式）
  - `--ring-zeta = 0,1`
  - `calc_iopp_params` 使用 `q = 2^162`

- `F_2^256`
  - `--mode field`
  - `--field-mod = 2`
  - `--field-F = F2_256_F`（脚本内置的 256 次模 2 不可约多项式）
  - `--field-zeta = 0,1`
  - `calc_iopp_params` 使用 `q = 2^256`

- `GR(2^2,162)`
  - `--mode ring`
  - `--ring-mod = 2^2 = 4`
  - `--ring-p = 2`
  - `--ring-F = RING_F_162`
  - `--ring-zeta = 0,1`
  - `calc_iopp_params` 使用 `q = 2^162`

- `64-bit prime`（`field-prime64-ext`）
  - `--mode field`
  - `--field-mod = 18446744073709551557`（`2^64 - 59`）
  - `--field-F = 1,1`
  - `--field-zeta = 0,1`
  - 使用三次扩域挑战（`--use-extension-challenges` + `E(U)=1+U+U^3`）
  - `calc_iopp_params` 使用 `q = p^(r*m) = p^(1*3) = p^3`

- `F_2^64`（`field-f2p64-ext`）
  - `--mode field`
  - `--field-mod = 2`
  - `--field-F = RING_F_64`（脚本内置 64 次模 2 不可约多项式）
  - `--field-zeta = 0,1`
  - 使用三次扩域挑战（`--use-extension-challenges` + 默认 `E(U)=1+U+U^3`）
  - `calc_iopp_params` 使用 `q = 2^(64*3) = 2^192`

- `128-bit prime`（`field-prime128-ext`）
  - `--mode field`
  - `--field-mod = 326594724262804054738278293730872375507`
  - `--field-F = 1,1`
  - `--field-zeta = 0,1`
  - 使用扩域挑战（`--use-extension-challenges` + 二次扩展挑战多项式）
  - `calc_iopp_params` 使用 `q = p^(r*m) = p^(1*2) = p^2`

- `F_2^128`（`field-f2p128-ext`）
  - `--mode field`
  - `--field-mod = 2`
  - `--field-F = x^128 + x^7 + x^2 + x + 1`（AES-GCM 常用不可约多项式）
  - `--field-zeta = 0,1`
  - 使用扩域挑战（`--use-extension-challenges` + 二次 Artin-Schreier 多项式 `U^2 + U + beta`，`Tr(beta)=1`）
  - `calc_iopp_params` 使用 `q = 2^(128*2) = 2^256`

- `GR(2^16,64)`（`ring-gr-2p16-64-ext`）
- `GR(2^16,128)`（`ring-gr-2p16-128-ext`）
- `GR(2^2,64)`（`ring-gr-2p2-64-ext`）
- `GR(2^2,128)`（`ring-gr-2p2-128-ext`）
  - `GR(*,128)` 的 `--ring-F` 使用上述 128 次多项式在 `Z/(2^k)` 上的 Hensel 提升（脚本中即同系数 0/1 表示）
  - `GR(*,64)` 使用三次扩环挑战（默认 `E(U)=1+U+U^3`）
  - `GR(*,128)` 使用二次 Artin-Schreier 扩环挑战（`U^2 + U + gamma`，`gamma mod 2 = beta` 且 `Tr(beta)=1`）
  - `calc_iopp_params` 的 `q`（按 `q=p^(r*m)`，且 ring 侧取 `p=2`）：
    - `ring-gr-2p16-64-ext` / `ring-gr-2p2-64-ext`：`q = 2^(64*3) = 2^192`
    - `ring-gr-2p16-128-ext` / `ring-gr-2p2-128-ext`：`q = 2^(128*2) = 2^256`

## 3) queries（每个 d、每个上下文单独推导）

脚本对每个 `(context, d)` 都会调用：

```bash
./build-release/calc_iopp_params \
  --d <d> --c 4 --k0 1 --lambda 128 \
  --p <context-specific-p> --r <context-specific-r> --m <context-specific-m> \
  --auto-gamma
```

并从输出解析：

- `gamma`
- `l_min_for_PCS`（作为该点的 `queries`）

说明：

- `calc_iopp_params` 统一使用 `q = p^(r*m)`。
- 对 ring context，脚本传入的是 `p=2`（残差域视角），因此这里的 `q` 不是环总基数 `(2^k)^(r*m)`，而是 `2^(r*m)`。

## 4) 运行方式

- 直接运行：

```bash
scripts/run_release_c4_lambda128.sh
```

- 常用环境变量：
  - `RUN_ID`：本次运行 ID（默认 `<timestamp>_pid<shell-pid>`），用于区分输出目录和（可选）构建目录
  - `D_MIN` / `D_MAX`：维度区间（默认 `3..29`）
  - `CONTEXTS`：选择上下文，默认 `all`
    - 可选值：`field-255,ring-gr-2p16-162,field-f2p256,ring-gr-2p2-162,field-prime64-ext,field-f2p64-ext,field-prime128-ext,field-f2p128-ext,ring-gr-2p16-64-ext,ring-gr-2p16-128-ext,ring-gr-2p2-64-ext,ring-gr-2p2-128-ext`
    - 示例：`CONTEXTS=field-prime128-ext` 或 `CONTEXTS=field-f2p128-ext,ring-gr-2p16-64-ext`
    - 兼容别名：`field-prime64 -> field-prime64-ext`，`field-f2p64 -> field-f2p64-ext`，`field-prime128 -> field-prime128-ext`，`field-f2p128 -> field-f2p128-ext`
  - `BENCH_THREADS`：单个 bench 进程内部线程数（默认 `8`）
    - 默认情况下脚本会设置 `OMP_NUM_THREADS=8`，并把
      `BASEFOLD_MERKLE_MAX_THREADS`、`BASEFOLD_VERIFY_QUERY_MAX_THREADS`
      设为同值（若你未手动设置）。
    - 若希望完全使用运行时默认线程策略，可显式设置 `BENCH_THREADS=0`。
  - `CPU_PIN_MODE`：`none` / `manual` / `slot`（默认 `none`）
    - `manual`：用 `CPU_SET` 手动指定 CPU 集合（如 `0-31`）
    - `slot`：按 `RUN_SLOT` 和 `RUN_SLOTS_TOTAL` 自动切分 CPU，适合并发多实例
  - `CPU_SET`：当 `CPU_PIN_MODE=manual` 时生效
  - `RUN_SLOT` / `RUN_SLOTS_TOTAL`：当 `CPU_PIN_MODE=slot` 时生效（`RUN_SLOT` 从 `0` 开始）
  - `USE_SMT_IN_SLOT`：`0` 或 `1`（默认 `0`）
    - `0`：按物理核切分（每核只取一个硬件线程）
    - `1`：按逻辑核切分（包含 SMT sibling）
  - `PIN_BUILD`：`0` 或 `1`（默认 `0`）
    - 设为 `1` 时，`cmake configure/build` 也会被 `taskset` 绑核
  - `ISOLATE_BUILD_DIR`：`0` 或 `1`（默认 `0`）
    - 设为 `1` 且未手动指定 `BUILD_DIR` 时，会使用 `build-release-<RUN_ID>`，避免并发实例共享同一个 build 目录
  - `BUILD_DIR`：显式指定构建目录（优先级高于 `ISOLATE_BUILD_DIR`）
  - `CMD_TIMEOUT_SEC`：单条 bench 超时秒数（默认 `0`，即不超时）
  - `CONTINUE_ON_ERROR`：遇到某个点失败后是否继续（默认 `1`）
  - `COMMIT_REPS` / `EVAL_REPS` / `SEED`

并发多实例建议（互不影响）：

```bash
# 实例 0
ISOLATE_BUILD_DIR=1 CPU_PIN_MODE=slot RUN_SLOT=0 RUN_SLOTS_TOTAL=2 \
CONTEXTS=ring-gr-2p16-64-ext scripts/run_release_c4_lambda128.sh

# 实例 1
ISOLATE_BUILD_DIR=1 CPU_PIN_MODE=slot RUN_SLOT=1 RUN_SLOTS_TOTAL=2 \
CONTEXTS=ring-gr-2p16-128-ext scripts/run_release_c4_lambda128.sh
```

说明：单个脚本内部仍按 `d` 串行推进；同一 `d` 下各 context 也串行执行。并行度主要来自单个 bench 进程内部线程；多脚本并发时建议使用上面的 CPU 分片与独立 build 目录。

## 5) 输出

脚本输出目录：`results/release_c4_lambda128_sweep_<RUN_ID>/`

- 明细 csv：`results.csv`
  - 每行一个 `(context, d)` 点，含 `gamma/queries/4项指标/status/error`
  - 其中 `proof_size_bytes` / `proof_size_kb` 来自 `bench_pcs_eval` 对真实 proof
    的 fixed-width counting 结果
- 汇总 markdown：`RESULTS.md`
- 原始日志：`logs/*.log`（来自 `calc_iopp_params`、`bench_pcs_commit`、
  `bench_pcs_eval`）
