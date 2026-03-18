# Release 实验参数（`c=4`, `lambda=128`, 维度 sweep）

本文档对应脚本 `scripts/run_release_c4_lambda128.sh`。

## 1) 目标场景

- 码率：`1/4`，即 `c=4`
- 安全性：`lambda=128`
- BaseFold release 的多项式维度：`2^3` 到 `2^29`
  - 在 bench 参数中对应 `d=3..29`
  - 输出表里的 `poly_dim` 记录 `k_d = k0*2^d`
    （默认 `k0=1` 时，`poly_dim = 2^d`）
- compiler suite 不再直接 sweep `d`
  - `compiler_eval`：用户给 `COMPILER_KAPPA` 和 `COMPILER_ELL_MIN/MAX`
  - 脚本内部统一算 `d = ell-kappa = ell'`
  - compiler 表里的 `poly_dim` 记录源问题大小 `2^ell`
- 构建：`Release`
- 指标：
  - `family=basefold` 的 `commit time`：`bench_basefold_pcs_commit` 的 `commit mean`
  - `family=basefold` 的 `prover/verifier/proof size`：`bench_basefold_pcs_eval`
    的 `prove-phase mean` / `verifier mean` / fixed-width payload counting
  - ring context 的 compiler full eval：`bench_z2k_*_eval`

关于 compiler bench 的参数来源与限制：

- `compiler_eval` 中，`compiler_kappa` 由脚本环境变量 `COMPILER_KAPPA`
  显式给定，`compiler_ell` 在 `COMPILER_ELL_MIN..COMPILER_ELL_MAX` 上 sweep。
- backend 侧维度统一按 `ell-kappa = d = ell'` 构造。
- 脚本会检查当前 ring context 的扩张次数满足 `deg(F) = 2^compiler_kappa`。
- 当前 ring-switch / Frobenius bench helper 仍然只构造 `k0=1` 的 BaseFold
  backend；因此脚本里当 `K0 != 1` 时，会把 compiler 行记为
  `unsupported_k0`。
- 当前脚本会跳过所有 `GR(2^2,r)` context 上的 Frobenius compiler eval，
  并把对应行记为 `disabled_gr2p2_context`；ring-switch 仍照常运行。
- 若 ring context 的 `deg(F)` 不是 2 的幂，或者虽为 2 的幂但与
  `COMPILER_KAPPA` 不匹配，则 compiler 行会记为
  `unsupported_context_degree`。

## 2) 实验上下文（默认 `all` 共 14 组）

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

- `F_3^40`（`field-f3p40-ext`）
  - `--mode field`
  - `--field-mod = 3`
  - `--field-F = F3_40_F`（脚本内置 40 次模 3 不可约多项式）
  - `--field-zeta = 0,1`
  - 使用三次扩域挑战（`--use-extension-challenges` + `E(U)=zeta+U+U^3`）
  - `calc_iopp_params` 使用 `q = 3^(40*3)`

- `F_3^81`（`field-f3p81-ext`）
  - `--mode field`
  - `--field-mod = 3`
  - `--field-F = F3_81_F`（脚本内置 81 次模 3 不可约多项式）
  - `--field-zeta = 0,1`
  - 使用二次扩域挑战（`--use-extension-challenges` + `E(U)=zeta+U+U^2`）
  - `calc_iopp_params` 使用 `q = 3^(81*2)`

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

## 3) queries（每个 d、每个上下文在当前 `K0` 下单独推导）

脚本对每个 `(context, d)` 都会调用：

```bash
./build-release/calc_iopp_params \
  --d <d> --c 4 --k0 <K0> --lambda 128 \
  --p <context-specific-p> --r <context-specific-r> --m <context-specific-m> \
  --auto-gamma
```

并从输出解析：

- `gamma`
- `l_min_for_IOPP`（作为该点的 `queries`）

说明：

- `calc_iopp_params` 统一使用 `q = p^(r*m)`。
- 对 ring context，脚本传入的是 `p=2`（残差域视角），因此这里的 `q` 不是环总基数 `(2^k)^(r*m)`，而是 `2^(r*m)`。

## 4) 运行方式

- 直接运行：

```bash
scripts/run_release_c4_lambda128.sh
```

- 默认 suite：`RUN_SUITE=basefold_release`
  - 运行 `calc_iopp_params`
  - 写入 `family=basefold` 的 BaseFold `commit + eval`
- 可选 suite：
  - `RUN_SUITE=compiler_eval_ring_switch`：只跑 ring-switch compiler full eval
  - `RUN_SUITE=compiler_eval_frobenius`：只跑 Frobenius compiler full eval

- 常用环境变量：
  - `RUN_SUITE`
    - `basefold_release`
    - `compiler_eval_ring_switch`
    - `compiler_eval_frobenius`
    - 兼容旧写法：`RUN_SUITE=compiler_eval`，但此时还需额外给
      `COMPILER_FAMILY=ring_switch|frobenius`
  - `RUN_ID`：本次运行 ID（默认 `<timestamp>_pid<shell-pid>`），用于区分输出目录和（可选）构建目录
  - `D_MIN` / `D_MAX`：只对 `basefold_release` 生效的维度区间（默认 `3..29`）
  - `K0`：基础消息维度 `k0`（默认 `1`，要求为 2 的幂）
    - `backend_eval_results.csv` 中的 `poly_dim = k_d = K0*2^d`
  - `COMPILER_KAPPA`：只对 compiler eval suite 生效；要求为正整数
  - `COMPILER_ELL_MIN` / `COMPILER_ELL_MAX`
    - 只对 compiler eval suite 生效
    - 要求满足 `COMPILER_ELL_MIN >= COMPILER_KAPPA`
    - 脚本内部按 `ell` sweep，并自动计算 `d = ell-kappa`
    - `compiler_eval_results.csv` 里的 `poly_dim = 2^ell`
  - `COMPILER_FAMILY`
    - 仅在兼容旧写法 `RUN_SUITE=compiler_eval` 时使用
    - 可选值：`ring_switch` / `frobenius`
  - `CONTEXTS`：选择上下文，默认 `all`
    - 可选值：`field-255,ring-gr-2p16-162,field-f2p256,ring-gr-2p2-162,field-prime64-ext,field-f2p64-ext,field-prime128-ext,field-f2p128-ext,field-f3p40-ext,field-f3p81-ext,ring-gr-2p16-64-ext,ring-gr-2p16-128-ext,ring-gr-2p2-64-ext,ring-gr-2p2-128-ext`
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
  - `COMMIT_WARMUP` / `COMMIT_REPS`
  - `EVAL_WARMUP` / `EVAL_REPS`
  - `SEED`

并发多实例建议（互不影响）：

```bash
# 实例 0
ISOLATE_BUILD_DIR=1 CPU_PIN_MODE=slot RUN_SLOT=0 RUN_SLOTS_TOTAL=2 \
CONTEXTS=ring-gr-2p16-64-ext scripts/run_release_c4_lambda128.sh

# 实例 1
ISOLATE_BUILD_DIR=1 CPU_PIN_MODE=slot RUN_SLOT=1 RUN_SLOTS_TOTAL=2 \
CONTEXTS=ring-gr-2p16-128-ext scripts/run_release_c4_lambda128.sh
```

若并发的是 compiler eval suite，请同时给出 `COMPILER_KAPPA` 和
`COMPILER_ELL_MIN/MAX`；若用了兼容旧写法 `RUN_SUITE=compiler_eval`，
还要额外给 `COMPILER_FAMILY`。

说明：`basefold_release` 在脚本内部按 `d` 串行推进；compiler eval suite
在脚本内部按 `ell` 串行推进，并统一使用 `d = ell-kappa`。单次运行只会执行一个 family（ring-switch 或 Frobenius）。并行度主要来自单个 bench 进程内部线程；多脚本并发时建议使用上面的 CPU 分片与独立 build 目录。

## 5) 输出

脚本输出目录：`results/release_c4_lambda128_sweep_<RUN_ID>/`

- `backend_eval_results.csv`
  - 统一表；每行一个 `(family, context, d)` 点
  - `family=basefold`：来自 `bench_basefold_pcs_commit + bench_basefold_pcs_eval`
- `compiler_eval_results.csv`
  - 仅对 ring context 生成；每行一个所选 family 的 `(family, context, ell)` 点
  - 包含 full compiler eval 的 `outer/backend/total` commit、prove、verify
    分解，以及 `outer_proof_size` / `proof_size`
  - 对 `family=frobenius`，凡是 `GR(2^2,r)` context 都不会实际执行 bench，
    而是直接写 `status=disabled_gr2p2_context`
- 汇总 markdown：`RESULTS.md`
- 原始日志：`logs/*.log`（来自 `calc_iopp_params`、`bench_basefold_pcs_commit`、
  `bench_basefold_pcs_eval`、`bench_z2k_*_eval`）
