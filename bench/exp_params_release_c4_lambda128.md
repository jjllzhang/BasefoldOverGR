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
  - `proof size`：`bench_pcs_proof_size` 的真实 proof size

## 2) 实验上下文（共 4 组）

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

## 3) queries（每个 d、每个上下文单独推导）

脚本对每个 `(context, d)` 都会调用：

```bash
./build-release/calc_iopp_params \
  --d <d> --c 4 --k0 1 --lambda 128 \
  --p <context-specific-p> --r <context-specific-r> --m 1 \
  --auto-gamma
```

并从输出解析：

- `gamma`
- `l_min_for_PCS`（作为该点的 `queries`）

## 4) 运行方式

- 直接运行：

```bash
scripts/run_release_c4_lambda128.sh
```

- 常用环境变量：
  - `D_MIN` / `D_MAX`：维度区间（默认 `3..29`）
  - `RUN_PROOF_SIZE`：`1` 或 `0`（默认 `1`）
  - `CMD_TIMEOUT_SEC`：单条 bench 超时秒数（默认 `0`，即不超时）
  - `CONTINUE_ON_ERROR`：遇到某个点失败后是否继续（默认 `1`）
  - `COMMIT_REPS` / `EVAL_REPS` / `SEED`

## 5) 输出

脚本输出目录：`results/release_c4_lambda128_sweep_<timestamp>/`

- 明细 csv：`results.csv`
  - 每行一个 `(context, d)` 点，含 `gamma/queries/4项指标/status/error`
- 汇总 markdown：`RESULTS.md`
- 原始日志：`logs/*.log`
