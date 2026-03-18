# BasefoldOverGR

本仓库包含一组在有限域和 Galois ring 上实现 BaseFold 风格编码、IOPP/PCS，以及 `Z_{2^k}` compiler bench 的 C++ 代码。

当前 bench 工作流以 [scripts/run_release_c4_lambda128.sh](/home/zjl/BasefoldOverGR/scripts/run_release_c4_lambda128.sh) 为中心。

## 基本目录结构

```text
.
├── CMakeLists.txt
├── README.md
├── README_zh.md
├── LICENSE
├── bench/
│   ├── bench_basefold_pcs_*.cpp
│   ├── bench_z2k_ring_switch_*.cpp
│   ├── bench_z2k_frobenius_*.cpp
│   ├── calc_iopp_params.cpp
│   └── exp_params_release_c4_lambda128.md
├── include/
├── src/
├── tests/
├── scripts/
│   ├── run_release_c4_lambda128.sh
│   └── plot_benchmark_results.py
├── results-legacy/
├── build/             # 生成目录
├── build-release/     # 生成目录
└── results/           # 生成的 bench 输出
```

## 如何使用 `scripts/run_release_c4_lambda128.sh`

这个脚本完全通过环境变量配置，没有 CLI 参数。默认会自动配置/构建 `build-release`，并把输出写到 `results/release_c4_lambda128_sweep_<RUN_ID>/`；如果显式设置了 `OUT_DIR`，则写到你指定的位置。

### 默认 BaseFold Release Sweep

```bash
scripts/run_release_c4_lambda128.sh
```

当前默认值：

- `RUN_SUITE=basefold_release`
- `CONTEXTS=all`
- `D_MIN=3`
- `D_MAX=29`
- `C=4`
- `K0=1`
- `LAMBDA=128`
- `COMMIT_WARMUP=1`, `COMMIT_REPS=3`
- `EVAL_WARMUP=1`, `EVAL_REPS=3`
- `BENCH_THREADS=8`

### 只跑一个 BaseFold Context 和较小的 `d` 区间

```bash
CONTEXTS=field-prime128-ext \
D_MIN=10 D_MAX=20 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

### 跑 Ring-Switch Compiler Eval

```bash
RUN_SUITE=compiler_eval_ring_switch \
CONTEXTS=ring-gr-2p16-64-ext \
COMPILER_KAPPA=6 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

### 跑 Frobenius Compiler Eval

```bash
RUN_SUITE=compiler_eval_frobenius \
CONTEXTS=ring-gr-2p16-64-ext \
COMPILER_KAPPA=6 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

### 最常用的环境变量

- `RUN_SUITE`：`basefold_release`、`compiler_eval_ring_switch`、`compiler_eval_frobenius`
- `CONTEXTS`：`all` 或逗号分隔的 context 子集
- `D_MIN`、`D_MAX`：`basefold_release` 下的 BaseFold 维度区间
- `COMPILER_KAPPA`、`COMPILER_ELL_MIN`、`COMPILER_ELL_MAX`：compiler suite 必填
- `K0`：BaseFold 消息基础维度，默认 `1`
- `BENCH_THREADS`：单个 bench 进程内部线程数
- `OUT_DIR`：显式指定输出目录
- `BUILD_DIR`：显式指定构建目录
- `ISOLATE_BUILD_DIR=1`：使用 `build-release-<RUN_ID>`，避免共享 `build-release`
- `CMD_TIMEOUT_SEC`：单条 bench 超时秒数，默认 `0` 表示不超时
- `CONTINUE_ON_ERROR`：单点失败后是否继续，默认 `1`

完整环境变量说明和支持的 context 列表见 [bench/exp_params_release_c4_lambda128.md](/home/zjl/BasefoldOverGR/bench/exp_params_release_c4_lambda128.md)。

## 输出文件

每次运行会输出：

- `backend_eval_results.csv`：来自 `bench_basefold_pcs_commit` 和 `bench_basefold_pcs_eval` 的 BaseFold release 行
- `compiler_eval_results.csv`：所选 family（`ring_switch` 或 `frobenius`）的 compiler-eval 行，包含 `outer_proof_size_*` 和 `total_proof_size_*`
- `RESULTS.md`：markdown 汇总表
- `logs/*.log`：`calc_iopp_params` 和 bench 二进制的原始日志

说明：

- `RUN_SUITE=compiler_eval_frobenius` 时，`GR(2^2,r)` context 不会实际跑 Frobenius full eval，而是直接写 `status=disabled_gr2p2_context`。
