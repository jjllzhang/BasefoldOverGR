# BasefoldOverGR

英文版见 [README.md](README.md)。

BasefoldOverGR 包含一组在有限域和 Galois ring 上实现 BaseFold 风格编码、
IOPP 和 PCS 的 C++ 代码，以及 `Z_{2^k}` compiler 的 ring-switch /
Frobenius benchmark。

当前主要 benchmark 入口是
[scripts/run_release_c4_lambda128.sh](scripts/run_release_c4_lambda128.sh)。
它会自动在 `build-release/` 下配置 Release 构建，并默认把新跑出的结果写到
`results/`；如果显式设置了 `OUT_DIR`，则写到指定目录。

## 构建与测试

依赖：

- CMake 3.16 或更新版本
- 支持 C++17 的编译器
- NTL 和 GMP
- OpenMP（可选，但多线程 bench 会使用）
- `python3` 和 `matplotlib`，用于
  [scripts/plot_benchmark_results.py](scripts/plot_benchmark_results.py)

手动构建与测试：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

`scripts/run_release_c4_lambda128.sh` 会单独维护自己的 `build-release/`
构建目录。`build/`、`build-release/`、`results/`、`.venv/` 这类目录都被
git ignore，可能出现在你的工作区里，但不属于仓库跟踪的源码布局。

## 仓库结构

```text
.
├── CMakeLists.txt
├── README.md
├── README_zh.md
├── LICENSE
├── bench/                      # benchmark 程序与 calc_iopp_params
├── include/                    # 公共头文件
├── src/                        # 库实现
├── tests/                      # 单元测试、CLI 测试和路径卫生检查
├── scripts/
│   ├── run_release_c4_lambda128.sh
│   ├── run_backend_eval_single_thread_7contexts.sh
│   └── plot_benchmark_results.py
├── FRI_Ligero-based_results.md # 跟踪的 FRI/Ligero baseline 表格
├── results-new/                # 当前格式的 CSV 与图
└── results-legacy/             # 历史 CSV 与旧图
```

当前仓库中跟踪的数据资产主要有：

- `results-new/1-thread/backend_eval_results.csv`
- `results-new/8-thread/backend_eval_results.csv`
- `results-new/8-thread/compiler_eval_results.csv`
- `results-new/fri_ligero_based_eval_results.csv`
- `results-new/figures/*.png`
- `results-legacy/` 中保留的旧 CSV 命名和 `*_vs_d.png` 图

## Benchmark 工作流

`scripts/run_release_c4_lambda128.sh` 完全通过环境变量配置。当前支持的
suite 有：

- `basefold_release`
- `compiler_eval_ring_switch`
- `compiler_eval_frobenius`
- `compiler_outer_commit_ring_switch`
- `compiler_outer_commit_frobenius`

当前默认值：

- `C=4`
- `K0=1`
- `LAMBDA=128`
- `D_MIN=3`，`D_MAX=29`
- `COMMIT_WARMUP=1`，`COMMIT_REPS=3`
- `EVAL_WARMUP=1`，`EVAL_REPS=3`
- `BENCH_THREADS=8`
- `CONTEXTS=all`（当前共 16 个 context）

当前 context id：

- `field-255`
- `ring-gr-2p16-162`
- `field-f2p256`
- `ring-gr-2p2-162`
- `field-prime64-ext`
- `field-f2p64-ext`
- `field-prime128-ext`
- `field-f2p128-ext`
- `field-f3p40-ext`
- `field-f3p81-ext`
- `ring-gr-2p16-64-ext`
- `ring-gr-2p16-128-ext`
- `ring-gr-2p32-64-ext`
- `ring-gr-2p32-128-ext`
- `ring-gr-2p2-64-ext`
- `ring-gr-2p2-128-ext`

当 `BENCH_THREADS > 0` 时，脚本还会把下面三个环境变量一起设为同一个值：

- `OMP_NUM_THREADS`
- `BASEFOLD_MERKLE_MAX_THREADS`
- `BASEFOLD_VERIFY_QUERY_MAX_THREADS`

如果你想保留运行时默认线程策略，可显式设置 `BENCH_THREADS=0`。

代表性命令：

```bash
# 默认 BaseFold release sweep
scripts/run_release_c4_lambda128.sh

# 只跑一个 BaseFold context，缩小 d 区间
CONTEXTS=field-prime128-ext \
D_MIN=10 D_MAX=20 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh

# 在 GR(2^32,64) 上跑 Ring-switch compiler full eval
RUN_SUITE=compiler_eval_ring_switch \
CONTEXTS=ring-gr-2p32-64-ext \
COMPILER_KAPPA=6 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh

# 在 GR(2^32,128) 上跑 Frobenius compiler full eval
RUN_SUITE=compiler_eval_frobenius \
CONTEXTS=ring-gr-2p32-128-ext \
COMPILER_KAPPA=7 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

当前需要注意的约束：

- compiler suite 只接受 ring context。
- 所选 ring context 必须满足 `deg(F) = 2^COMPILER_KAPPA`。
- 当前 compiler bench 只支持 `K0=1`。
- Frobenius compiler suite 会跳过 `GR(2^2,r)` context，并写出
  `status=disabled_gr2p2_context`。
- sweep CSV 里的 `queries` 是通过 `calc_iopp_params` 输出中的
  `l_min_for_IOPP` 解析得到的。

常用环境变量：

- `CONTEXTS`：选择 context 子集
- `COMPILER_KAPPA`、`COMPILER_ELL_MIN`、`COMPILER_ELL_MAX`：compiler suite
  参数
- `OUT_DIR`、`OUT_ROOT`、`RUN_ID`、`TIMESTAMP`：输出路径控制
- `BUILD_DIR`、`ISOLATE_BUILD_DIR`、`PIN_BUILD`：构建目录控制
- `CPU_PIN_MODE`、`CPU_SET`、`RUN_SLOT`、`RUN_SLOTS_TOTAL`、
  `USE_SMT_IN_SLOT`：绑核控制
- `CMD_TIMEOUT_SEC`、`CONTINUE_ON_ERROR`：失败处理

更完整的 context 参数表和环境变量说明见
[bench/exp_params_release_c4_lambda128.md](bench/exp_params_release_c4_lambda128.md)。

### 并发单线程 Backend Sweep

[scripts/run_backend_eval_single_thread_7contexts.sh](scripts/run_backend_eval_single_thread_7contexts.sh)
会并发启动 7 个单线程 `basefold_release` 子任务，把每个子任务固定到一个 CPU，
最后合并出一个总的 `backend_eval_results.csv`。

默认示例：

```bash
scripts/run_backend_eval_single_thread_7contexts.sh
```

这个辅助脚本依赖 `bash`、`lscpu`、`python3` 和 `taskset`。

## 输出文件

默认输出目录是 `results/release_c4_lambda128_sweep_<RUN_ID>/`。根据
`RUN_SUITE` 不同，目录中可能包含：

- `backend_eval_results.csv`：`basefold_release` 的输出
- `compiler_eval_results.csv`：`compiler_eval_ring_switch` /
  `compiler_eval_frobenius` 的输出
- `compiler_outer_commit_results.csv`：
  `compiler_outer_commit_ring_switch` /
  `compiler_outer_commit_frobenius` 的输出
- `RESULTS.md`：markdown 汇总表
- `logs/*.log`：`calc_iopp_params` 和各 bench 程序的原始日志

当前 `backend_eval_results.csv` 的主要指标列包括：

- `commit_mean_ms`
- `open_mean_ms`
- `prove_mean_ms`
- `verifier_mean_ms`
- `proof_size_kb`
- `proof_size_bytes`

当前 `compiler_eval_results.csv` 的主要指标列包括：

- `outer_commit_mean_ms`
- `backend_commit_mean_ms`
- `commit_total_mean_ms`
- `open_total_mean_ms`
- `prove_total_mean_ms`
- `verify_total_mean_ms`
- `outer_proof_size_kb`
- `total_proof_size_kb`

## 绘图

`scripts/plot_benchmark_results.py` 可以读取一个或多个 benchmark CSV，
并绘制 `commit`、`open`、`prover`、`verifier`、`proof_size` 相对于
约束数量的图。

如果不传 `-o`，脚本默认把 PNG 写到 `result/plots`。若想重生成仓库里跟踪的图，
请显式传 `-o results-new/figures`。

示例：

```bash
# 重生成当前跟踪的 1-thread backend 图
python3 scripts/plot_benchmark_results.py \
  results-new/1-thread/backend_eval_results.csv \
  -o results-new/figures \
  --prefix PCSoverGaloisRing_1_thread

# 重生成当前跟踪的 8-thread compiler 图
python3 scripts/plot_benchmark_results.py \
  results-new/8-thread/compiler_eval_results.csv \
  -o results-new/figures \
  --prefix PCSoverZ2K_8_threads
```

常用绘图选项：

- `--metrics commit prover verifier proof_size`
- `--proof-size-column <column>`
- `--legend-mode auto|inline|split`
- 多输入模式，可把多个 CSV 叠加到一张图里
