# BasefoldOverGR

基于 [NTL (Number Theory Library)](https://libntl.org/) 的一组 C++ 代码，用于在 **Galois Ring** `GR(p^k, s)`（直观上可理解为“模 `p^k` 的系数环上，再做次数为 `s` 的多项式扩张”）上实现/辅助一些常用计算：求逆、Hensel 提升、插值，以及不同表示之间的转换。

> 说明：NTL 的 `ZZ_p`/`ZZ_pX`/`ZZ_pE` 等类型依赖全局模数上下文（例如 `ZZ_p::init(mod)`、`ZZ_pE::init(F)`）。调用本项目函数前请确保相关上下文已正确初始化。

此外，本仓库还包含一份 `(c, k0, d)`-foldable linear code 的编码过程实现，支持在有限域 `F_{p^s}` 与 Galois ring `GR(p^k,s)` 上运行（见 `include/BaseFold/FoldableCode.hpp`）。该编码实现可用于 BaseFold 论文中的 IOPP 与 PCS 构造。

同时，仓库中也实现了 BaseFold 论文里的：

- BaseFold IOPP（folding + query，一并支持有限域与 Galois ring，见 `include/BaseFold/IOPP.hpp`）。
- 一个最小化的、基于 **Merkle + Fiat–Shamir** 的非交互 BaseFold PCS 单点求值证明（目前实现的是 `k0==1` 版本，见 `include/BaseFold/BaseFoldPCS.hpp`）。

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
│       └── BaseFoldPCS.cpp
├── bench
│   └── bench_pcs_commit.cpp
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
  - `FindPrimitiveElement(ZZ p, long k, long s)`：尝试在 `GR(p^k,s)` 上构造一个候选的本原元素。

### `src/GaloisRing/PrimitiveElement.cpp`

- `FindPrimitiveElement` 的实现：
  - 先在 `Z_p[x]` 中给定一个（需人工指定系数的）多项式 `F`，再用 `ZZ_pE::init(F)` 初始化扩张。
  - 以 `b = x (mod F)` 为基础返回 `b^(p^{k-1})`（常见于 Teichmüller 代表/单位根相关构造）。
  - 注意：当前 `F` 的系数是硬编码示例，实际使用时应替换为与你的 `p,s` 匹配的不可约/primitive 多项式。

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

### `include/BaseFold/BaseFoldPCS.hpp` / `src/BaseFold/BaseFoldPCS.cpp`

- 一个最小化的非交互 BaseFold PCS 单点求值证明（Protocol 4 + Merkle + Fiat–Shamir）：
  - `BaseFoldPCSCommit/ProveEval/VerifyEval`
  - 目前实现限制为 `k0==1`（对应 `Basefold_over_GR.pdf` 的 Protocol 4 版本）。

### `bench/bench_pcs_commit.cpp`

- 编码性能基准：测量 **去掉校验后的纯编码时间**（`EncodeFoldableUnchecked`），支持由命令行分别指定有限域与 Galois ring 的上下文参数。

## 依赖

- NTL（以及其底层依赖 GMP）。编译/链接方式因系统环境而异。

## 测试

项目提供如下测试：

- `tests/test_galois_ring_basic.cpp`：覆盖主要工具函数、求逆、插值，以及 Hensel 提升/本原元素的基本 smoke test。
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

构建后可运行：

```bash
./build/bench_pcs_commit --help
```

示例：

```bash
# GF(2^2) with F(x)=x^2+x+1 and zeta=x
./build/bench_pcs_commit --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --d 16 --reps 5 --warmup 1

# GR(4,2) with the same extension polynomial and zeta=x
./build/bench_pcs_commit --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 16 --reps 5 --warmup 1
```
