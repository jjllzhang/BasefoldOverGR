# BasefoldOverGR

基于 [NTL (Number Theory Library)](https://libntl.org/) 的一组 C++ 代码，用于在 **Galois Ring** `GR(p^k, s)`（直观上可理解为“模 `p^k` 的系数环上，再做次数为 `s` 的多项式扩张”）上实现/辅助一些常用计算：求逆、Hensel 提升、插值，以及不同表示之间的转换。

> 说明：NTL 的 `ZZ_p`/`ZZ_pX`/`ZZ_pE` 等类型依赖全局模数上下文（例如 `ZZ_p::init(mod)`、`ZZ_pE::init(F)`）。调用本项目函数前请确保相关上下文已正确初始化。

## 目录结构

```
.
├── include
│   └── GaloisRing
│       ├── utils.hpp
│       ├── Inverse.hpp
│       ├── HenselLift.hpp
│       └── PrimitiveElement.hpp
├── src
│   └── GaloisRing
│       ├── utils.cpp
│       ├── Inverse.cpp
│       ├── HenselLift.cpp
│       └── PrimitiveElement.cpp
└── tests
    └── test_galois_ring.cpp
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

## 依赖

- NTL（以及其底层依赖 GMP）。编译/链接方式因系统环境而异。

## 测试

项目提供了一个简单的测试程序 `tests/test_galois_ring.cpp`，覆盖了主要工具函数、求逆、插值，以及 Hensel 提升/本原元素的基本 smoke test。

在安装好 NTL/GMP 后，可在项目根目录编译并运行：

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -Iinclude \
  tests/test_galois_ring.cpp \
  src/GaloisRing/utils.cpp src/GaloisRing/Inverse.cpp src/GaloisRing/HenselLift.cpp src/GaloisRing/PrimitiveElement.cpp \
  -lntl -lgmp -lm -o test_galois_ring

./test_galois_ring
```

也可以使用 CMake（推荐 out-of-source 构建）：

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```
