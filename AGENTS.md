# AGENTS.md - NBODY6++GPU 开发指南

本文档为 AI 代理和开发者提供 NBODY6++GPU 项目的全面概述，帮助理解代码结构、开发流程和最佳实践。

---

## 1. 项目概述

**NBODY6++GPU** 是一个用于天体物理数值模拟的 N 体星团演化代码，由 Rainer Spurzem 团队维护。

### 核心特性
- **直接 N 体积分**：使用 Hermite Scheme 和 Block time step
- **GPU 加速**：利用 CUDA 进行规则力和势能计算
- **混合并行化**：MPI + OpenMP + GPU + SIMD (SSE/AVX)
- **恒星演化**：SSE/BSE 模型，包含质量损失、合并、潮汐效应
- **双星物理**：KS 正则化、公共包层演化、引力波反冲
- **星系环境**：模拟星团受绕转星系的潮汐力，包括多种可选模型
- **灵活输出**：支持 HDF5 / 二进制 / ANSI 格式，包含较少的 Python 读取接口

### 适用规模
- **粒子数**：任意
- **最佳性能**：N > 50,000（GPU 加速生效）

### 科学应用
- 球状星团演化
- 开放星团动力学
- 核星团模拟
- 银河系潮汐效应研究
- 引力波源形成

---

## 2. 仓库结构

```
nbody-fork/
├── src/
│   ├── Main/          # 主要源代码 (~500+ Fortran/CUDA 文件)
│   └── Tools/         # 工具程序 (如 MWpotential.f)
├── include/           # 头文件和参数配置 (~25 文件)
│   ├── params.h       # 最大粒子数等参数配置
│   ├── common6.h      # 全局公共块和物理常数
│   ├── kspars.h       # KS 正则化参数
│   └── cuda_*.h       # GPU 接口头文件
├── extra_inc/         # 额外包含文件（MPI、CUDA）
├── build/             # 构建输出目录
├── examples/          # 示例输入文件和 Jupyter 分析脚本
│   ├── input_files/   # 测试输入 (N10k, N100k)
│   └── *.ipynb        # 数据分析 Jupyter 笔记本
├── doc/               # 文档
│   └── nbody6++_manual.txt  # 指向 Overleaf 手册
├── macro/             # Autoconf M4 宏
├── config.examples/   # 配置示例
├── configure.ac       # Autoconf 配置脚本
├── Makefile.in        # Makefile 模板
├── README.md          # 用户文档
└── .github/
    └── workflows/     # GitHub Actions CI/CD
```

---

## 3. 编程语言和技术栈

### 主要语言

| 语言 | 用途 | 文件数量 | 关键文件 |
|------|------|----------|----------|
| **Fortran 77** | 核心 N 体积分、恒星演化 | ~400+ | `nbody6.F`, `intgrt.F`, `setup.F`, `hrplot.F` |
| **CUDA C++** | GPU 加速力计算 | 3 | `gpunb.gpu.cu`, `gpupot.gpu.cu`, `gpunb.velocity.cu` |
| **C++** | SIMD 矢量化 (SSE/AVX) | ~15 | `pot.sse.cpp`, `reg.avx.cpp`, `irr.sse.cpp` |
| **Python** | 数据分析、输入转换 | ~5 | Jupyter 笔记本, `restore_mtime.py` |
| **Shell** | 构建脚本 | ~3 | `configure`, `@input-transform` |

### 技术栈
- **构建系统**: GNU Autoconf/Automake
- **编译器**: GNU Fortran (gfortran), NVCC (CUDA), g++ (C++)
- **并行库**: OpenMPI, OpenMP
- **GPU 框架**: CUDA 10.0+
- **数据格式**: HDF5, 二进制, 纯文本
- **版本控制**: Git

---

## 4. 核心组件详解

### 4.1 主程序循环 (`src/Main/nbody6.F`)
- 主积分循环
- 时间步管理
- 输出控制
- MPI 进程协调

### 4.2 力计算模块
| 模块 | 文件 | 功能 |
|------|------|------|
| **规则力** | `regcor_gpu.F`, `nbint.F` | 邻居列表、规则力修正 |
| **不规则力** | `fpoly*.F` | 不规则力多项式外推 |
| **GPU 加速** | `gpunb.gpu.cu`, `gpupot.gpu.cu` | GPU 势能和力计算 |
| **SIMD 优化** | `pot.{sse,avx}.cpp`, `reg.{sse,avx}.cpp` | CPU 矢量化力计算 |

### 4.3 双星和正则化
| 类型 | 关键文件 | 说明 |
|------|----------|------|
| **KS 双星** | `ksinit.F`, `kspreg.F`, `ksreg.F`, `ksterm.F` | Kustaanheimo-Stiefel 正则化方法 |
| **链系统** | `chain.f`, `chlist.f`, `chinit.f` | 三体及更高层次系统 |
| **双星演化** | `binev.f`, `binpop.F`, `tides.f` | 双星演化、潮汐圈化 |

### 4.4 恒星演化
| 组件 | 文件 | 模型 |
|------|------|------|
| **单星演化** | `hrplot.F`, `star.f` | Hurley et al. SSE |
| **双星演化** | `binev.f`, `mix.f` | Hurley et al. BSE |
| **公共包层** | `comenv.f` | 包层演化物理 |
| **巨星物理** | `giant.f`, `giant2.f`, `giant3.f` | 巨星阶段质量损失 |
| **合并** | `coal.f`, `merge.f`, `merge2.f` | 碰撞和并合事件 |
| **GW 反冲** | `kickGW.f`, `recoil.f` | 引力波合并后的反冲速度 |

### 4.5 银河系势场
| 文件 | 功能 |
|------|------|
| `mwpotinit.F` | MWPotential2014 初始化 |
| `fmwpot.f` | 银河系势场力计算 |
| `pmwpot.f` | 银河系势能计算 |
| `MWpotential.f` | 工具函数 (src/Tools) |

### 4.6 输入/输出
| 类型 | 文件 | 格式 |
|------|------|------|
| **标准输出** | `output.F` | 纯文本或fortran二进制格式 |
| **HDF5 输出** | `custom_output.F`, `custom_output_facility.F` | HDF5 快照文件 |
| **诊断输出** | `energy.F`, `lagr.f`, `binout.f`, `stdout` | 能量、拉格朗日半径、双星统计 |

---

## 5. 构建系统

### 5.1 配置步骤

```bash
./configure [选项]
```

#### 关键配置选项

| 选项 | 说明 | 建议值 |
|------|------|--------|
| `--enable-mcmodel=ARG` | 内存模型 (small/medium/large) | `large` |
| `--with-par=ARG` | 最大粒子数 (1k/10k/100k/1m/b1m/b4m/b10m) | `b1m` (100万) |
| `--disable-gpu` | 禁用 GPU | 仅当 N<50k 或无 GPU 时 |
| `--enable-simd=ARG` | SIMD 优化 (sse/avx/no) | `avx` (如支持) |
| `--enable-hdf5` | HDF5 输出 | 推荐但需手动配置 |
| `--disable-mpi` | 禁用 MPI | **不推荐**（仅调试用） |
| `--enable-debug` | 调试模式 | 开发时使用 |
| `--prefix=PATH` | 安装路径 | 自定义 |

#### 快速开始配置

```bash
# 个人电脑测试（无GPU，无需MPI）
./configure --enable-mcmodel=large --with-par=b1m --disable-gpu --disable-mpi

# 典型生产环境
./configure --enable-mcmodel=large --with-par=b1m
```

### 5.2 编译

```bash
make clean && make -j
```

生成的可执行文件位于 `build/`，命名格式：
```
nbody6++.[simd].[mpi].[gpu].[hdf5]
```

例如：`nbody6++.avx.mpi.gpu`

### 5.3 安装

```bash
make install
```

将可执行文件安装到 `--prefix` 指定的路径。

---

## 6. 开发工作流程

### 6.1 获取代码

```bash
# 稳定版 (stable 分支)
git clone https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing

# 开发版 (dev 分支，包含最新更新和 bugfix)
git clone -b dev https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing
```

### 6.2 开发流程

1. **创建特性分支**
   ```bash
   git checkout -b feature/your-feature
   ```

2. **修改代码**
   - 源代码位于 `src/Main/`
   - 修改参数配置：`include/params.h`
   - GPU 代码：`.cu` 文件
   - SIMD 代码：`.cpp` 文件

3. **重新配置和编译**
   ```bash
   ./configure [选项]
   make clean && make -j
   ```

4. **测试**
   - 使用 `examples/input_files/N10k_noDat10.inp` 快速测试
   - 检查输出文件完整性

5. **提交代码**
   ```bash
   git add .
   git commit -m "描述性的提交信息"
   git push origin feature/your-feature
   ```

### 6.3 调试技巧

#### Fortran 调试
```bash
# 使用调试模式编译
./configure --enable-debug --disable-gpu --disable-mpi
make clean && make -j

# 使用 gdb
gdb build/nbody6++.debug
(gdb) run < examples/input_files/N10k_noDat10.inp
```

#### GPU 调试
```bash
# 使用 cuda-memcheck
cuda-memcheck ./build/nbody6++.avx.mpi.gpu < input.inp

# 使用 cuda-gdb
cuda-gdb ./build/nbody6++.avx.mpi.gpu
```

#### 内存检查
```bash
# 设置环境变量
export OMP_STACKSIZE=4096M
ulimit -s unlimited
```

---

## 7. 测试和示例

### 7.1 示例输入文件

| 文件 | 粒子数 | 模拟时长 | 用途 |
|------|--------|----------|------|
| `N10k_noDat10.inp` | 10,000 | 2 Myr | 快速测试（Plummer 模型） |
| `N100k.inp` | 100,000 | 1 Gyr | 生产运行示例 |
| `dat.10` | - | - | 预生成初始条件 |

### 7.2 运行示例

```bash
# 设置环境变量
export OMP_STACKSIZE=4096M
export OMP_NUM_THREADS=16  # 根据 CPU 核心数调整，最大32
ulimit -s unlimited

# 复制输入文件
cp examples/input_files/N10k_noDat10.inp ./

# 运行模拟
./build/nbody6++.avx.mpi.gpu < N10k_noDat10.inp
```

### 7.3 输出文件

| 文件 | 内容 |
|------|------|
| `OUT3` | 主输出文件（轨道、能量、碰撞等） |
| `snap.40_*` | HDF5 快照文件 |
| `sev.83_*` | 恒星演化数据 |
| `lagr.7` | 拉格朗日半径 |
| `bev.82` | 双星演化事件 |

### 7.4 数据分析

使用 Jupyter 笔记本分析输出。可用笔记本：
- `01_Basics.ipynb` - 基础数据读取
- `02_Hertzsprung–Russell_diagram.ipynb` - HR 图绘制
- `03_HDF5_Basics.ipynb` - HDF5 文件处理
- `readhdf5.ipynb` - HDF5 单文件读取示例

---

## 8. CI/CD 流程

### GitHub Actions 工作流

位于 `.github/workflows/`：

| 工作流 | 触发条件 | 功能 |
|--------|----------|------|
| `autotest.yml` | Push/PR | 自动测试构建 |
| `release_on_tag.yml` | 标签推送 | 发布版本 |
| `sync-stardisk.yml` | 定时/手动 | 同步特定分支 |

### 测试状态
[![autotest status](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/actions/workflows/autotest.yml/badge.svg)](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/actions/workflows/autotest.yml)

---

## 9. 代码规范和最佳实践

### 9.1 Fortran 代码规范

- **格式**: Fortran 77 固定格式（6 列缩进）或 极少数文件使用 Fortran 90 自由格式
- **大小写**: 混合使用（历史代码为大写，新代码可用小写）
- **注释**: 使用 `C` 或 `!` 开头
- **公共块**: 通过 `COMMON` 块共享变量（见 `include/common6.h`）
- **预处理器**: 使用 `.F` 扩展名启用预处理器指令（如 `#ifdef`）

### 9.2 命名约定

| 类型 | 约定 | 示例 |
|------|------|------|
| 子程序 | 小写，描述性名称 | `intgrt`, `ksinit`, `output` |
| 函数 | 大写，简短 | `ENERGY`, `LAGR` |
| 变量 | 大写（公共块），小写（局部） | `TIME`, `TPHYS`, `body(i)` |
| 参数 | 大写 | `NMAX`, `KMAX` |

### 9.3 性能优化指南

1. **向量化**: 使用 SSE/AVX SIMD 指令
2. **GPU 卸载**: 将规则力计算移至 GPU
3. **OpenMP 并行**: 在力计算循环中使用 `!$OMP PARALLEL DO`
4. **邻居列表**: 优化邻居搜索半径
5. **时间步控制**: 调整 `DTMIN`, `DTMAX` 参数

### 9.4 常见陷阱

1. **内存溢出**: 
   - 使用 `--enable-mcmodel=large`
   - 设置 `export OMP_STACKSIZE=4096M`
   - 执行 `ulimit -s unlimited`

2. **GPU 效率低下**:
   - N < 50k 时禁用 GPU (`--disable-gpu`)
   - 多 GPU 节点需手动绑定

3. **双星过多**:
   - 双星比例 >5% 时性能显著下降
   - 考虑调整 `KMAX` 参数

4. **OpenMP 线程数**:
   - 不要超过 32 个线程
   - 避免使用全部逻辑核心（超线程）

---

## 10. 文档资源

### 官方文档
- **完整手册**: [Overleaf Manual](https://www.overleaf.com/read/hcmxcyffjkzq)
- **README**: [GitHub README](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/blob/stable/README.md)
- **讨论区**: [GitHub Discussions](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions)

### 参考文献

| 主题 | 文献 | NASA ADS |
|------|------|----------|
| NBODY1-6 综述 | Aarseth 1999 | [1999PASP..111.1333A](https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/) |
| NBODY6++ | Spurzem 1999 | [1999JCoAM.109..407S](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/) |
| SSE/BSE | Hurley et al. 2005 | [2005MNRAS.363..293H](https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/) |
| NBODY6GPU | Nitadori & Aarseth 2012 | [2012MNRAS.424..545N](https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/) |
| NBODY6++GPU | Wang et al. 2015 | [2015MNRAS.450.4070W](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/) |
| 恒星演化更新 | Kamlah et al. 2022 | [2022MNRAS.511.4060K](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/) |

### 相关工具
- **dragon3_pipeline**: [完整的、大量的数据分析Python脚本（建设中）](https://github.com/kaiwu-astro/dragon3_pipeline)
- **McLuster**: [初始条件生成器](https://github.com/agostinolev/mcluster)
- **Galpy**: [银河系势场库](https://github.com/jobovy/galpy)

---

## 11. 已知问题

| 问题 | 状态 | 解决方案 |
|------|------|----------|
| 多 GPU 节点 GPU 分配 | 进行中 | 手动绑定 GPU bus ID |
| 双星比例 >5% 性能差 | 研究中 | 减少双星数量或使用更大硬件 |
| `--disable-simd` + `--enable-omp` 不工作 | 进行中 | 使用 SSE/AVX |
| 许多参数编译时硬编码 | 改进中 | 将来使用完整 Namelist 输入 |
| `KZ(7) >= 4` 输出错误 | 已知 | 使用 `KZ(7) <= 3` |
| HDF5 配置选项失效 | 已知 | 手动编辑 `build/Makefile` |

---

## 附录 A: 快速参考

### 环境变量

```bash
export OMP_STACKSIZE=4096M       # OpenMP 栈大小
export OMP_NUM_THREADS=16        # OpenMP 线程数 (<=32)
ulimit -s unlimited              # 栈大小无限制
```

### 常用命令

```bash
# 配置
./configure --enable-mcmodel=large --with-par=b1m --enable-simd=avx

# 编译
make clean && make -j

# 运行
./build/nbody6++.avx.mpi.gpu < input.inp

# 安装
make install
```

### 重要文件路径

- 主程序: `src/Main/nbody6.F`
- 参数配置: `include/params.h`
- 公共块: `include/common6.h`
- GPU 代码: `src/Main/gpunb.gpu.cu`, `gpupot.gpu.cu`
- 示例输入: `examples/input_files/`

---

**最后更新**: 2026-01-20  
**文档版本**: 1.0  
**适用代码版本**: NBODY6++GPU Beijing (stable/dev)
