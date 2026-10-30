# FSILBM3D CUDA 并行说明

本文说明三件事：

1. 相对 CPU 版本（`src` 下除 `cuda_lbm.cuf` 外的原始模块），哪些步骤改成了 GPU 计算，算法/结果是否有差别，哪些部分仍由 CPU 并行完成；
2. CUDA 实现的线程、内存、stream/event 细节；
3. `cuda_lbm.cuf` 模块的内容，以及为接入 GPU 而对其它 `src` 模块所做的修改。

当前版本的关键特性：

- 求解核心（碰撞、传播、宏观量、边界、IBM、父子块通信）全部在 GPU 上；
- IBM 罚力使用「先算力、再修正」的两遍 kernel，与 CPU 的 `PenaltyForce_` 语义一致；
- GPU 宏观量按需计算（lazy macro），避免步末无条件计算及不必要的下载；CPU 路径保留每步结束后的统一宏观量更新；
- IBM 元素数组使用 pinned 主机内存 + 独立 upload stream 异步上传；
- flux、探针、流场输出均在 GPU 上打包，只回传实际需要的数据。

> 验证基线：在 3 block、细化比 2/4、BC 101/103/301、`iCollidModel=1` 的球体算例上，GPU 与 CPU（gfortran 构建）连续 100 步的 `Group001_forces.dat` 中 **Fx/Fy 逐位一致**，Fz 仅 `~1e-13`；`FluidFlux.dat` 完全一致，探针差 `~1e-12` 或更小。

---

## 1. GPU 化的步骤

### 1.1 逐阶段对照

| 阶段 | CPU 版本 | CUDA 版本 | 结果是否一致 |
|---|---|---|---|
| 宏观量 | `calculate_macro_quantities_`，OpenMP 按 x 并行 | `macro_kernel`；显式调用 `cuda_backend_macro`，按需调用 `cuda_backend_macro_if_needed` | 公式相同（含 `0.5*volumeForce*dh`），一致 |
| 体积力清零 / 叠加 | `ResetVolumeForce_`、`add_volume_force_` | `reset_force_kernel`、`add_uniform_force_kernel` | 一致 |
| 碰撞 | `collision_`，支持模型 1/2/3/11/12/13/14/15 | 模型 1：`collision_srt_kernel`（单遍标量 SRT）；模型 2：`collision_kernel`（TRT）；模型 3：`collision_mrt_kernel`；模型 11/14/15：`collision_les_kernel` | 已移植的 1/2/3/11/14/15 公式与 CPU 相同；12/13 不支持，见 1.2 |
| 半程反弹缓存 | `halfwayBCset_` 写入 `fIn_hwx*` 等主机数组 | `halfway_boundary_kernel`（mode=0）写入设备数组 `hwx*` | 一致 |
| 传播 | `streaming_`，按方向原地换位（`swapzy` + `swapx`） | `stream_kernel`，双缓冲 gather，周期邻居用 `mod` 折返 | 数学等价，见 1.2 |
| 边界条件 | `set_boundary_conditions_`，逐面循环 | `cuda_backend_boundary` 分派到 `general_boundary_kernel` / `static_boundary_kernel` / `halfway_boundary_kernel` / `moving_halfway_apply_kernel` | 支持的 BC 公式相同；当前算例的 101/103/301 已验证 |
| IBM 罚力迭代 | `PenaltyForce_`：先用未修正的 `u` 算完所有单元力，再统一修正 `u` | `ibm_penalty_force_kernel` + `ibm_penalty_apply_kernel`，同一条 compute stream 顺序执行 | 两遍语义与 CPU 一致 |
| IBM 力扩散 | `FluidVolumeForce_` 中的流体侧 `do x=-1,2 ...` | `ibm_spread_force_kernel` | 贡献集合一致；原子加顺序只带来舍入差 |
| 湍流统计 | `calculate_turbulent_statistic_` | `turbulent_statistic_kernel` | 公式相同；`invStep` 精度不同，见 1.2 |
| 父块 → 子块边界采样 | `extract_interpolate_layer`，父块 6 个面打包下载到主机 | `extract_face_kernel`，父块设备数组直接写入子块设备缓冲 | 一致，全程 D2D |
| 时间插值 | `extract_interpolate_layer` 中 `0.5*(t1+t2)` | `average_face_f_kernel` / `average_face_tau_kernel` | 一致 |
| 父块 → 子块空间插值 | `interpolation_father_to_son` 主机循环 + `fIn_GridTransform` | `interpolate_face_kernel` + `transform_f_lbm` | 一致 |
| 子块 → 父块 | `deliver_son_to_father` 主机循环 + `fIn_GridTransform` | `deliver_face_kernel` + `transform_f_lbm` | 一致 |
| 通量诊断 | `write_fluid_flux` 下载整场 `den/uuu` 后求和 | `pack_flux_planes_kernel` 只打包 3 个 x 截面的 `u/rho`（6 个平面），主机按相同顺序求和 | 数值一致；只减少传输量 |
| 探针 | `write_fluid_information` 使用整场速度插值 | `pack_probes_kernel` 一次打包所有 probe 的 2×2×2×3 模板，主机用相同系数/顺序插值 | 数值一致；只减少 launch/同步次数 |
| 流场输出 | 下载整场 `den/uuu`，主机计算 `outtmp` 后写盘 | `outputtype=1` 时由 `pack_flow_output_kernel` 在 GPU 上直接生成 `real(4)` 输出并回传；`outputtype>=2` 仍走主机路径 | 输出文件一致（早期输出逐位相同） |
| 输出 / 续算 | 主机数组直接写出 | 按需从设备下载后再写出；continue 读入后会重新上传 `fIn` | 一致 |

GPU 侧支持的边界条件为 `0, 1, 101, 102, 103, 104, 201, 202, 203, 204, 301, 302`，与 CPU 版本定义的全部边界类型一一对应。

### 1.2 与 CPU 算法**确定**有区别的地方

1. **碰撞模型范围**
   CPU 支持 SRT(1)、TRT(2)、MRT(3)、Smagorinsky(11)、RBGK(12)、ELBM(13)、WALE(14)、Vreman(15)。
   GPU 已实现：
   - 1：`collision_srt_kernel`（单遍标量 SRT，默认快路径）；
   - 2：`collision_kernel`（generic TRT）；
   - 3：`collision_mrt_kernel`，两个 19×19 的 `M_COLLID/M_FORCE` 矩阵由主机 `calculate_MRT_params` 计算后、每个 block 上传一次（`cuda_backend_set_mrt`），kernel 内做 19 项点积；
   - 11：`collision_les_kernel`，Smagorinsky 只用本点非平衡量，无需邻域；
   - 14/15：`collision_les_kernel`，WALE/Vreman 读取 `u` 的一阶差分；为与 CPU 逐位一致，保留了 CPU 的 `invdh=dh` 写法（不是 `1/dh`）和边界单侧差分规则；
   - 12：RBGK **未移植到 GPU**。CUDA 构建在 `FluidDomain` 的 `collision_` 中、调用 GPU 后端之前显式 `stop`，提示该模型不支持；
   - 13：ELBM **未移植到 GPU**，与模型 12 一样，在调用 GPU 后端之前显式 `stop`，提示该模型不支持。
   未支持的模型不会再有“静默跳过碰撞”的行为。

   验证：模型 1 与 CPU 对照 100 步；模型 2/3/11/14/15 对照 25 步。`Group001_forces` 的 **Fx/Fy 差值均为 0**，Fz 仅 `~1e-12`。

2. **湍流统计的 `invStep` 精度**
   CPU：`invStep = 1/real(step-step_s+1)`，按默认单精度计算后赋给 `real(8)`；
   GPU：`1.0d0/real(step-step_s+1,8)`，全程双精度。累计量在末几位会有差别（仅 `outputtype>=2` 使用）。

3. **浮点累加顺序**
   IBM 罚力的「施加修正」kernel（`ibm_penalty_apply_kernel`）和力扩散 kernel（`ibm_spread_force_kernel`）对 `u`、`force` 使用双精度 `atomicAdd`，累加顺序不确定；CPU 是固定顺序的 OpenMP/顺序循环。**单元力的计算已经拆成独立的第一遍 kernel，所有单元都基于修正前的 `u`，与 CPU 的两遍语义一致**；剩余的原子加顺序只造成 `1e-16` 级舍入差，不再产生系统性偏差。异步上传和 IBM 串行化 event 只会改变 GPU 调度，不改变算法。

> IBM 收敛判据已与 CPU 一致：GPU 现在每轮检查（`IBM_CHECK_INTERVAL=1`），不再每 4 轮才检查。Zhu 板算例（`ntolLBM=10`）验证：t≤1 的 `nodeAverage` 位移差由每 4 轮检查时的 ~3.5e-2 降到 ~2e-4；步长约 2.4–3.2 ms，代价可忽略。t≥2 后板振动的相位仍会因浮点/原子顺序的微小差异漂移，这是混沌敏感性问题，不是检查频率造成的系统性偏差。

> 已确认**没有**算法差别、只是实现方式不同的部分：传播（in-place 换位 vs gather）、BCSymmetric（`mirror_x/y/z` 与 CPU 镜像一致）、宏观量（`SUM` vs 顺序累加）、父子块通信（设备端实现）、probe/flux/flow 输出打包、IBM 面积上传/单元力清零方式。

### 1.3 仍在 CPU 上执行（OpenMP 并行）的部分

| 部分 | 说明 |
|---|---|
| block 树与调度 | `build_block_tree`、`tree_collision_streaming_IBM_FEM` 递归、`CommPair` 通信拓扑 |
| 体积力系数 | `update_volume_force_`，每个 block 计算一个标量值后再交给 GPU 叠加 |
| IBM 预处理 | `UpdateElmtInterp`（4×4×4 插值索引/权重）、`UpdateElmt*` 中的虚拟体位置/速度/面积更新 |
| 固体求解 | Newmark/FEM 求解器、节点载荷 `lodFlow` 计算、固体输出 |
| 边界调度 | 判断每个面用哪种边界、每轮 IBM 收敛判断、虚拟体槽位管理 |
| 输出与后处理 | `write_flow_` 的文件写入、`write_fluid_flux` 的主机求和、探针文件、`ComputeFieldStat`、续算文件 |
| 文件 I/O | 输入解析、结果文件写入 |

这些部分沿用原有的 `!$OMP PARALLEL DO` 指令，由编译器打开 OpenMP（nvfortran 用 `-mp`，gfortran 用 `-fopenmp`）。

> 注意：`set_boundary_conditions_`、`halfwayBCset_`、`interpolation_father_to_son`、`deliver_son_to_father` 中的主机循环在 CUDA 构建下不再执行（子程序开头直接走 GPU 分支返回），但保留给 CPU 构建使用。

### 1.4 CPU/GPU 对拍时的编译选项（FMA）

gfortran 的 CPU 构建默认面向通用 x86-64（`-O3`，未指定 `-march`），**不生成 FMA 指令**；nvfortran 的 GPU 构建在 `-O3` 下默认生成 FMA。实测默认 GPU 可执行文件含 1194 条 FMA，CPU 可执行文件为 0 条。两者每步会产生 ~1e-16 的舍入差；对敏感 FSI 算例，这个差会被迅速放大：

- Zhu 板算例，`ntolLBM=10, dtolLBM=1e-12`（每步固定 10 轮罚力）：
  - 只把 GPU 的 FMA 关掉、代码和参数完全不变，与默认 FMA 的 GPU 结果相比：t=0.75 的 |Δu|=0.117，t=1.0 为 0.419，t=1.5 为 0.232；
  - CPU 与默认 GPU 的差为：t=0.75 0.117，t=1.0 0.410，t=1.5 0.190。两者量级相同，说明差异来自算术求值，而不是算法或移植错误。
- 关掉 FMA 后（`make GPU=1 STRICT=1`），同一组参数下 CPU 与 GPU 的 |Δu| 降为：t=0.25 1.0e-7、t=0.5 1.3e-5、t=0.75 1.2e-4、t=1.0 1.4e-3、t=1.5 4.5e-3；t=0~4.4 的最大值为 1.4e-2（t=2.83），而 |u|≈1.5。逐步输出中第 10 步差 4e-14、第 64 步差 6.3e-9（默认 FMA 为 3.1e-7）。

CPU/GPU 对拍时，应明确区分两种构建，并在切换前清理旧编译产物：

```bash
# CPU 基线：CMP 应与 master 的对照构建一致，CPU 回归不要求 STRICT=1
make clean && make GPU=0 CMP=gfortran
# 保存 CPU 可执行文件或完成 CPU 测试后，再构建 GPU 对照版本
make clean && make GPU=1 STRICT=1
```

当前 Makefile 的 `CMP` 可以选择 CPU 编译器；上面显式选择 `gfortran` 是为了对应本节的 gfortran 基线。使用 ifort/ifx 对照时，新旧 CPU 版本也应使用相同编译器和参数，不能把 gfortran 的 FMA 行为直接套用到 ifort/ifx。

`STRICT=1` 在 nvfortran 构建中同时加 `-Mnofma`（主机代码）和 `-gpu=cc80,nofma`（设备代码），用于减少与 gfortran CPU 基线之间的 FMA 舍入差异，不保证逐位一致。对于 `GPU=0` 的 gfortran/ifort/ifx 路径，`STRICT=0/1` 不改变编译参数，因此无需添加 `STRICT=1`；若显式使用 `GPU=0 CMP=nvfortran`，该开关仍会影响其编译参数。生产计算仍使用默认 FMA 构建即可；本算例中两者单步耗时接近（10 轮罚力时约 13 ms/步）。

罚力迭代参数的选取建议：

- `dtolLBM=0.01`（原设置）时残差紧贴阈值，迭代次数会在 1/2 之间翻转（实测第 224 步首次出现 CPU/GPU 迭代次数不同，到 t=0.5 有 519/3200 步不同），不适合作为对拍设置。
- `dtolLBM=1e-12, ntolLBM=10` 时迭代次数恒为 10，计数确定；但罚力耦合更刚，对舍入差更敏感，分叉反而更早（t≈0.7）。
- **纯代码对拍推荐 `ntolLBM=1`**：每步固定 1 轮，GPU 两次运行 t=0..7 逐位相同；配合 `STRICT=1` 时分叉最晚。
- 物理生产不要用 `ntolLBM=1`（IBM 无滑移欠收敛）。生产可用 `ntolLBM=10, dtolLBM=1e-3~1e-4` 或其它"残差远离阈值"的设置，让迭代次数稳定。

> 即使关掉 FMA，CPU/GPU 仍会因求和顺序、原子加顺序、主机/设备编译器差异保留 ~1e-16 的种子；该 FSI 算例在 t≈1.5–2 之后仍会逐渐分离。长时间验证应比较统计量（振动幅值、频率/相位、平均或均方根力、能量），而不是逐点曲线。

---

## 2. CUDA 并行细节

### 2.1 线程分配

**块内一维展平（宏观量、碰撞、传播、体积力、输出打包）**

```
BLOCK_SIZE = 128
tid = (blockIdx%x - 1) * BLOCK_SIZE + threadIdx%x
z   = mod(tid - 1, zDim) + 1
y   = mod((tid - 1) / zDim, yDim) + 1
x   = (tid - 1) / (zDim * yDim) + 1
grid = ceil(zDim*yDim*xDim / BLOCK_SIZE)
```

- 数组按 `(z,y,x,q)` 列优先存放，相邻 `tid` 对应相邻 `z`；对同一个 `q`，一个 warp 内对 `f(z,y,x,q)`、`rho`、`u` 的访问地址连续，满足合并访问；
- 整个三维块只做一次尾部补齐（`tid > npoints` 的线程直接返回）；
- 一个 block 的线程可能跨越多个 `(x,y)` 行，所有下标由 `tid` 反算，不依赖 `threadIdx%x` 与 `z` 的对应关系；
- 该 `tid → (z,y,x)` 映射沿用原版本经实测确定的效率最优分配方式，CUDA 移植没有改动这一映射；`BLOCK_SIZE` 只是独立的启动参数。

**面型 kernel**（边界条件、半程反弹、父子块通信）

| 面 | b | a | 线程数 |
|---|---|---|---|
| 1, 2（x 面） | z | y | `zDim*yDim` |
| 3, 4（y 面） | z | x | `zDim*xDim` |
| 5, 6（z 面） | y | x | `yDim*xDim` |

展平方式为 `b = mod(tid-1,B)+1`、`a = (tid-1)/B+1`，B 为该面的第一维长度。

**IBM kernel**

- `ibm_penalty_force_kernel`：一个线程对应一个虚拟体单元，做 4×4×4 = 64 个点的 gather，计算 `delta`、`local_force`，累加 `elem_force` 并把本次增量写入 `elem_dforce`；
- `ibm_penalty_apply_kernel`：一个线程对应一个单元，读 `elem_dforce`，对 64 个点做 `atomicAdd` 修正 `u`；
- `ibm_spread_force_kernel`：一个线程对应一个单元，对 64 个点做 `atomicAdd` 叠加 `force`；
- 几何、速度、单元力常驻显存，迭代过程中不再从主机读取。

**父子块 kernel**

- `extract_face_kernel`：一个线程对应父块面上属于该子块的一个节点；
- `average_face_*`：一个线程对应一个面节点；
- `interpolate_face_kernel`：一个线程对应一个细网格面节点；
- `deliver_face_kernel`：一个线程对应一个粗网格节点（每隔 `m_gridDelta` 取一个节点）。

**打包 kernel**

- `pack_flux_planes_kernel`：一个线程对应 `(z,y)` 平面上的一个点，写 6 个平面值；
- `pack_probes_kernel`：一个线程对应一个 probe 的一个模板分量（`24*N` 个线程）；
- `pack_flow_output_kernel`：一个线程对应输出窗口中的一个网格点，写 `p/u/v/w` 的 `real(4)` 值；
- `zero_ibm_force_kernel` / `zero_scalar_kernel`：设备端清零。

### 2.2 共享内存

| kernel | 共享内存 | 说明 |
|---|---|---|
| `collision_kernel` / `collision_srt_kernel` | `block_rho(128)`、`block_u(128,3)`、`block_force(128,3)`，约 7 KB | 每个线程只读写自己下标对应的元素，没有跨线程复用；`syncthreads()` 只是把数据搬到共享内存再读回。这一步可以去掉以进一步提高占用率 |
| `stream_kernel` | 无 | 见下 |
| 其它 kernel | 无 | — |

`stream_kernel` 不使用共享内存 halo，原因：

1. 展平后的 `threadIdx%x` 并不总等于 `z`（当 `zDim` 不能整除 `BLOCK_SIZE`，或一个 block 跨越多个 `(y,x)` 行时），`halo(local_z±1)` 会取到其它行的数据，结果错误；
2. halo 数组会显著增加共享内存占用、压低占用率。

去掉 halo 后，z 方向邻居依然由相邻线程访问相邻地址，硬件合并访问即可覆盖。

### 2.3 显存分配与复用

**每个 fluid block 常驻的设备数组**（`cuda_backend_init` 中一次性分配，整个运行期不释放）：

| 数组 | 大小 |
|---|---|
| `f`, `f_next` | `19 * N * 8` 字节，各一份 |
| `rho` | `N * 8` |
| `u` | `3 * N * 8` |
| `force` | `3 * N * 8` |
| `uuu_ave` | `9 * N * 8`（仅 `outputtype >= 2` 时实际使用） |
| `tau_all` | `N * 8` |
| `d_tol` | 1 个 `real(8)`（IBM 收敛容差累加器） |
| `mrt_collid`, `mrt_force` | `19×19` 个 `real(8)`（MRT 矩阵，每个 block 一份，约 5.8 KB） |
| `probebuf` / `probe_idx` | `8*3*Nprobe` 个 `real(8)` / `3*Nprobe` 个 `integer`（按需分配） |
| `fluxbuf` | `6*zDim*yDim` 个 `real(8)` |
| `flowbuf` | 输出窗口的 `real(4)` 缓冲（仅 `outputtype=1` 使用，按需分配） |

其中 `N = zDim*yDim*xDim`。每点合计约 368 字节（含统计量时约 440 字节）。

**双缓冲**：`f` 和 `f_next` 用 `device, pointer` 声明，传播结束后在主机侧交换指针（`tmp=>state%f; state%f=>state%f_next; ...`），没有任何数据拷贝。

**IBM 设备缓冲**：每个 carrier block 维护一个 `cuda_ibm_slot` 数组，按承载的虚拟体数量分配槽位；每个槽位包含：

- `elem_index`（`integer(2)`，12×nElem）
- `elem_weight`（`real(4)`，12×nElem）
- `elem_area`（`real(8)`，nElem）
- `elem_velocity`（`real(8)`，3×nElem）
- `elem_force`（`real(8)`，3×nElem，跨罚力迭代累加）
- `elem_dforce`（`real(8)`，3×nElem，本轮增量，仅设备使用）
- `area_loaded` 标志（刚体只上传一次面积）

- 几何与速度每个时间步只上传一次（`cuda_backend_ibm_prepare`），不在罚力迭代内部上传；
- `elem_area` 对刚体只在首次/重新分配时上传；对弹性体（`iBodyModel=2`）每个时间步上传一次；
- 单元力在设备上清零（`zero_ibm_force_kernel`），迭代中在设备上累积，全部结束后由 `cuda_backend_ibm_download_force` 一次性回传；
- 收敛容差是设备标量，`cuda_backend_ibm_tolerance` 同步取回并清零。

**父子块缓冲**：每个子块按面维护 6 组 `cuda_face_source`（`f1/f2/tau1/tau2`），第一次使用时分配，之后不再变化。父块 → 子块的采样、时间平均、空间插值全部读写这些设备缓冲。

**半程反弹缓存**：`hwx1..hwz2` 共 6 个设备数组，只在对应面确实使用 203/204 边界时分配。

**主机内存**：`fIn`、`fIn_Fx1t1` 等主机数组仍然保留（CPU 构建和断点续算用）；CUDA 构建下父子块通信不再访问它们。IBM 的 `v_Ei/v_Ew/v_Ea/v_Evel/v_Eforce` 在 CUDA 构建下为 **pinned allocatable**，供异步上传使用。

### 2.4 主机 ↔ 设备传输

只在真正需要主机数据时发生：

| 传输 | 频率/方式 |
|---|---|
| 初始化上传 `f/rho/u/force/tau_all` | 每个 block 一次 |
| 宏观量下载（`den/uuu`） | 输出、探针、场统计、断点续算时；lazy macro 只有消费者需要时才计算/下载 |
| 分布函数下载（`f`） | 断点续算写盘时 |
| IBM 几何/速度上传 | 每个时间步、每个虚拟体一次；独立 upload stream + pinned，异步 |
| IBM 面积上传 | 刚体首次一次；弹性体每个时间步一次 |
| IBM 单元力清零 | 每个时间步在 **GPU 上**清零，不传零数组 |
| IBM 单元力下载 | 每个时间步、每个虚拟体一次 |
| 通量平面下载 | info 输出时，6 个 x 平面 |
| 探针模板下载 | info 输出时，`8*3*Nprobe` 个值一次 |
| 流场输出下载 | `outputtype=1` 时回传 `real(4)` 输出窗口；`outputtype>=2` 回传整场 rho/u |
| 湍流统计量下载 | 输出时（仅 `outputtype>=2`） |

块与块之间的父子插值**没有任何主机传输**：父块采样、时间平均、空间插值、fine-to-coarse 变换都在设备间完成。

### 2.5 stream 与事件

- 每个 block 有 compute stream 和 transfer stream；IBM 元素上传另有一个 upload stream；
- 宏观量、碰撞、传播、边界、IBM、父子块 kernel 都排在 compute stream 上；
- D2H/H2D 传输排在 transfer stream 上，传输前用 `cudaEventRecord` + `cudaStreamWaitEvent` 等待 compute stream，传输后用 `cudaEventSynchronize` 等待数据就绪；
- IBM 元素上传：
  1. upload stream 只等待 `ibm_event`（上一个 spread kernel，最后一个读这些缓冲区的 kernel）；
  2. 上传完成后记录 `upload_event`，compute stream 在 penalty 前等待它；
  3. 因此上传可以与上一步的流体 kernel 重叠，但不会早于/晚于错误的数据依赖；
- IBM 串行化：模块级 `ibm_serial_event` 在 penalty 前等待、在 spread 后记录，使各 block 的 penalty/spread 按主机处理顺序串行，降低原子加调度的随机性；
- 父子块之间跨 stream 的数据依赖由每个 state 上的 `comm_event` 交接：
  1. 父块采样 kernel 跑在父块 stream 上，结束后记录事件，子块插值 kernel 在子块 stream 上等待；
  2. `deliver_son_to_father` 先在子块 stream 上记录事件，父块 stream 等待后才执行 `deliver_face_kernel`；
- 依赖边始终指向更早入队的事件，不会成环；时间步内不需要 `cudaDeviceSynchronize`。

---

## 3. `cuda_lbm.cuf` 模块简介

`src/cuda_lbm.cuf` 定义模块 `cuda_lbm_backend`，是整个 CUDA 后端。主要内容：

### 3.1 数据类型

| 类型 | 内容 |
|---|---|
| `cuda_lbm_state` | 一个 fluid block 的设备状态：`f/f_next`（双缓冲）、`rho/u/force/uuu_ave/tau_all`、边界缓存 `hwx*/hwy*/hwz*`、IBM 槽位 `ibm(:)`、父子块面缓冲 `facesrc(:)`、`d_tol`、打包缓冲 `probebuf/probe_idx/fluxbuf/flowbuf`、stream/event、`macro_valid` 等 |
| `cuda_ibm_slot` | 一个虚拟体的设备数据：`elem_index/elem_weight/elem_area/elem_velocity/elem_force/elem_dforce`、`nElem`、`loaded/area_loaded` |
| `cuda_face_source` | 一个子块面的父块采样缓冲：`f1/f2/tau1/tau2`、`bF/aF` |

### 3.2 kernel

| 类别 | kernel |
|---|---|
| 宏观量 | `macro_kernel`、`turbulent_statistic_kernel` |
| 碰撞 | `collision_kernel`（generic，模型 1/2）、`collision_srt_kernel`（模型 1 单遍标量）、`collision_mrt_kernel`（模型 3）、`collision_les_kernel`（模型 11/14/15） |
| 传播 | `stream_kernel` |
| 体积力 | `reset_force_kernel`、`add_uniform_force_kernel` |
| 边界 | `general_boundary_kernel`（101/102/202）、`static_boundary_kernel`（103/104/201/203/302）、`halfway_boundary_kernel`、`moving_halfway_apply_kernel` |
| IBM | `ibm_penalty_force_kernel`、`ibm_penalty_apply_kernel`、`ibm_spread_force_kernel` |
| 父子块 | `extract_face_kernel`、`average_face_f_kernel`、`average_face_tau_kernel`、`interpolate_face_kernel`、`deliver_face_kernel` |
| 打包/清零 | `pack_flux_planes_kernel`、`pack_probes_kernel`、`pack_flow_output_kernel`、`zero_scalar_kernel`、`zero_ibm_force_kernel` |
| device 辅助 | `center_diff`、`onesid_diff`（WALE/Vreman 的差分函数） |

### 3.3 主机端接口

| 类别 | 子程序 |
|---|---|
| 生命周期 | `cuda_backend_init`、`cuda_backend_finalize`、`cuda_backend_device_info` |
| 宏观量 | `cuda_backend_macro`、`cuda_backend_macro_batch`、`cuda_backend_macro_if_needed`、`cuda_backend_download_macro` |
| 流体步 | `cuda_backend_collision`、`cuda_backend_stream`、`cuda_backend_boundary`、`cuda_backend_halfway`、`cuda_backend_moving_halfway`、`cuda_backend_reset_force`、`cuda_backend_add_uniform_force` |
| 碰撞模型 | `cuda_backend_set_mrt`（上传 MRT 矩阵）；`cuda_backend_collision` 按 model 分派 1/2/3/11/14/15，其他模型由通用检查 `stop`；12/13 已在 `FluidDomain` 的 `collision_` 中提前拦截 |
| IBM | `cuda_backend_ibm_setup`、`cuda_backend_ibm_prepare`、`cuda_backend_ibm_penalty`、`cuda_backend_ibm_tolerance`、`cuda_backend_ibm_download_force`、`cuda_backend_spread_ibm_force`、`cuda_backend_ibm_clear` |
| 父子块 | `cuda_backend_extract_face`、`cuda_backend_average_face`、`cuda_backend_interpolate_face`、`cuda_backend_deliver_face`、`cuda_backend_record_event`、`cuda_backend_wait_event` |
| 诊断/输出 | `cuda_backend_download_flux_planes`、`cuda_backend_download_probes`、`cuda_backend_pack_flow_output`、`cuda_backend_turbulent_statistic` |
| 辅助 | `cuda_backend_download_f`、`cuda_backend_upload_f`、`cuda_boundary_supported`、`cuda_backend_ibm_clear_slot`、`cuda_backend_device_info` |

### 3.4 实现约定

- 所有物理量（`f/rho/u/force/tau`、IBM 面/速度/力）都是 `real(8)`；IBM 插值权重是 `real(4)`（与 CPU 相同）；流场输出是 `real(4)`（与 CPU 相同）；
- 一维展平 kernel 中一个线程对应一个网格点；IBM kernel 中一个线程对应一个单元；面 kernel 一个线程对应一个面节点；
- 罚力分两个 kernel：force 只读 `u`，apply 只写 `u`，与 CPU 的两遍循环一致；
- MRT 矩阵使用设备全局内存（每个 block 一份），不是 constant memory——不同 block 的 `Omega/Omega2/params` 可能不同，不能共用同一份常量；19×19 矩阵很小，访问会被 L1/常量缓存覆盖；
- LES 的 `collision_les_kernel` 保留 CPU 的局部数组写法（`feq/flb/old/fneq`），因此会比模型 1 的融合 SRT 慢；这是为与 CPU 逐位一致而做的取舍；
- 该 NVHPC 21.2 构建中，Fortran 类型化数组的 `cudaMemcpyAsync` 第三个参数是**元素个数**（内部再乘元素大小）；`cudaMemsetAsync` 对 device 数组无法解析，因此设备清零使用 `zero_*_kernel`；
- 所有 kernel 通过 `check_cuda(cudaGetLastError(), ...)` 检查启动错误；
- `cuda_backend_device_info` 在启动时打印 `CUDA_VISIBLE_DEVICES`、可见设备数，以及每个可见设备的 PCI bus:device。配合 `nvidia-smi --query-gpu=index,pci.bus_id` 可以确认进程实际使用的物理 GPU；例如 `CUDA_VISIBLE_DEVICES=7` 时打印 `CUDA devices: 1`、`PCI bus:device = 202:0`，对应 `nvidia-smi` 中 bus 为 `CA`（十六进制 202）的 GPU 7。

---

## 4. 其它 `src` 模块的修改

为接入 CUDA，只修改了以下文件（其余文件没有 GPU 专属改动）：

### 4.1 `FluidDomain.f90`

- `use cuda_lbm_backend`：引入 GPU 后端接口；
- `allocate_fluid_` / `initialise_`：`#ifdef USE_CUDA` 下调用 `cuda_backend_init`；
- `calculate_macro_quantities_`：GPU 分支调用 `cuda_backend_macro`，并根据 `cuda_boundary_supported` 决定是否回传 `den/uuu`；
- `calculate_macro_quantities_blocks`：初始化和需要时调用逐 block 的 `cuda_backend_macro`；
- 新增 `ensure_macro_blocks`：lazy macro 的消费者入口，只在 `macro_valid` 为假时重算；
- `initialise_`：模型 3 在 `cuda_backend_init` 之后调用 `cuda_backend_set_mrt` 上传 `M_COLLID/M_FORCE`；
- `collision_` / `streaming_` / `set_boundary_conditions_` / `halfwayBCset_`：`#ifdef USE_CUDA` 下直接调用对应 GPU 例程并返回，主机循环保留给 CPU 构建；
- GPU `collision_` 允许模型 3/11/14/15，并把 `this%tau`、`flow%nu` 传给后端；模型 12/13 在入口显式 `stop`；
- `write_flow_blocks`：GPU 路径先 `ensure_macro_blocks`；
- `write_flow_`：`outputtype=1` 时调用 `cuda_backend_pack_flow_output` 直接得到 `real(4)` 输出；`outputtype>=2` 保留原来的整场下载/主机计算；
- `write_fluid_flux`：GPU 路径调用 `cuda_backend_download_flux_planes`，只取 6 个截面平面；
- 新增 `write_fluid_probe_gpu`：批量计算 probe 索引/系数，调用 `cuda_backend_download_probes`；
- `computeFieldStat_blocks` / `calculate_turbulent_statistic_`：按需 `ensure macro`；
- `check_is_continue`：读取/排序/插值 continue 数据后，GPU 路径调用 `cuda_backend_upload_f` 重新上传 `fIn`；
- `cuda_boundary_supported`：判断 BC 是否全部在 GPU 支持列表内。

### 4.2 `LBMBlockComm.f90`

- 原 CPU 版本的 `IBM_FEM` 将流固相互作用力计算和固体求解放在同一子程序内；当前拆为 `IBM` 和 `FEM`，CPU 与 GPU 构建共用这一调度结构：
  - `IBM(carrierFluidId,time_IBM)`：调用 `FSInteraction_force`，完成罚力迭代、固体节点载荷与流体受力更新，并累计 IBM 耗时。CPU 分支使用主机数组；GPU 分支额外传入该流体块的 `cuda` 状态，使用设备端 IBM 后端；
  - `FEM(carrierFluidId,time_FEM,time)`：仍在 CPU 上执行固体求解。每次调用按 `flow%numsubstep` 推进固体子步，子步长为 `dh/numsubstep`，并累计 FEM 耗时；
- `tree_collision_streaming_IBM_FEM` 对每个块的调用顺序为：更新体积力参数 → 计算宏观量 → 清零流体受力 → `IBM` → 添加体积力 → 提取旧时刻通信层 → 碰撞 → 缓存半程反弹数据 → 传播 → `FEM` → 设置边界 → 提取新时刻通信层，然后递归推进子块并交换父子块数据。
- 相比 CPU 版本，FEM 从碰撞前移至传播后，CPU 路径也采用这个顺序。IBM 完成后，固体求解使用已得到的载荷，流体碰撞、传播使用当前流场和受力，不依赖此次 FEM 更新后的固体位置；因此当前依赖关系允许调整它们的先后顺序。CPU 顺序执行；GPU 路径先提交设备上的碰撞、传播任务，再由 CPU 求解固体，为两者重叠执行提供机会，实际重叠仍受同步依赖约束。
- `use cuda_lbm_backend`：引入父子块通信接口；
- `extract_interpolate_layer`：`#ifdef USE_CUDA` 下走 `extract_interpolate_layer_gpu`，父块采样/平均/写入子块面缓冲全部在设备上；
- 新增 `extract_interpolate_layer_gpu`：按子块面调用 `cuda_backend_extract_face` / `cuda_backend_average_face`；
- `interpolation_father_to_son`：GPU 分支走 `interpolation_father_to_son_gpu`（`cuda_backend_interpolate_face`）；
- `deliver_son_to_father`：GPU 分支走 `deliver_son_to_father_gpu`（`cuda_backend_record_event` + `cuda_backend_deliver_face`）；
- 块树结构、父子块递归顺序及细网格 subcycling 次数不变；块内 IBM/FEM 的拆分与调度调整如上所述。

### 4.3 `Solidbody.f90`

- `#ifdef USE_CUDA` 下 `use cudafor`，并把 `v_Ei/v_Ew/v_Ea/v_Evel/v_Eforce` 声明为 **pinned allocatable**；
- `calculate_interaction_force` / `FSInteraction_force_`：GPU 分支接收 `cuda_state`，调用 `cuda_backend_ibm_setup` / `cuda_backend_ibm_prepare`；`v_Eforce` 的清零只在 CPU 路径执行（GPU 在设备上清零）；
- `PenaltyForce_`：GPU 分支调用 `cuda_backend_ibm_penalty` 后直接返回；CPU 两遍循环保留；
- `PenaltyForce_` 的罚力迭代：GPU 与 CPU 一样每轮检查一次收敛（`IBM_CHECK_INTERVAL=1`）；
- `FluidVolumeForce_`：GPU 分支先调用 `cuda_backend_spread_ibm_force` 叠加 `force`，随后仍由主机循环计算梁载荷/力矩 `lodFlow`；
- `UpdateElmtInterp_` 保持不变（仍由主机计算 4×4×4 索引/权重）。

### 4.4 `main.f90`

- 初始化的 `calculate_macro_quantities_blocks()` 保留；
- CPU 路径（`#ifndef USE_CUDA`）在每步流场推进后、湍流统计前调用 `calculate_macro_quantities_blocks()`，保证统计和输出读取当前时刻的密度、速度；
- GPU 路径省略步末统一更新，流场输出等消费者通过 `ensure_macro_blocks` 或直接调用 `cuda_backend_macro_if_needed` 按需更新；需要主机数据时再下载；
- info 分支：GPU 路径调用 `write_fluid_probe_gpu`；CPU 路径调用 `write_fluid_information`，使用步末统一更新后的速度场，不再重复调用 `calculate_macro_quantities_iblock`。

### 4.5 `Makefile`

- `GPU=1` 为默认构建：`FC=nvfortran`，`FFLAGS` 增加 `-mp -cuda -gpu=cc80 -DUSE_CUDA`；
- `cuda_lbm.cuf` 单独编译为 `src/cuda_lbm.o`，并作为 `Solidbody.o`、`FluidDomain.o`、`LBMBlockComm.o` 的依赖，保证模块先编译；
- `make GPU=0` 可生成纯 CPU 版本（gfortran/ifort/ifx），用于对照和验证。

### 4.6 未改动的模块

`ConstParams.f90`、`FlowCondition.f90`、`SolidSolver.f90`、`Util.f90`、`forkthread.c`、`velocitymap.cpp` 中没有 GPU 专属代码；它们与所提供 CPU master 的差异来自版本差异，而不是 CUDA 移植改动。GPU 版本直接复用了这些模块的 CPU 实现（固体求解、输入解析、插值工具、fork 输出等）。
