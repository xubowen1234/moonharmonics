# MoonHarmonics

MoonHarmonics 是 MoonBit 的实球谐分析与谱算子库，用于把单位球面上的标量场转换为可复用的有限阶系数，并重建、平滑、卷积、旋转或求解球面 Poisson 方程。核心使用 MoonBit 和 `moonbitlang/core`，运行时不调用 Python、SciPy 或外部 C 数值库。

图形程序可以处理方向辐亮度及漫反射辐照度；空间音频离线工具可以编码方向、转换系数约定和旋转方向场；科学数据处理程序可以分析球面采样、检查谱能量并求解零均值势场。这三类完整工作流见下方可运行示例。

## 环境与安装

本项目使用 `moonc 0.10.22+124d61d05`，高于最低的 0.10.16。CI 固定同一官方工具链，覆盖 `wasm-gc`、`wasm`、`js`、Linux `native`。JS 运行需要 Node.js；Windows native 需要 MSVC C++ 工具和 Windows SDK。本机缺少该工具链，不能把 Windows native 构建计为本地通过。

MoonCakes 当前已发布 `0.1.0`；本文所述 0.2.0 修订（含 DH1）已在本地完成实现和验证，尚未发布。发布版本可安装：

```sh
moon add xubowen1234/moonharmonics@0.1.0
```

在消费包的 `moon.pkg` 添加：

```moonbit
import {
  "xubowen1234/moonharmonics" @h,
}
options(is_main: true)
```

写入 `main.mbt`，在消费模块根执行 `moon run . --target wasm-gc`：

```moonbit
fn main {
  let grid = @h.SphereGrid::new(4).unwrap()
  let spectrum = grid.project(fn(d) {
    let (_, _, z) = d.cartesian()
    2.0 + z
  }).unwrap()
  let value = spectrum.evaluate(@h.Direction::new(0.0, 0.0).unwrap()).unwrap()
  println(value) // 约为 3
}
```

`unwrap()` 仅用于参数已知正确的示例。处理外部输入时，应匹配 `Ok` / `Err` 并处理 `HarmonicError`。

## 已实现的能力

- 0–64 阶实球谐、归一化关联 Legendre 稳定递推及明确的有符号阶布局；64 阶共 4225 个系数。
- 经过验证的球面角度及笛卡尔方向，不可变 `Spectrum`，显式补零、截断和算术运算。
- Gauss–Legendre × 等间隔经度网格、可分离正反变换、过采样及流式加权投影。
- Driscoll–Healy DH1 等角网格正反变换；纬度与经度各 2(L+1) 点，带限阶数 0–64。
- 分阶功率、均值、RMS、Parseval 缺陷、L2 距离和梯度能量。
- 热方程平滑、带通及余弦渐消、Laplacian、零均值 Poisson 和屏蔽 Poisson 求解。
- 任意 SO(3) 主动旋转（ZYZ 欧拉角）、方向编码、纬圈对称核卷积、包含高偶数阶的 Lambertian 核。
- 有极点极限的笛卡尔切向梯度、精确经度导数。
- 不规则样本的带列主元 Householder QR 加权拟合，严格版本化 JSON，以及三种归一化、两种相位、三种排列的 18 种组合双向转换（含 SN3D、N3D/ACN）。

## 数据模型与数学约定

角度单位是弧度。`theta` 为从正 z 轴开始的余纬，范围 `[0, pi]`；`phi` 为经度，归一化到 `[0, 2*pi)`。经度及旋转输入绝对值最多 `1e12`，避免无限大角度约化。向量构造先缩放再归一化，拒绝零向量与非有限输入。

系数使用 `index = l*l + l + m`，`-l <= m <= l`，包括从 0 到声明最大阶的所有带。实基正交归一，包含 Condon–Shortley 相位：`m>0` 为余弦项，`m<0` 为正的正阶复球谐虚部对应的正弦项。因此一阶基按索引排列为 `sqrt(3/(4*pi)) * [-y,z,-x]`。不要直接把其他音频/地学约定的数组当作本库系数。

最小网格有 `L+1` 个按 `cos(theta)` 升序排列的纬圈及 `2*L+1` 个经度，平面索引 `ring*nphi+longitude`。球面积分权重已包含 `2*pi/nphi`，不可再乘 `sin(theta)`。在精确算术下，此网格对不超过 L 阶的带限场可精确恢复；一般非带限输入存在混叠，应使用 `SphereGrid::oversampled` 并做收敛检查。

`Spectrum::rotate_z(alpha)` 定义主动作用 `g(theta,phi)=f(theta,phi-alpha)`。`Spectrum::rotate_zyz(alpha,beta,gamma)` 定义右手主动旋转 `R=Rz(alpha)*Ry(beta)*Rz(gamma)`，最右侧先作用，旋转后 `g(d)=f(R^-1*d)`。方向向量可调用同名方法；逆旋转参数为 `(-gamma,-beta,-alpha)`。`laplacian()` 的特征值是 `-l*(l+1)`；`poisson(tolerance)` 解 `Delta(u)=f` 并固定 `u00=0`。容差是源系数 `a00` 的绝对容差，超过容差会拒绝，不会静默移除均值；需要移除时调用 `without_mean()`。

## API 用法

```moonbit
let grid = @h.SphereGrid::new(8).unwrap()
let field = grid.project(fn(d) {
  let (_, _, z) = d.cartesian()
  z
}).unwrap().without_mean()
let smooth = field.heat(0.2).unwrap()
let potential = field.poisson(0.0).unwrap()
let turned = field.rotate_zyz(0.4, 0.8, -1.2).unwrap()
let irradiance = field.convolve(@h.ZonalKernel::lambertian(8).unwrap()).unwrap()
let persisted = @h.Spectrum::from_json_string(field.to_json_string()).unwrap()
let dh = @h.DriscollHealyGrid::new(8).unwrap()
let equiangular_coefficients = dh.project(fn(d) {
  let (_, _, z) = d.cartesian()
  z
}).unwrap()
let dh_samples = equiangular_coefficients.synthesize_dh(dh).unwrap()
```

`DriscollHealyGrid` 使用 DH1 等间隔网格，节点为 `theta_j=pi*j/N`、`phi_k=2*pi*k/N`，其中 `N=2(L+1)`；按纬圈优先存储，包含北极、不重复存储南极。分析使用 Driscoll–Healy 纬向求积权重，不把等角采样误当作等面积采样。

### 与 SHTOOLS 的能力对照

SHTOOLS 提供 DH 等角与 Gauss–Legendre 网格的球谐展开/重建，也支持复球谐和更广泛的地学分析。本库当前本地修订实现实标量 DH1（对应 SHTOOLS `sampling=1` 的等采样网格）与 Gauss–Legendre 网格正反变换，最高 64 阶；经度方向采用直接三角求和，不使用 SHTOOLS 的 FFT 路径，也不提供复球谐、`sampling=2` 网格、SHTOOLS 文件适配器或其全部分析功能。网格采样布局虽可对应，系数约定仍需显式对齐。详细边界见仓库 `docs/shtools-comparison.md`。

`Spectrum::new` 和 `values()` 分别复制输入、输出，调用者不能通过数组修改已验证的系数。`filter_bands` 和卷积要求完整相同阶数，不隐式截断。发生非有限输入、形状错误、数值溢出或拟合秩不足时返回错误。

`fit_weighted(L, directions, values, weights)` 求最小化 `sum w_i*(f_i-y_i)^2`。支持 `L<=12`、最多 4096 行，权重有限且非负。QR 对秩不足返回 `SingularSystem`，不返回伪逆最小范数解。`diagonal_ratio()` 是 QR 对角线比值诊断，不是严格条件数。用于统计拟合的权重和流式投影的立体角权重具有不同含义。

`from_convention` / `to_convention` 显式选择 `Orthonormal`、`FourPi`、`Schmidt`，`CondonShortley` / `PhaseFree`，以及 `Signed`、`Interleaved`、`CosineThenSine` 排列。`FourPi + PhaseFree + Signed` 对应 N3D/ACN，`Schmidt + PhaseFree + Signed` 对应 SN3D/ACN。旧接口 `to_sn3d()` 返回无 CS 相位的 SN3D **展开系数**，不是基函数采样值；打包次序保持不变。JSON 固定格式 `moonharmonics-real-v1`，携带 `orthonormal-cs-l*l+l+m` 约定，拒绝重复键、未知字段及错误形状。格式和旋转资料见仓库 `docs/`。

## 三个可运行工作流

从本仓库根目录运行：

```sh
moon run examples/graphics --target wasm-gc
moon run examples/audio --target wasm-gc
moon run examples/science --target wasm-gc
```

1. 图形示例从 `2+z` 辐亮度采样生成系数，卷积半球余弦核，验证北极辐照度 `8*pi/3`，再乘反照率除以 pi 得到反射辐亮度。
2. 音频工具示例编码两个带权方向，从 SN3D/ACN 导入，做三个角均明确的任意空间旋转，验证旋转前后对应方向值一致，导出无 CS 相位 N3D 并改变排列，最后重新导入核验。它是角向标量处理示例，不提供实时设备播放或 HRTF。
3. 科学示例处理 `z+x*y` 场，求解 `u=-z/2-x*y/6`，检查 Laplacian 残差，并验证热平滑使能量下降。

## 检查、测试与独立参考

```sh
moon info
moon fmt --check
moon check --target wasm-gc --deny-warn
moon build --target wasm-gc --deny-warn
moon test --target wasm-gc --deny-warn
```

将目标替换为 `wasm`、`js` 或 `native` 可复核其他后端。当前 31 个测试用例。测试包括全部 4225 模式索引、解析一阶/二阶基、64 阶加法定理及反演、Gauss 求积多项式精度、DH1 面积与 0/1/4/12/32/64 阶正反变换、Parseval、独立解析场、旋转相位、SO(3) 逆变换/组合/逐带能量及 64 阶稳定性、18 种约定往返和独立带标记、Poisson 符号、梯度极点极限、病态拟合拒绝、JSON 歧义、溢出和错误输入。测试用例内部遍历大量模式，测试数量不等于数值比较数量。

开发验证脚本另外使用 SciPy 作为独立参考（不是运行时依赖）：

```sh
python -m pip install numpy==2.3.4 scipy==1.16.3
python scripts/reference_check.py --target wasm-gc
```

它独立比较 64 阶实基、256 节点求积、16 阶 Gauss–Legendre 与 DH1 正反变换、解析导数、64 阶 Lambertian 系数和带噪声加权拟合，并用独立三维旋转后的 SciPy 球谐验证 64 阶旋转，使用稠密 Padé 矩阵指数对照最高阶系数，逐项验证 18 种约定的导入与导出。每个后端执行 112,152 次标量比较。拟合参考使用 SVD，区别于本库 QR。脚本打印实际比较数量及最大绝对误差；阈值和输入在脚本中固定，不能用一次往返测试替代独立参考。详细输入与实测误差见仓库 `docs/numerical-validation.md`。

## 边界、复杂度与许可

支持单位球面上的实标量有限阶场，使用 Double。网格变换使用可分离直接傅里叶求和，不是 FFT；对于最小网格，主要工作量为 `O(L^3)`，存储约为 `O(L^2)`。任意点基求值为 `O(L^2)`；密集拟合矩阵为 `O(N*(L+1)^2)` 存储，QR 工作量约 `O(N*(L+1)^4)`。输入非常大或条件很差时会报错，库不承诺任意尺度下的高相对精度。

SO(3) 旋转按阶应用角动量反对称生成元的指数，采用每步范数小于 1/2 的 18 项 Taylor 作用；不重采样，工作量 `O(L^3*|beta|)`，约化后 `|beta|<=pi`，总存储 `O(L^2)`。它面向离线数值处理，不承诺实时音频延迟。

不提供 HEALPix、spin-weighted 球谐、FuMa 全格式、未归一化阶乘型系数、GPU、实时音频引擎、地球重力径向模型或 SHTOOLS 全功能兼容。Lambertian 卷积包括所有支持的偶数阶；低阶截断仍可能对尖锐照明产生振铃，不保证重建场非负。

本项目为公式驱动的独立 MoonBit 实现，采用 MIT 许可证。数学和验证参考包括 SciPy（BSD-3-Clause）及 SHTOOLS（BSD-3-Clause）；未复制其实现代码，也不声称完整移植。运行时依赖 `moonbitlang/core`（Apache-2.0）。归属与来源边界见 `THIRD_PARTY.md`，AI 辅助开发说明见 `AI_USAGE.md`，版本变更见 `CHANGELOG.md`。

[GitHub](https://github.com/xubowen1234/moonharmonics) · [MoonCakes 模块](https://mooncakes.io/docs/xubowen1234/moonharmonics) · [持续集成](https://github.com/xubowen1234/moonharmonics/actions)
