# 与 SHTOOLS 的能力边界

核对日期：2026-11-17。对照对象为 SHTOOLS 公开文档中的球谐变换和球面网格接口，不是逐函数兼容声明。

| 能力 | MoonHarmonics | SHTOOLS / pyshtools |
|---|---|---|
| 带限实标量正反变换 | 支持，0–64 阶 | 支持 |
| Gauss–Legendre 求积网格 | 支持可分离分析/合成 | 支持 GLQ 展开与合成 |
| Driscoll–Healy 等角网格 | 支持 DH1：每轴 2(L+1) 点，对应 SHTOOLS `sampling=1` 的等采样尺寸 | 支持 `sampling=1`（N×N）与 `sampling=2`（N×2N） |
| 经度变换 | 直接三角求和，不是 FFT | DH 变换使用 FFT 型经度求和路径 |
| 系数模型 | 实数、正交归一、CS 相位，`l*l+l+m` | 有多种归一化、相位和系数存储选项 |
| 转换 | 归一化、相位、带内排列 18 种组合双向转换 | 更广泛的地学与语言 API |
| 其他球谐类别 | 不支持复球谐、自旋加权球谐 | 有复数变换、多窗谱分析及地球/磁场分析等额外接口 |

MoonHarmonics 的 DH1 采样点数、纬度间距、经度间距及省略南极行的布局对应 SHTOOLS `sampling=1`；不支持其每纬圈经度数加倍的 `sampling=2` 布局。即使采样网格可对应，交换系数仍须对齐归一化、相位和带内排列。当前库不提供 SHTOOLS 文件读写器或完整 API 适配层。

本项目采用 MoonBit 的可移植实现与通用直接求和，最高 64 阶。此范围不与 SHTOOLS 的性能或规模竞争；选择 DH 网格的主要价值是接受常见的等角采样，而 GL 网格对相同带限模型提供较少纬圈的求积采样。对非带限输入，两类离散投影都可能发生混叠，应提高采样密度并检查收敛。

参考文档：

- [SHTOOLS 球谐变换接口](https://shtools.github.io/SHTOOLS/spherical-harmonic-transforms.html)
- [SHTOOLS 网格格式与 DH1/DH2 采样](https://shtools.github.io/SHTOOLS/fortran-grid-formats.html)
- [SHTOOLS 实现细节（经度 FFT 与纬向求积）](https://shtools.github.io/SHTOOLS/fortran-implementation-details.html)
- [SHTOOLS 项目能力与许可证](https://shtools.github.io/SHTOOLS/shtools.html)
- [SHTOOLS 项目与许可证](https://github.com/SHTOOLS/SHTOOLS)
