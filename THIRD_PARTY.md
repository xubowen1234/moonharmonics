# 参考、依赖与许可证

MoonHarmonics 的生产实现是依据数学公式独立编写的 MoonBit 代码，MIT 许可证，Copyright (c) 2026 xubowen1234。没有复制 SciPy、SHTOOLS 或相邻 MoonBit 项目的源码。归属不表示上述项目为本库背书。

| 项目或资料 | 许可证/用途 | 本仓库边界 |
|---|---|---|
| moonbitlang/core，https://github.com/moonbitlang/core | Apache-2.0；数学、JSON、基础容器与调试接口 | 工具链提供的运行依赖，未 vendoring 源码 |
| SciPy，https://github.com/scipy/scipy | BSD-3-Clause；开发参考 `sph_harm_y`、`roots_legendre` 与 SVD `lstsq` | CI 按精确版本安装；本库无 Python 运行时依赖，未复制实现 |
| SHTOOLS，https://github.com/SHTOOLS/SHTOOLS | BSD-3-Clause；成熟球面谱处理能力对照 | 阅读 README、许可证以界定用途与边界，未复制代码，不兼容其所有约定/API |
| NIST DLMF 14.30，https://dlmf.nist.gov/14.30 | 球谐数学定义的参考入口 | 本轮抓取正文失败，不能把链接当成已读取全文的证据；公式由解析测试及独立数值参考验证 |

Gauss–Legendre、归一化关联 Legendre 递推、Householder QR、Laplace–Beltrami 特征值与球面余弦卷积使用通用数学算法。测试中的解析数值和实时计算的 SciPy 对照结果不构成对上游实现代码的复制。没有来源不明的数据包或训练数据随发布物分发。
