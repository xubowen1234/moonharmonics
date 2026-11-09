# MoonHarmonics 选题与严格语义查重

检查日期：2026-11-09，Asia/Shanghai。用户授权新项目自动选题与完整交付。
固定比赛要求沿用本地 skill，不访问在线章程。此前约定的 2026-10-30 24:00 已过去，本文不声称延期或已报名。

## 方法与范围

执行已安装 osc2026-guide 的 Project Research Guide：先 `moon search <词> --limit 100`，再读取 MoonCakes 文档和所链接仓库的 README、API、源码、测试。
使用 MoonBit 0.10.18 / moon 0.1.20261028。通过项目绑定包装器读取 GitHub 仓库 API、树及源文件，搜索 `language:MoonBit`，每页上限100，保存 total_count/incomplete_results。
研究原始输出仅保存在本地验收目录，不包含凭据；下面列出可公开复现的范围。

四个正式候选涉及网页归档、数字波形、数据序列化、球面数值分析四个领域。每个候选均查 MoonCakes、GitHub 和永久注册表，避免按名称推断。
扩词包括 WARC/web archive/warcio/crawl/archive；VCD/waveform/verilog；Avro/schema/serialization；spherical harmonics/real spherical/associated legendre/gauss legendre/spherical transform/ambisonic/ambisonics/SH transform/spectral poisson/real_basis/sphere projection。
另探索 NRRD、NIfTI、raster volume、medical imaging、metaimage；因相邻体数据工具较多，未选入最终候选。

`archive`、`serialization`、`http`、`digest` 返回达到100的上限，不能视为完整列表。其余球谐相关词结果通常0—24条，按实际行为阅读，不以零结果证明唯一。
GitHub 仓库检索 spherical harmonics、real spherical、ambisonic、spherical spectrum、NRRD、WARC、Avro、VCD；文本元数据搜索可能漏掉藏在多功能库内的实现，因此补查相邻库完整 API、目录和源码。
官方 awesome 清单 README 地址不可用，记录为限制；MoonCakes 以及 moonbit-community 相关项目通过 CLI 搜索覆盖，不能声称穷尽所有官方清单。
NIST DLMF 在线正文访问失败；数学约定以可访问的 SciPy/Pyshtools 材料与独立解析公式交叉核验。

## 需求落地与完整工作流

| 候选 | 真实用户及输入→处理→输出 | 成熟依据 | 可完成核心与验收 |
|---|---|---|---|
| 网页归档 WARC | 爬虫保存 HTTP 记录；档案机构检查帧；研究者建立索引 | warcio / ISO28500 | 字节安全分帧、读写、字段约束、截断测试；但存在直接替代 |
| 数字波形 VCD | RTL 仿真读信号；FPGA工具追踪事件；测试系统分析时序 | IEEE VCD / GTKWave | 头部/值变更解码和事件查询；但存在直接替代 |
| Apache Avro | 消息系统编码；ETL读对象容器；schema工具检查演进 | Apache Avro | 模式、datum、容器的独立对照；但存在直接替代 |
| 通用实球谐分析 | 图形引擎输入方向光照→加权投影/余弦卷积→可重建漫反射方向场；空间音频工具输入方向样本→实基编码/z轴旋转→方向系数与重建值；科学程序输入球面标量场→谱分析/热滤波/Poisson求解→平滑场或势函数 | SHTOOLS/pyshtools BSD-3-Clause；SciPy BSD-3-Clause；正交球谐与Gauss–Legendre公式 | 有限阶正逆变换、谱模型、拟合及验证，解析低阶值、SciPy数值对照、正交性、往返、Parseval和算子残差 |

没有用户访谈证据；三类任务来自成熟跨生态库的使用模型，并非虚构参赛者经验。复用接口是 MoonBit 公共数值 API，不绑定某个图形/音频/地理引擎。

## 相邻项目行为矩阵

| 项目及检查入口 | 实际能力及源码证据 | 重合判断与决定 |
|---|---|---|
| [moon-warc](https://github.com/xiaojing012/moon-warc) 0.1.0；[MoonCakes](https://mooncakes.io/docs/xiaojing012/moon-warc) | README、decoder.mbt、decoder_test.mbt 已读；WARC解析、写入、验证、索引 | 直接重合，淘汰 WARC；某次pkg.generated.mbti路径404，已有核心源码/测试足以判定 |
| [moonbit-vcd](https://github.com/lqhy1234567/moonbit-vcd) 0.3.1；[MoonCakes](https://mooncakes.io/docs/lqhy1234567/moonbit-vcd) | README、完整公开接口、parser_test.mbt；含查询、断言、协议解码 | 直接重合，淘汰 VCD；MiniVCDMBT进一步佐证 |
| [moon-avro](https://github.com/yugonlian/moon-avro) 0.3.0；[MoonCakes](https://mooncakes.io/docs/yugonlian/moon-avro) | README、根接口、codec源码/测试、容器目录 | 直接重合，淘汰 Avro |
| [moonbit-geofield](https://github.com/2515050243qyf/moonbit-geofield) 0.2.0；[MoonCakes](https://mooncakes.io/docs/2515050243qyf/moonbit-geofield) | harmonics.mbt、完整接口、advanced_test.mbt、README；带地球半径/海拔径向项的有限地理模型 evaluate_spherical_harmonic | 高度相邻。未见通用采样正逆变换与球面谱算子；本库用无量纲单位球、明确归一化、网格/系数模型供地理及非地理系统复用。放进该地理库会引入不必要地球模型约束，独立底层模块更合适 |
| [symbit](https://github.com/CAIMEOX/symbit) 0.5.10 | spherical_harmonics.mbt、polynomials.mbt、tests、README；ynm/znm仅构造符号树，并明确未实现归一化/显式展开 | 同方向不同功能，不是数值采样变换 |
| [moon-geometry](https://github.com/walkzzz/moon-geometry) 0.251.0 | README、完整接口、bounding_sphere源码/测试和完整树 | 几何包围球/GIS变换，不是球面谱投影；可消费本库 |
| [mbmath](https://github.com/CAIMEOX/mbmath) / CAIMEOX/moon_floating 0.3.1 | README、目录、quadrature_test.mbt；高精度特殊函数/求积，搜索的harmonic主要为调和数 | 通用数值底座，未见球面场变换；初按模块名构造仓库路径404，随后使用真实链接重新检查 |
| [moonbit-dsp](https://github.com/hwlxmm/moonbit-dsp) 0.1.0 | README、目录、spatial_audio公开接口/测试；FFT/DCT、双耳panning，接口只有compute_binaural_panning | 时间/平面谱及声像增益不同于球谐系数分析，可集成使用 |
| [moontide](https://github.com/WYT-ai-nb/moontide) 0.3.0 | README及MoonCakes实际说明，API文件路径404未声称读到 | NOAA时间分潮调和预测，不是二维球面正交变换 |
| [moonnifti](https://github.com/wangjiale6036-dotcom/moonnifti) 0.2.0 | README、接口、目录；体素/仿射世界坐标变换 | NRRD探索的强近邻，未选择另一体数据工具 |

## 永久注册表与保留边界

读取永久注册表项目标题清单，按 harmonic/spherical/legendre/raster/medical/volume/transform 及任务检索，并阅读匹配完整指纹和 MoonOT、MoonPhone、MoonCite、MoonRecoder、MoonConstruct、MoonIon、MoonWell、MoonARFF 的完整相关记录。
MoonBench性能门禁、MoonContract HTTP契约、MoonRecur民用日历递推、MoonShard内容分块、MoonPhone号码、MoonCite引用排版均视为保留，不复用其核心实现。
MoonOT的最优传输耦合、MoonWFST的自动机、MoonDiffraction衍射、MoonLoudness时间音频响度与球面实正交展开的数据模型和主流程不同；不以换名规避保留题目。
注册表中其他 reserved/completed/rejected 指纹同样保留；本项目不复用旧仓库或历史。

## 排序与边界

评分顺序：广度、用户多样性、复用、低重合风险、深度、可验证性、可行性、维护性（各1–5）。
MoonHarmonics：5/5/5/4/4/5/4/4=36，需求与查重门禁通过。
WARC：5/5/5/1/4/5/4/4=33；VCD：4/4/5/1/4/5/4/4=31；Avro：5/5/5/1/5/5/3/3=32；三者直接重合为硬淘汰，不能由总分抵消。
选择最高分且通过门禁的 MoonHarmonics，2026-11-09 用注册表同目录排他锁重新检查并原子追加预留及 external-checked，没有覆盖旧记录。

交付核心：l≤64、实正交归一化、Condon–Shortley相位、弧度、Gauss–Legendre×均匀经度网格、可复现数值分析。非目标：HEALPix、自旋加权球谐、任意SO(3)旋转、GPU、完整SHTOOLS能力、真实音频设备或地球重力模型。
验收在后续质量文档记录实测值，不把预定目标写成已完成结果。

结论：**在已记录检索范围内未发现成熟直接替代。** 未索引、私有、离线、不同标签或检索上限外项目不能排除；发布前再反向复核核心行为。
