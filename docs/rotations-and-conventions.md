# 空间旋转与系数互操作

## 主动 SO(3) 旋转

`Spectrum::rotate_zyz(alpha,beta,gamma)` 作用矩阵为 `Rz(alpha) Ry(beta) Rz(gamma)`。采用右手坐标系，参数为弧度，矩阵从右往左作用。新的场满足 `g(d)=f(R^-1 d)`。这包括任意三维正旋转，不包括反射。

```moonbit
let rotated = field.rotate_zyz(0.4, 0.8, -1.2).unwrap()
let target = direction.rotate_zyz(0.4, 0.8, -1.2).unwrap()
// rotated.evaluate(target) 约等于 field.evaluate(direction)
let restored = rotated.rotate_zyz(1.2, -0.8, -0.4).unwrap()
```

`rotate_z` 保留直接的正弦/余弦二元操作。`rotate_y` 在每个 l 带中把实系数转换为复展开系数：

- `a_0=c_0`；
- `a_m=(c_m-i*c_-m)/sqrt(2)`；
- `a_-m=(-1)^m*conj(a_m)`。

y 旋转的生成元为实反对称三对角矩阵，`G[m,m+1]=sqrt((l-m)(l+m+1))/2`，下对角为其负转置。计算 `exp(beta*G) a`，每步 `h=beta/ceil(2*(l+1)*|beta|)`，因此 `||hG||_infinity<1/2`。每步累加到 18 阶 Taylor 项；每带先按系数最大绝对值缩放，避免不必要的中间溢出。这是有限精度数值旋转，未使用采样投影或不稳定的大阶乘相消公式。

不改变最大阶数，不混合不同 l 带；精确算术下逐带能量不变。64 阶测试检查逆变换、组合、任意点求值与能量；独立参考另外用 SciPy 空间旋转后的球谐及稠密 Padé 指数验证。工作量为 `O(L^3*|beta|)`，约化后 `|beta|<=pi`，包含输出在内存储为 `O(L^2)`。不提供实时处理时限保证。

所有角必须有限且绝对值不超过 `1e12`，超范围返回错误。极大合法角在浮点约化中仍可能损失精度；高精度应用应先提供接近一个周期的角度。

## 三种归一化、两种相位、三种排列

`from_convention` 导入外部展开系数；`to_convention` 导出。两个函数均显式指定以下三个参数，不猜测来源。

| 归一化 | 外部基相对正交归一实基 R_lm 的因子 k_l | 单个基的球面积分平方 |
|---|---|---|
| `Orthonormal` | 1 | 1 |
| `FourPi` | sqrt(4*pi) | 4*pi |
| `Schmidt` | sqrt(4*pi/(2*l+1)) | 4*pi/(2*l+1) |

`CondonShortley` 保留本库的 CS 相位；`PhaseFree` 相对本库基再乘 `s_m=(-1)^|m|`。所以外部基 `B_lm=k_l*s_m*R_lm`，外部**展开系数**导入时乘 `k_l*s_m`，导出时除以该因子；不要将这个方向直接用于基函数采样值。

`FourPi + PhaseFree + Signed` 表示 N3D/ACN，`Schmidt + PhaseFree + Signed` 表示 SN3D/ACN。`from_sn3d` / `to_sn3d` 是后一组合的便利接口。其他组合须把三项参数连同数组一起保存在调用方的数据协议中。

每带从 `l*l` 开始，各排列中的 l=2 示例为：

| 排列 | m 顺序 | 对规范带 `[a_-2,a_-1,a_0,a_1,a_2]` 的导出结果（不改变相位及归一化） |
|---|---|---|
| `Signed` | -2,-1,0,+1,+2 | `[a_-2,a_-1,a_0,a_1,a_2]` |
| `Interleaved` | 0,+1,-1,+2,-2 | `[a_0,a_1,a_-1,a_2,a_-2]` |
| `CosineThenSine` | 0,+1,+2,-1,-2 | `[a_0,a_1,a_2,a_-1,a_-2]` |

```moonbit
let field = @h.Spectrum::from_convention(
  1, [1.0, 2.0, 3.0, 4.0], @h.Schmidt, @h.PhaseFree, @h.Signed,
).unwrap()
let n3d_acn = field.to_convention(@h.FourPi, @h.PhaseFree, @h.Signed).unwrap()
// n3d_acn 约为 [1,2/sqrt(3),3/sqrt(3),4/sqrt(3)]
```

支持 0..64 阶完整实系数带，长度必须恰为 `(L+1)^2`。形状错误、非有限值和转换溢出会返回错误；构造和访问均防御复制。未支持 FuMa 的整体通道/权重规则、任意厂家声道标记、未归一化阶乘型基或复球谐输入。JSON v1 仍只存规范内部约定，未改变格式或自动接受外部标签。
