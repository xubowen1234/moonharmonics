///|
/// Cartesian tangent gradient on the unit sphere; unlike coordinate derivatives it is defined at poles.
pub struct TangentGradient {
  x : Double
  y : Double
  z : Double
}

///|
pub fn TangentGradient::components(
  self : TangentGradient,
) -> (Double, Double, Double) {
  (self.x, self.y, self.z)
}

///|
pub fn TangentGradient::norm(self : TangentGradient) -> Double {
  let scale = @cmp.maximum(
    self.x.abs(),
    @cmp.maximum(self.y.abs(), self.z.abs()),
  )
  if scale == 0.0 {
    return 0.0
  }
  let x = self.x / scale
  let y = self.y / scale
  let z = self.z / scale
  scale * (x * x + y * y + z * z).sqrt()
}

///|
/// Evaluate the gradient using differentiated recurrences, including exact polar limits.
/// Recurrences differentiate with respect to theta directly, avoiding division by sin(theta).
pub fn Spectrum::gradient(
  self : Spectrum,
  direction : Direction,
) -> Result[TangentGradient, HarmonicError] {
  let x = @math.cos(direction.theta)
  let s = @math.sin(direction.theta)
  if direction.theta == 0.0 || direction.theta == @math.PI {
    let mut gx = 0.0
    let mut gy = 0.0
    for l = 1; l <= self.degree; l = l + 1 {
      let mut factor = -((2 * l + 1).to_double() *
      (l * (l + 1)).to_double() /
      (8.0 * @math.PI)).sqrt()
      if direction.theta == @math.PI && l % 2 == 0 {
        factor = -factor
      }
      gx = gx + factor * self.coefficients[index_of(l, 1)]
      gy = gy + factor * self.coefficients[index_of(l, -1)]
    }
    if !finite(gx) || !finite(gy) {
      return Err(InvalidInput("polar gradient overflow"))
    }
    return Ok({ x: gx, y: gy, z: 0.0, })
  }
  let q = legendre_factors(self.degree, x, s)
  let derivative = Array::make(self.size(), 0.0)
  for m = 0; m <= self.degree; m = m + 1 {
    if m > 0 {
      let md = m.to_double()
      let factor = -((2.0 * md + 1.0) / (2.0 * md)).sqrt()
      derivative[index_of(m, m)] = factor *
        (x * q[index_of(m - 1, m - 1)] + s * derivative[index_of(m - 1, m - 1)])
    }
    if m < self.degree {
      let factor = (2.0 * m.to_double() + 3.0).sqrt()
      derivative[index_of(m + 1, m)] = factor *
        (-s * q[index_of(m, m)] + x * derivative[index_of(m, m)])
    }
    for l = m + 2; l <= self.degree; l = l + 1 {
      let ld = l.to_double()
      let md = m.to_double()
      let a = ((4.0 * ld * ld - 1.0) / (ld * ld - md * md)).sqrt()
      let b = ((2.0 * ld + 1.0) *
      ((ld - 1.0) * (ld - 1.0) - md * md) /
      ((2.0 * ld - 3.0) * (ld * ld - md * md))).sqrt()
      derivative[index_of(l, m)] = a *
        (-s * q[index_of(l - 1, m)] + x * derivative[index_of(l - 1, m)]) -
        b * derivative[index_of(l - 2, m)]
    }
  }
  let mut dt = 0.0
  let mut dp = 0.0
  for l = 0; l <= self.degree; l = l + 1 {
    dt = dt + self.coefficients[index_of(l, 0)] * derivative[index_of(l, 0)]
    for m = 1; m <= l; m = m + 1 {
      let angle = m.to_double() * direction.phi
      let c = @math.cos(angle)
      let sn = @math.sin(angle)
      let pos = self.coefficients[index_of(l, m)]
      let neg = self.coefficients[index_of(l, -m)]
      dt = dt + 2.0.sqrt() * derivative[index_of(l, m)] * (pos * c + neg * sn)
      dp = dp +
        2.0.sqrt() *
        (q[index_of(l, m)] / s) *
        m.to_double() *
        (-pos * sn + neg * c)
    }
  }
  let cphi = @math.cos(direction.phi)
  let sphi = @math.sin(direction.phi)
  let gx = x * cphi * dt - sphi * dp
  let gy = x * sphi * dt + cphi * dp
  let gz = -s * dt
  if !finite(gx) || !finite(gy) || !finite(gz) {
    Err(InvalidInput("gradient overflow"))
  } else {
    Ok({ x: gx, y: gy, z: gz, })
  }
}
