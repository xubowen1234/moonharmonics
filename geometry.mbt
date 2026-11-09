///|
/// A point on the unit sphere: colatitude theta in [0,pi], longitude phi in [0,2pi).
pub struct Direction {
  theta : Double
  phi : Double
}

///|
fn finite(x : Double) -> Bool {
  x == x && x - x == 0.0
}

///|
/// Validate angles in radians. Longitudes up to 1.0e12 radians are reduced modulo 2pi.
pub fn Direction::new(
  theta : Double,
  phi : Double,
) -> Result[Direction, HarmonicError] {
  if !finite(theta) || theta < 0.0 || theta > @math.PI {
    return Err(InvalidInput("colatitude must be finite and in [0,pi]"))
  }
  if !finite(phi) || phi.abs() > 1.0e12 {
    return Err(InvalidInput("longitude must be finite and bounded by 1.0e12"))
  }
  let period = 2.0 * @math.PI
  let reduced = phi - (phi / period).floor() * period
  let canonical = if reduced >= period {
    0.0
  } else if reduced < 0.0 {
    reduced + period
  } else {
    reduced
  }
  Ok({ theta, phi: canonical, })
}

///|
pub fn Direction::theta(self : Direction) -> Double {
  self.theta
}

///|
pub fn Direction::phi(self : Direction) -> Double {
  self.phi
}

///|
/// Return the Cartesian unit vector, with theta measured from the positive z axis.
pub fn Direction::cartesian(self : Direction) -> (Double, Double, Double) {
  let s = @math.sin(self.theta)
  (s * @math.cos(self.phi), s * @math.sin(self.phi), @math.cos(self.theta))
}

///|
/// Normalize a nonzero Cartesian vector. Scale first to avoid squared-length overflow.
pub fn Direction::from_cartesian(
  x : Double,
  y : Double,
  z : Double,
) -> Result[Direction, HarmonicError] {
  if !finite(x) || !finite(y) || !finite(z) {
    return Err(InvalidInput("Cartesian coordinates must be finite"))
  }
  let scale = @cmp.maximum(x.abs(), @cmp.maximum(y.abs(), z.abs()))
  if scale == 0.0 {
    return Err(InvalidInput("zero vector has no direction"))
  }
  let a = x / scale
  let b = y / scale
  let c = z / scale
  let norm = (a * a + b * b + c * c).sqrt()
  let cos_theta = @cmp.maximum(-1.0, @cmp.minimum(1.0, c / norm))
  Direction::new(@math.acos(cos_theta), @math.atan2(b, a))
}

///|
/// Numerically stable great-circle distance in radians, including antipodes.
pub fn Direction::distance(self : Direction, other : Direction) -> Double {
  let (a, b, c) = self.cartesian()
  let (x, y, z) = other.cartesian()
  let cross_x = b * z - c * y
  let cross_y = c * x - a * z
  let cross_z = a * y - b * x
  let sine = (cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).sqrt()
  @math.atan2(sine, a * x + b * y + c * z)
}
