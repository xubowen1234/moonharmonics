///|
/// Active rotation about +z: rotated(theta,phi)=original(theta,phi-angle).
/// This is an exact two-by-two rotation in each signed-order pair, not general SO(3).
pub fn Spectrum::rotate_z(
  self : Spectrum,
  angle : Double,
) -> Result[Spectrum, HarmonicError] {
  if !finite(angle) || angle.abs() > 1.0e12 {
    return Err(
      InvalidInput("rotation angle must be finite and bounded by 1.0e12"),
    )
  }
  let reduced = Direction::new(0.0, angle).unwrap().phi
  let values = self.values()
  for m = 1; m <= self.degree; m = m + 1 {
    let a = m.to_double() * reduced
    let c = @math.cos(a)
    let s = @math.sin(a)
    for l = m; l <= self.degree; l = l + 1 {
      let positive = self.coefficients[index_of(l, m)]
      let negative = self.coefficients[index_of(l, -m)]
      values[index_of(l, m)] = c * positive - s * negative
      values[index_of(l, -m)] = s * positive + c * negative
    }
  }
  Spectrum::new(self.degree, values)
}

///|
/// Encode a weighted point direction: coefficients are amplitude*R_lm(direction).
/// This is a truncated delta distribution, not a continuous unit-energy point source.
pub fn encode_direction(
  degree : Int32,
  direction : Direction,
  amplitude : Double,
) -> Result[Spectrum, HarmonicError] {
  if !finite(amplitude) {
    return Err(InvalidInput("encoding amplitude must be finite"))
  }
  let values = match real_basis(degree, direction) {
    Ok(v) => v
    Err(e) => return Err(e)
  }
  for i = 0; i < values.length(); i = i + 1 {
    values[i] = values[i] * amplitude
  }
  Spectrum::new(degree, values)
}
