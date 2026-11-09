///|
/// Reconstruct one scalar using compensated summation. Overflow is an error.
pub fn Spectrum::evaluate(
  self : Spectrum,
  direction : Direction,
) -> Result[Double, HarmonicError] {
  let basis = real_basis(self.degree, direction).unwrap()
  let mut sum = 0.0
  let mut correction = 0.0
  for i = 0; i < self.size(); i = i + 1 {
    let term = self.coefficients[i] * basis[i] - correction
    let next = sum + term
    correction = next - sum - term
    sum = next
  }
  if !finite(sum) {
    Err(InvalidInput("synthesis overflow"))
  } else {
    Ok(sum)
  }
}

///|
/// Reconstruct independent arbitrary directions, preserving input order.
pub fn Spectrum::evaluate_many(
  self : Spectrum,
  directions : Array[Direction],
) -> Result[Array[Double], HarmonicError] {
  let output = Array::make(directions.length(), 0.0)
  for i = 0; i < directions.length(); i = i + 1 {
    output[i] = match self.evaluate(directions[i]) {
      Ok(v) => v
      Err(e) => return Err(e)
    }
  }
  Ok(output)
}

///|
/// Separable synthesis on a tensor grid. Legendre factors are reused per ring.
pub fn Spectrum::synthesize(
  self : Spectrum,
  grid : SphereGrid,
) -> Result[Array[Double], HarmonicError] {
  if self.degree > grid.degree {
    return Err(InvalidInput("grid band limit is below spectrum degree"))
  }
  let result = Array::make(grid.size(), 0.0)
  let root_two = 2.0.sqrt()
  for ring = 0; ring < grid.rings(); ring = ring + 1 {
    let q = normalized_legendre(self.degree, grid.quadrature.nodes[ring]).unwrap()
    let cos_amplitudes = Array::make(self.degree + 1, 0.0)
    let sin_amplitudes = Array::make(self.degree + 1, 0.0)
    for l = 0; l <= self.degree; l = l + 1 {
      cos_amplitudes[0] = cos_amplitudes[0] +
        self.coefficients[index_of(l, 0)] * q[index_of(l, 0)]
      for m = 1; m <= l; m = m + 1 {
        let factor = root_two * q[index_of(l, m)]
        cos_amplitudes[m] = cos_amplitudes[m] +
          factor * self.coefficients[index_of(l, m)]
        sin_amplitudes[m] = sin_amplitudes[m] +
          factor * self.coefficients[index_of(l, -m)]
      }
    }
    for longitude = 0; longitude < grid.nphi; longitude = longitude + 1 {
      let phi = 2.0 * @math.PI * longitude.to_double() / grid.nphi.to_double()
      let mut value = cos_amplitudes[0]
      for m = 1; m <= self.degree; m = m + 1 {
        let angle = m.to_double() * phi
        value = value +
          cos_amplitudes[m] * @math.cos(angle) +
          sin_amplitudes[m] * @math.sin(angle)
      }
      if !finite(value) {
        return Err(InvalidInput("grid synthesis overflow"))
      }
      result[ring * grid.nphi + longitude] = value
    }
  }
  Ok(result)
}
