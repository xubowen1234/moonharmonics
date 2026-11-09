///|
/// Orthogonal projection of grid samples onto the grid's declared band limit.
/// Longitude Fourier sums are reused for each ring, then integrated in cos(theta).
/// Non-bandlimited input may alias: use an oversampled grid and convergence checks.
pub fn SphereGrid::analyze(
  self : SphereGrid,
  values : Array[Double],
) -> Result[Spectrum, HarmonicError] {
  if values.length() != self.size() {
    return Err(ShapeMismatch(self.size(), values.length()))
  }
  for value in values {
    if !finite(value) {
      return Err(InvalidInput("non-finite analysis sample"))
    }
  }
  let result = Array::make((self.degree + 1) * (self.degree + 1), 0.0)
  let corrections = Array::make(result.length(), 0.0)
  let dphi = 2.0 * @math.PI / self.nphi.to_double()
  for ring = 0; ring < self.rings(); ring = ring + 1 {
    let cosine = Array::make(self.degree + 1, 0.0)
    let sine = Array::make(self.degree + 1, 0.0)
    let cc = Array::make(self.degree + 1, 0.0)
    let sc = Array::make(self.degree + 1, 0.0)
    for longitude = 0; longitude < self.nphi; longitude = longitude + 1 {
      let phi = dphi * longitude.to_double()
      let v = values[ring * self.nphi + longitude] * dphi
      for m = 0; m <= self.degree; m = m + 1 {
        let a = m.to_double() * phi
        let ct = v * @math.cos(a) - cc[m]
        let cn = cosine[m] + ct
        cc[m] = cn - cosine[m] - ct
        cosine[m] = cn
        let st = v * @math.sin(a) - sc[m]
        let sn = sine[m] + st
        sc[m] = sn - sine[m] - st
        sine[m] = sn
      }
    }
    let q = normalized_legendre(self.degree, self.quadrature.nodes[ring]).unwrap()
    let weight = self.quadrature.weights[ring]
    for l = 0; l <= self.degree; l = l + 1 {
      for m = -l; m <= l; m = m + 1 {
        let order = if m < 0 { -m } else { m }
        let fourier = if m < 0 { sine[order] } else { cosine[order] }
        let normalization = if m == 0 { 1.0 } else { 2.0.sqrt() }
        let i = index_of(l, m)
        let term = weight * q[index_of(l, order)] * normalization * fourier -
          corrections[i]
        let next = result[i] + term
        corrections[i] = next - result[i] - term
        result[i] = next
      }
    }
  }
  Spectrum::new(self.degree, result)
}

///|
/// Sample then project a user field. The callback is evaluated once per grid node.
pub fn SphereGrid::project(
  self : SphereGrid,
  field : (Direction) -> Double,
) -> Result[Spectrum, HarmonicError] {
  let values = match self.sample(field) {
    Ok(v) => v
    Err(e) => return Err(e)
  }
  self.analyze(values)
}
