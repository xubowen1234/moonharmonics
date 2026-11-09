///|
/// Zonal kernel represented by degree multipliers k_l=2pi*integral K(t)P_l(t)dt.
pub struct ZonalKernel {
  priv gains : Array[Double]
}

///|
pub fn ZonalKernel::new(
  gains : Array[Double],
) -> Result[ZonalKernel, HarmonicError] {
  if gains.length() < 1 || gains.length() > MAX_DEGREE + 1 {
    return Err(InvalidInput("kernel needs 1..65 gains"))
  }
  for v in gains {
    if !finite(v) {
      return Err(InvalidInput("non-finite kernel gain"))
    }
  }
  Ok({ gains: gains.copy(), })
}

///|
pub fn ZonalKernel::degree(self : ZonalKernel) -> Int32 {
  self.gains.length() - 1
}

///|
pub fn ZonalKernel::gains(self : ZonalKernel) -> Array[Double] {
  self.gains.copy()
}

///|
/// Construct from orthonormal m=0 coefficients, which are not convolution multipliers.
pub fn ZonalKernel::from_zonal_coefficients(
  coefficients : Array[Double],
) -> Result[ZonalKernel, HarmonicError] {
  let gains = coefficients.copy()
  for l = 0; l < gains.length(); l = l + 1 {
    gains[l] = gains[l] * (4.0 * @math.PI / (2 * l + 1).to_double()).sqrt()
  }
  ZonalKernel::new(gains)
}

///|
/// Project a smooth kernel K(dot) by quadrature. Increase nodes and check convergence for nonsmooth K.
pub fn ZonalKernel::project(
  degree : Int32,
  nodes : Int32,
  kernel : (Double) -> Double,
) -> Result[ZonalKernel, HarmonicError] {
  if !valid_degree(degree) {
    return Err(InvalidDegree(degree))
  }
  if nodes < degree + 1 {
    return Err(InvalidInput("kernel quadrature needs at least degree+1 nodes"))
  }
  let q = match gauss_legendre(nodes) {
    Ok(q) => q
    Err(e) => return Err(e)
  }
  let gains = Array::make(degree + 1, 0.0)
  for i = 0; i < q.size(); i = i + 1 {
    let x = q.nodes[i]
    let v = kernel(x)
    if !finite(v) {
      return Err(InvalidInput("non-finite kernel callback"))
    }
    for l = 0; l <= degree; l = l + 1 {
      let (p, _) = legendre_pair(l, x)
      gains[l] = gains[l] + 2.0 * @math.PI * q.weights[i] * v * p
    }
  }
  ZonalKernel::new(gains)
}

///|
/// Exact degree multipliers for K(t)=max(0,t), including all supported even degrees.
pub fn ZonalKernel::lambertian(
  degree : Int32,
) -> Result[ZonalKernel, HarmonicError] {
  if !valid_degree(degree) {
    return Err(InvalidDegree(degree))
  }
  let gains = Array::make(degree + 1, 0.0)
  gains[0] = @math.PI
  if degree >= 1 {
    gains[1] = 2.0 * @math.PI / 3.0
  }
  for l = 2; l <= degree; l = l + 2 {
    let (p, _) = legendre_pair(l, 0.0)
    gains[l] = -2.0 * @math.PI * p / ((l - 1) * (l + 2)).to_double()
  }
  ZonalKernel::new(gains)
}

///|
/// Convolution integral over solid angle; no extra pi or albedo factor is applied.
pub fn Spectrum::convolve(
  self : Spectrum,
  kernel : ZonalKernel,
) -> Result[Spectrum, HarmonicError] {
  if self.degree != kernel.degree() {
    return Err(ShapeMismatch(self.degree + 1, kernel.gains.length()))
  }
  self.filter_bands(kernel.gains)
}
