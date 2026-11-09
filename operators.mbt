///|
/// Laplace-Beltrami operator on the unit sphere: eigenvalue -l(l+1).
pub fn Spectrum::laplacian(self : Spectrum) -> Result[Spectrum, HarmonicError] {
  let gains = Array::make(self.degree + 1, 0.0)
  for l = 1; l <= self.degree; l = l + 1 {
    gains[l] = -(l * (l + 1)).to_double()
  }
  self.filter_bands(gains)
}

///|
/// Solve Delta(u)=self with zero-mean gauge. Reject a nonzero source mean.
/// mean_tolerance is an absolute tolerance on a00; accepted residual is removed explicitly.
pub fn Spectrum::poisson(
  self : Spectrum,
  mean_tolerance : Double,
) -> Result[Spectrum, HarmonicError] {
  if !finite(mean_tolerance) || mean_tolerance < 0.0 {
    return Err(InvalidInput("invalid mean tolerance"))
  }
  if self.coefficients[0].abs() > mean_tolerance {
    return Err(InvalidInput("Poisson source must have zero mean"))
  }
  let gains = Array::make(self.degree + 1, 0.0)
  for l = 1; l <= self.degree; l = l + 1 {
    gains[l] = -1.0 / (l * (l + 1)).to_double()
  }
  self.filter_bands(gains)
}

///|
/// Solve (identity - time*Delta)u=f. Finite nonnegative time keeps the solve coercive.
pub fn Spectrum::screened_poisson(
  self : Spectrum,
  time : Double,
) -> Result[Spectrum, HarmonicError] {
  if !finite(time) || time < 0.0 {
    return Err(InvalidInput("screening time must be finite and nonnegative"))
  }
  let gains = Array::make(self.degree + 1, 1.0)
  for l = 1; l <= self.degree; l = l + 1 {
    gains[l] = 1.0 / (1.0 + time * (l * (l + 1)).to_double())
  }
  self.filter_bands(gains)
}

///|
/// Azimuthal coordinate derivative, represented exactly within the same band limit.
pub fn Spectrum::differentiate_phi(
  self : Spectrum,
) -> Result[Spectrum, HarmonicError] {
  let values = Array::make(self.size(), 0.0)
  for l = 1; l <= self.degree; l = l + 1 {
    for m = 1; m <= l; m = m + 1 {
      values[index_of(l, m)] = m.to_double() *
        self.coefficients[index_of(l, -m)]
      values[index_of(l, -m)] = -m.to_double() *
        self.coefficients[index_of(l, m)]
    }
  }
  Spectrum::new(self.degree, values)
}
