///|
/// Sum of squared coefficients per degree; orthonormal normalization makes this energy.
pub fn Spectrum::band_power(
  self : Spectrum,
) -> Result[Array[Double], HarmonicError] {
  let power = Array::make(self.degree + 1, 0.0)
  for l = 0; l <= self.degree; l = l + 1 {
    let mut sum = 0.0
    let mut correction = 0.0
    for m = -l; m <= l; m = m + 1 {
      let v = self.coefficients[index_of(l, m)]
      let term = v * v - correction
      let next = sum + term
      correction = next - sum - term
      sum = next
    }
    if !finite(sum) {
      return Err(InvalidInput("band energy overflow"))
    }
    power[l] = sum
  }
  Ok(power)
}

///|
/// Integral of f^2 over the unit sphere (not a sample average).
pub fn Spectrum::energy(self : Spectrum) -> Result[Double, HarmonicError] {
  let power = match self.band_power() {
    Ok(p) => p
    Err(e) => return Err(e)
  }
  let mut total = 0.0
  for p in power {
    total = total + p
  }
  if !finite(total) {
    Err(InvalidInput("total energy overflow"))
  } else {
    Ok(total)
  }
}

///|
pub fn Spectrum::mean(self : Spectrum) -> Double {
  self.coefficients[0] / (4.0 * @math.PI).sqrt()
}

///|
/// RMS with respect to uniform solid angle, not unweighted latitude samples.
pub fn Spectrum::rms(self : Spectrum) -> Result[Double, HarmonicError] {
  match self.energy() {
    Ok(e) => Ok((e / (4.0 * @math.PI)).sqrt())
    Err(e) => Err(e)
  }
}

///|
/// H1 seminorm squared, integral |surface-gradient f|^2.
pub fn Spectrum::gradient_energy(
  self : Spectrum,
) -> Result[Double, HarmonicError] {
  let power = match self.band_power() {
    Ok(p) => p
    Err(e) => return Err(e)
  }
  let mut total = 0.0
  for l = 1; l <= self.degree; l = l + 1 {
    total = total + (l * (l + 1)).to_double() * power[l]
  }
  if !finite(total) {
    Err(InvalidInput("gradient energy overflow"))
  } else {
    Ok(total)
  }
}

///|
/// L2 distance between equal-bandlimit spectra.
pub fn Spectrum::distance(
  self : Spectrum,
  other : Spectrum,
) -> Result[Double, HarmonicError] {
  if self.degree != other.degree {
    return Err(ShapeMismatch(self.size(), other.size()))
  }
  let mut scale = 0.0
  let mut sum = 1.0
  // Scaled sum of squares avoids intermediate squaring overflow.
  for i = 0; i < self.size(); i = i + 1 {
    let value = (self.coefficients[i] - other.coefficients[i]).abs()
    if !finite(value) {
      return Err(InvalidInput("distance subtraction overflow"))
    }
    if value != 0.0 {
      if scale < value {
        let ratio = scale / value
        sum = 1.0 + sum * ratio * ratio
        scale = value
      } else {
        let ratio = value / scale
        sum = sum + ratio * ratio
      }
    }
  }
  let answer = scale * sum.sqrt()
  if !finite(answer) {
    Err(InvalidInput("distance overflow"))
  } else {
    Ok(answer)
  }
}

///|
/// Relative Parseval defect for samples believed to represent this bandlimited spectrum.
pub fn Spectrum::parseval_defect(
  self : Spectrum,
  grid : SphereGrid,
  values : Array[Double],
) -> Result[Double, HarmonicError] {
  if values.length() != grid.size() {
    return Err(ShapeMismatch(grid.size(), values.length()))
  }
  let squared = values.map(fn(v) { v * v })
  let sampled = match grid.integrate(squared) {
    Ok(e) => e
    Err(e) => return Err(e)
  }
  let exact = match self.energy() {
    Ok(e) => e
    Err(e) => return Err(e)
  }
  Ok((sampled - exact).abs() / @cmp.maximum(1.0, exact))
}
