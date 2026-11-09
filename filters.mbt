///|
/// Apply one scalar gain to each complete degree band.
pub fn Spectrum::filter_bands(
  self : Spectrum,
  gains : Array[Double],
) -> Result[Spectrum, HarmonicError] {
  if gains.length() != self.degree + 1 {
    return Err(ShapeMismatch(self.degree + 1, gains.length()))
  }
  for gain in gains {
    if !finite(gain) {
      return Err(InvalidInput("non-finite band gain"))
    }
  }
  let values = self.values()
  for l = 0; l <= self.degree; l = l + 1 {
    for m = -l; m <= l; m = m + 1 {
      let i = index_of(l, m)
      values[i] = values[i] * gains[l]
    }
  }
  Spectrum::new(self.degree, values)
}

///|
/// Solve the heat equation for finite nonnegative time: exp(-t*l*(l+1)).
pub fn Spectrum::heat(
  self : Spectrum,
  time : Double,
) -> Result[Spectrum, HarmonicError] {
  if !finite(time) || time < 0.0 {
    return Err(InvalidInput("heat time must be finite and nonnegative"))
  }
  let gains = Array::make(self.degree + 1, 1.0)
  for l = 1; l <= self.degree; l = l + 1 {
    gains[l] = @math.exp(-time * (l * (l + 1)).to_double())
  }
  self.filter_bands(gains)
}

///|
/// Keep degrees in an inclusive interval. Preserve the declared coefficient shape.
pub fn Spectrum::bandpass(
  self : Spectrum,
  low : Int32,
  high : Int32,
) -> Result[Spectrum, HarmonicError] {
  if low < 0 || high < low || high > self.degree {
    return Err(InvalidInput("invalid inclusive pass band"))
  }
  let gains = Array::make(self.degree + 1, 0.0)
  for l = low; l <= high; l = l + 1 {
    gains[l] = 1.0
  }
  self.filter_bands(gains)
}

///|
/// Smooth low pass: unit gain through start, cosine roll-off to zero at stop.
pub fn Spectrum::cosine_taper(
  self : Spectrum,
  start : Int32,
  stop : Int32,
) -> Result[Spectrum, HarmonicError] {
  if start < 0 || stop <= start || stop > self.degree {
    return Err(InvalidInput("taper requires 0<=start<stop<=degree"))
  }
  let gains = Array::make(self.degree + 1, 0.0)
  for l = 0; l <= self.degree; l = l + 1 {
    if l <= start {
      gains[l] = 1.0
    } else if l < stop {
      gains[l] = 0.5 *
        (
          1.0 +
          @math.cos(
            @math.PI * (l - start).to_double() / (stop - start).to_double(),
          )
        )
    }
  }
  self.filter_bands(gains)
}

///|
/// Remove the mean explicitly; useful before a Poisson solve when desired by the caller.
pub fn Spectrum::without_mean(self : Spectrum) -> Spectrum {
  let values = self.values()
  values[0] = 0.0
  { degree: self.degree, coefficients: values, }
}
