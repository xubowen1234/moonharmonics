///|
/// Immutable complete real spectrum. Constructors copy arrays; accessors never expose storage.
pub struct Spectrum {
  degree : Int32
  coefficients : Array[Double]
}

///|
/// Construct from exactly (degree+1)^2 finite coefficients, copying the input.
pub fn Spectrum::new(
  degree : Int32,
  coefficients : Array[Double],
) -> Result[Spectrum, HarmonicError] {
  let n = match coefficient_count(degree) {
    Ok(n) => n
    Err(e) => return Err(e)
  }
  if coefficients.length() != n {
    return Err(ShapeMismatch(n, coefficients.length()))
  }
  for v in coefficients {
    if !finite(v) {
      return Err(InvalidInput("non-finite coefficient"))
    }
  }
  Ok({ degree, coefficients: coefficients.copy(), })
}

///|
pub fn Spectrum::zeros(degree : Int32) -> Result[Spectrum, HarmonicError] {
  let n = match coefficient_count(degree) {
    Ok(n) => n
    Err(e) => return Err(e)
  }
  Ok({ degree, coefficients: Array::make(n, 0.0), })
}

///|
pub fn Spectrum::degree(self : Spectrum) -> Int32 {
  self.degree
}

///|
pub fn Spectrum::size(self : Spectrum) -> Int32 {
  self.coefficients.length()
}

///|
pub fn Spectrum::values(self : Spectrum) -> Array[Double] {
  self.coefficients.copy()
}

///|
pub fn Spectrum::coefficient(
  self : Spectrum,
  degree : Int32,
  order : Int32,
) -> Result[Double, HarmonicError] {
  let i = match coefficient_index(degree, order) {
    Ok(i) => i
    Err(e) => return Err(e)
  }
  if degree > self.degree {
    return Err(InvalidDegree(degree))
  }
  Ok(self.coefficients[i])
}

///|
/// Return a new spectrum with a single coefficient changed.
pub fn Spectrum::with_coefficient(
  self : Spectrum,
  degree : Int32,
  order : Int32,
  value : Double,
) -> Result[Spectrum, HarmonicError] {
  let i = match coefficient_index(degree, order) {
    Ok(i) => i
    Err(e) => return Err(e)
  }
  if degree > self.degree {
    return Err(InvalidDegree(degree))
  }
  if !finite(value) {
    return Err(InvalidInput("non-finite coefficient"))
  }
  let values = self.values()
  values[i] = value
  Ok({ degree: self.degree, coefficients: values, })
}

///|
/// Resize by truncation or exact zero-padding, never inferring a band limit from an array.
pub fn Spectrum::resize(
  self : Spectrum,
  degree : Int32,
) -> Result[Spectrum, HarmonicError] {
  let n = match coefficient_count(degree) {
    Ok(n) => n
    Err(e) => return Err(e)
  }
  let values = Array::make(n, 0.0)
  for i = 0; i < @cmp.minimum(n, self.size()); i = i + 1 {
    values[i] = self.coefficients[i]
  }
  Ok({ degree, coefficients: values, })
}

///|
/// Scale with overflow rejection; the input spectrum remains unchanged on failure.
pub fn Spectrum::scale(
  self : Spectrum,
  factor : Double,
) -> Result[Spectrum, HarmonicError] {
  if !finite(factor) {
    return Err(InvalidInput("non-finite scale"))
  }
  let values = self.values()
  for i = 0; i < values.length(); i = i + 1 {
    values[i] = values[i] * factor
  }
  Spectrum::new(self.degree, values)
}

///|
/// Add spectra at equal band limit, rejecting accidental truncation and overflow.
pub fn Spectrum::add(
  self : Spectrum,
  other : Spectrum,
) -> Result[Spectrum, HarmonicError] {
  if self.degree != other.degree {
    return Err(ShapeMismatch(self.size(), other.size()))
  }
  let values = self.values()
  for i = 0; i < values.length(); i = i + 1 {
    values[i] = values[i] + other.coefficients[i]
  }
  Spectrum::new(self.degree, values)
}
