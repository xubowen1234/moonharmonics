///|
/// Versioned JSON with explicit normalization, phase and signed-order layout.
pub fn Spectrum::to_json(self : Spectrum) -> Json {
  Json::object({
    "format": Json::string("moonharmonics-real-v1"),
    "convention": Json::string("orthonormal-cs-l*l+l+m"),
    "degree": Json::number(self.degree.to_double()),
    "coefficients": Json::array(
      self.coefficients.map(fn(v) { Json::number(v) }),
    ),
  })
}

///|
pub fn Spectrum::to_json_string(self : Spectrum) -> String {
  self.to_json().stringify()
}

///|
/// Strict decoding: reject unknown/missing fields and convention mismatches.
/// JSON callers cannot omit the phase tag or rely on an inferred coefficient shape.
pub fn Spectrum::from_json(value : Json) -> Result[Spectrum, HarmonicError] {
  let object = match value {
    Object(o) => o
    _ => return Err(InvalidInput("spectrum must be a JSON object"))
  }
  if object.length() != 4 {
    return Err(InvalidInput("spectrum JSON needs exactly four fields"))
  }
  match object.get("format") {
    Some(String("moonharmonics-real-v1")) => ()
    _ => return Err(InvalidInput("unsupported spectrum format"))
  }
  match object.get("convention") {
    Some(String("orthonormal-cs-l*l+l+m")) => ()
    _ => return Err(InvalidInput("unsupported coefficient convention"))
  }
  let degree = match object.get("degree") {
    Some(Number(d, ..)) => {
      if !finite(d) || d < 0.0 || d > MAX_DEGREE.to_double() || d.floor() != d {
        return Err(InvalidInput("invalid JSON degree"))
      }
      d.to_int()
    }
    _ => return Err(InvalidInput("JSON degree must be a number"))
  }
  let raw = match object.get("coefficients") {
    Some(Array(a)) => a
    _ => return Err(InvalidInput("JSON coefficients must be an array"))
  }
  let n = (degree + 1) * (degree + 1)
  if raw.length() != n {
    return Err(ShapeMismatch(n, raw.length()))
  }
  let coefficients = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    match raw[i] {
      Number(d, ..) => coefficients[i] = d
      _ => return Err(InvalidInput("coefficient must be numeric"))
    }
  }
  Spectrum::new(degree, coefficients)
}

///|
/// Parse bounded input and reject duplicate keys rather than accepting ambiguous metadata.
pub fn Spectrum::from_json_string(
  text : String,
) -> Result[Spectrum, HarmonicError] {
  if text.length() > 2_000_000 {
    return Err(InvalidInput("spectrum JSON exceeds two million UTF-16 units"))
  }
  let json = @json.parse(text, max_nesting_depth=8, reject_duplicate_keys=true) catch {
    _ => return Err(InvalidInput("invalid spectrum JSON"))
  }
  Spectrum::from_json(json)
}

///|
/// Convert from phase-free SN3D real harmonics in the same packed signed-order layout.
/// Input coefficients expand S_lm = (-1)^|m| sqrt(4pi/(2l+1)) R_lm.
/// This converts expansion coefficients, not sampled basis values.
pub fn Spectrum::from_sn3d(
  degree : Int32,
  coefficients : Array[Double],
) -> Result[Spectrum, HarmonicError] {
  let base = match Spectrum::new(degree, coefficients) {
    Ok(s) => s
    Err(e) => return Err(e)
  }
  let values = base.values()
  for l = 0; l <= degree; l = l + 1 {
    let factor = (4.0 * @math.PI / (2 * l + 1).to_double()).sqrt()
    for m = -l; m <= l; m = m + 1 {
      let sign = if m % 2 == 0 { 1.0 } else { -1.0 }
      let i = index_of(l, m)
      values[i] = values[i] * factor * sign
    }
  }
  Spectrum::new(degree, values)
}

///|
/// Expansion coefficients for phase-free SN3D, with identical signed-order indexing.
pub fn Spectrum::to_sn3d(self : Spectrum) -> Array[Double] {
  let values = self.values()
  for l = 0; l <= self.degree; l = l + 1 {
    let factor = ((2 * l + 1).to_double() / (4.0 * @math.PI)).sqrt()
    for m = -l; m <= l; m = m + 1 {
      let sign = if m % 2 == 0 { 1.0 } else { -1.0 }
      let i = index_of(l, m)
      values[i] = values[i] * factor * sign
    }
  }
  values
}
