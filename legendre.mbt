///|
/// Fully normalized associated Legendre functions q_lm for nonnegative order.
/// Entries with negative order are zero; the layout is shared with real coefficients.
/// The Condon-Shortley phase is included. No factorials are formed.
pub fn normalized_legendre(
  degree : Int32,
  x : Double,
) -> Result[Array[Double], HarmonicError] {
  if !valid_degree(degree) {
    return Err(InvalidDegree(degree))
  }
  if !finite(x) || x < -1.0 || x > 1.0 {
    return Err(InvalidInput("Legendre argument outside [-1,1]"))
  }
  let values = Array::make((degree + 1) * (degree + 1), 0.0)
  values[0] = 1.0 / (4.0 * @math.PI).sqrt()
  let s = @cmp.maximum(0.0, (1.0 - x) * (1.0 + x)).sqrt()
  for m = 0; m <= degree; m = m + 1 {
    if m > 0 {
      let md = m.to_double()
      values[index_of(m, m)] = -((2.0 * md + 1.0) / (2.0 * md)).sqrt() *
        s *
        values[index_of(m - 1, m - 1)]
    }
    if m < degree {
      values[index_of(m + 1, m)] = (2.0 * m.to_double() + 3.0).sqrt() *
        x *
        values[index_of(m, m)]
    }
    for l = m + 2; l <= degree; l = l + 1 {
      let ld = l.to_double()
      let md = m.to_double()
      let denom = ld * ld - md * md
      let a = ((4.0 * ld * ld - 1.0) / denom).sqrt()
      let b = ((2.0 * ld + 1.0) *
      ((ld - 1.0) * (ld - 1.0) - md * md) /
      ((2.0 * ld - 3.0) * denom)).sqrt()
      values[index_of(l, m)] = a * x * values[index_of(l - 1, m)] -
        b * values[index_of(l - 2, m)]
    }
  }
  Ok(values)
}
