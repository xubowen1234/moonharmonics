///|
/// A rejected argument or an operation that cannot satisfy its mathematical contract.
pub(all) enum HarmonicError {
  InvalidDegree(Int32)
  InvalidOrder(Int32, Int32)
  InvalidInput(String)
  ShapeMismatch(Int32, Int32)
  SingularSystem(Int32)
} derive(Eq, @debug.Debug)

///|
pub extend HarmonicError with Eq::{equal, not_equal}

///|
pub extend HarmonicError with @debug.Debug::{to_repr}

///|
/// Largest supported band limit. Storage contains (degree + 1)^2 real coefficients.
pub const MAX_DEGREE : Int32 = 64

///|
fn valid_degree(degree : Int32) -> Bool {
  degree >= 0 && degree <= MAX_DEGREE
}

///|
/// ACN-like contiguous layout, with signed order: index = l*l + l + m.
pub fn coefficient_index(
  degree : Int32,
  order : Int32,
) -> Result[Int32, HarmonicError] {
  if !valid_degree(degree) {
    Err(InvalidDegree(degree))
  } else if order < -degree || order > degree {
    Err(InvalidOrder(degree, order))
  } else {
    Ok(index_of(degree, order))
  }
}

///|
fn index_of(degree : Int32, order : Int32) -> Int32 {
  degree * degree + degree + order
}

///|
/// Convert a packed index back to (l,m), rejecting indices outside the supported model.
pub fn degree_order(index : Int32) -> Result[(Int32, Int32), HarmonicError] {
  if index < 0 || index >= (MAX_DEGREE + 1) * (MAX_DEGREE + 1) {
    return Err(InvalidInput("coefficient index out of range"))
  }
  let mut degree = 0
  while (degree + 1) * (degree + 1) <= index {
    degree = degree + 1
  }
  Ok((degree, index - degree * degree - degree))
}

///|
/// Exact number of coefficients including all bands from zero to degree.
pub fn coefficient_count(degree : Int32) -> Result[Int32, HarmonicError] {
  if valid_degree(degree) {
    Ok((degree + 1) * (degree + 1))
  } else {
    Err(InvalidDegree(degree))
  }
}
