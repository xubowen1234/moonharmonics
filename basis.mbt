///|
/// Evaluate the real orthonormal basis at one direction.
/// Positive orders are cosine terms, negative orders sine terms, both with CS phase.
pub fn real_basis(
  degree : Int32,
  direction : Direction,
) -> Result[Array[Double], HarmonicError] {
  if !valid_degree(degree) {
    return Err(InvalidDegree(degree))
  }
  let q = legendre_factors(
    degree,
    @math.cos(direction.theta),
    @math.sin(direction.theta),
  )
  let root_two = 2.0.sqrt()
  for m = 1; m <= degree; m = m + 1 {
    let angle = m.to_double() * direction.phi
    let c = root_two * @math.cos(angle)
    let s = root_two * @math.sin(angle)
    for l = m; l <= degree; l = l + 1 {
      let value = q[index_of(l, m)]
      q[index_of(l, m)] = value * c
      q[index_of(l, -m)] = value * s
    }
  }
  Ok(q)
}

///|
/// Evaluate one real harmonic. Uses the same signed-order normalization as real_basis.
pub fn real_harmonic(
  degree : Int32,
  order : Int32,
  direction : Direction,
) -> Result[Double, HarmonicError] {
  let i = match coefficient_index(degree, order) {
    Ok(i) => i
    Err(e) => return Err(e)
  }
  match real_basis(degree, direction) {
    Ok(v) => Ok(v[i])
    Err(e) => Err(e)
  }
}

///|
/// Addition-theorem check value: sum over one band's squared basis values.
pub fn basis_band_norm(
  degree : Int32,
  direction : Direction,
) -> Result[Double, HarmonicError] {
  let b = match real_basis(degree, direction) {
    Ok(v) => v
    Err(e) => return Err(e)
  }
  let mut sum = 0.0
  for m = -degree; m <= degree; m = m + 1 {
    let v = b[index_of(degree, m)]
    sum = sum + v * v
  }
  Ok(sum)
}
