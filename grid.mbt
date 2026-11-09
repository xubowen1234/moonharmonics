///|
/// Tensor quadrature: ascending cos(theta) rings and increasing uniform longitude.
/// Flat index is ring*nphi+longitude. All sample arrays use exactly this order.
pub struct SphereGrid {
  priv degree : Int32
  priv quadrature : Quadrature
  priv nphi : Int32
}

///|
/// The minimal product grid for bandlimited degree L uses L+1 rings and 2L+1 longitudes.
pub fn SphereGrid::new(degree : Int32) -> Result[SphereGrid, HarmonicError] {
  SphereGrid::oversampled(degree, degree + 1, 2 * degree + 1)
}

///|
/// Use additional nodes for non-bandlimited fields; exactness is not guaranteed for them.
pub fn SphereGrid::oversampled(
  degree : Int32,
  rings : Int32,
  longitudes : Int32,
) -> Result[SphereGrid, HarmonicError] {
  if !valid_degree(degree) {
    return Err(InvalidDegree(degree))
  }
  if rings < degree + 1 ||
    rings > 256 ||
    longitudes < 2 * degree + 1 ||
    longitudes > 1024 {
    return Err(
      InvalidInput("grid needs L+1..256 rings and 2L+1..1024 longitudes"),
    )
  }
  let quadrature = match gauss_legendre(rings) {
    Ok(q) => q
    Err(e) => return Err(e)
  }
  Ok({ degree, quadrature, nphi: longitudes, })
}

///|
pub fn SphereGrid::degree(self : SphereGrid) -> Int32 {
  self.degree
}

///|
pub fn SphereGrid::rings(self : SphereGrid) -> Int32 {
  self.quadrature.size()
}

///|
pub fn SphereGrid::longitudes(self : SphereGrid) -> Int32 {
  self.nphi
}

///|
pub fn SphereGrid::size(self : SphereGrid) -> Int32 {
  self.rings() * self.nphi
}

///|
pub fn SphereGrid::direction(
  self : SphereGrid,
  index : Int32,
) -> Result[Direction, HarmonicError] {
  if index < 0 || index >= self.size() {
    return Err(InvalidInput("grid index out of range"))
  }
  let x = self.quadrature.nodes[index / self.nphi]
  Direction::new(
    @math.acos(x),
    2.0 * @math.PI * (index % self.nphi).to_double() / self.nphi.to_double(),
  )
}

///|
pub fn SphereGrid::weight(
  self : SphereGrid,
  index : Int32,
) -> Result[Double, HarmonicError] {
  if index < 0 || index >= self.size() {
    return Err(InvalidInput("grid index out of range"))
  }
  Ok(
    self.quadrature.weights[index / self.nphi] *
    2.0 *
    @math.PI /
    self.nphi.to_double(),
  )
}

///|
/// Sample a callback in documented flat order. Non-finite callback output is rejected.
pub fn SphereGrid::sample(
  self : SphereGrid,
  field : (Direction) -> Double,
) -> Result[Array[Double], HarmonicError] {
  let values = Array::make(self.size(), 0.0)
  for i = 0; i < self.size(); i = i + 1 {
    let value = field(self.direction(i).unwrap())
    if !finite(value) {
      return Err(InvalidInput("sample callback returned non-finite value"))
    }
    values[i] = value
  }
  Ok(values)
}

///|
/// Integral over unit sphere using compensated weighted summation.
pub fn SphereGrid::integrate(
  self : SphereGrid,
  values : Array[Double],
) -> Result[Double, HarmonicError] {
  if values.length() != self.size() {
    return Err(ShapeMismatch(self.size(), values.length()))
  }
  let mut sum = 0.0
  let mut correction = 0.0
  for i = 0; i < values.length(); i = i + 1 {
    if !finite(values[i]) {
      return Err(InvalidInput("non-finite grid value"))
    }
    let term = values[i] * self.weight(i).unwrap() - correction
    let next = sum + term
    correction = next - sum - term
    sum = next
  }
  if !finite(sum) {
    Err(InvalidInput("sphere integral overflow"))
  } else {
    Ok(sum)
  }
}
