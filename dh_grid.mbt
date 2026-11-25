///|
/// Driscoll-Healy equiangular samples, with 2(L+1) rings and longitudes.
/// theta_j=pi*j/(2(L+1)); phi_k=pi*k/(L+1). The south pole is omitted.
pub struct DriscollHealyGrid {
  priv degree : Int32
  priv weights : Array[Double]
  priv n : Int32
}

///|
/// Create a DH sampling grid for bandlimited fields of degree at most L.
pub fn DriscollHealyGrid::new(
  degree : Int32,
) -> Result[DriscollHealyGrid, HarmonicError] {
  if !valid_degree(degree) {
    return Err(InvalidDegree(degree))
  }
  let n = degree + 1
  let weights = Array::make(2 * n, 0.0)
  for j = 0; j < 2 * n; j = j + 1 {
    let theta = @math.PI * j.to_double() / (2 * n).to_double()
    let mut sum = 0.0
    let mut correction = 0.0
    for k = 0; k < n; k = k + 1 {
      let odd = (2 * k + 1).to_double()
      let term = @math.sin(odd * theta) / odd - correction
      let next = sum + term
      correction = next - sum - term
      sum = next
    }
    let weight = 2.0 / n.to_double() * @math.sin(theta) * sum
    if !finite(weight) || weight < -1.0e-14 {
      return Err(InvalidInput("Driscoll-Healy quadrature construction failed"))
    }
    weights[j] = if weight < 0.0 { 0.0 } else { weight }
  }
  Ok({ degree, weights, n, })
}

///|
pub fn DriscollHealyGrid::degree(self : DriscollHealyGrid) -> Int32 {
  self.degree
}

///|
pub fn DriscollHealyGrid::rings(self : DriscollHealyGrid) -> Int32 {
  2 * self.n
}

///|
pub fn DriscollHealyGrid::longitudes(self : DriscollHealyGrid) -> Int32 {
  2 * self.n
}

///|
pub fn DriscollHealyGrid::size(self : DriscollHealyGrid) -> Int32 {
  4 * self.n * self.n
}

///|
pub fn DriscollHealyGrid::direction(
  self : DriscollHealyGrid,
  index : Int32,
) -> Result[Direction, HarmonicError] {
  if index < 0 || index >= self.size() {
    return Err(InvalidInput("grid index out of range"))
  }
  let ring = index / self.longitudes()
  let longitude = index % self.longitudes()
  Direction::new(
    @math.PI * ring.to_double() / (2 * self.n).to_double(),
    @math.PI * longitude.to_double() / self.n.to_double(),
  )
}

///|
pub fn DriscollHealyGrid::weight(
  self : DriscollHealyGrid,
  index : Int32,
) -> Result[Double, HarmonicError] {
  if index < 0 || index >= self.size() {
    return Err(InvalidInput("grid index out of range"))
  }
  let ring = index / self.longitudes()
  Ok(self.weights[ring] * 2.0 * @math.PI / self.longitudes().to_double())
}

///|
pub fn DriscollHealyGrid::sample(
  self : DriscollHealyGrid,
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
pub fn DriscollHealyGrid::integrate(
  self : DriscollHealyGrid,
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

///|
/// Analyze DH samples using the exact equiangular sampling quadrature.
pub fn DriscollHealyGrid::analyze(
  self : DriscollHealyGrid,
  values : Array[Double],
) -> Result[Spectrum, HarmonicError] {
  if values.length() != self.size() {
    return Err(ShapeMismatch(self.size(), values.length()))
  }
  for value in values {
    if !finite(value) {
      return Err(InvalidInput("non-finite analysis sample"))
    }
  }
  let result = Array::make((self.degree + 1) * (self.degree + 1), 0.0)
  let corrections = Array::make(result.length(), 0.0)
  let nphi = self.longitudes()
  let dphi = @math.PI / self.n.to_double()
  for ring = 0; ring < self.rings(); ring = ring + 1 {
    let cosine = Array::make(self.degree + 1, 0.0)
    let sine = Array::make(self.degree + 1, 0.0)
    let cc = Array::make(self.degree + 1, 0.0)
    let sc = Array::make(self.degree + 1, 0.0)
    for longitude = 0; longitude < nphi; longitude = longitude + 1 {
      let phi = dphi * longitude.to_double()
      let v = values[ring * nphi + longitude] * dphi
      for m = 0; m <= self.degree; m = m + 1 {
        let a = m.to_double() * phi
        let ct = v * @math.cos(a) - cc[m]
        let cn = cosine[m] + ct
        cc[m] = cn - cosine[m] - ct
        cosine[m] = cn
        let st = v * @math.sin(a) - sc[m]
        let sn = sine[m] + st
        sc[m] = sn - sine[m] - st
        sine[m] = sn
      }
    }
    let theta = @math.PI * ring.to_double() / (2 * self.n).to_double()
    let q = normalized_legendre(self.degree, @math.cos(theta)).unwrap()
    for l = 0; l <= self.degree; l = l + 1 {
      for m = -l; m <= l; m = m + 1 {
        let order = if m < 0 { -m } else { m }
        let fourier = if m < 0 { sine[order] } else { cosine[order] }
        let norm = if m == 0 { 1.0 } else { 2.0.sqrt() }
        let i = index_of(l, m)
        let term = self.weights[ring] * q[index_of(l, order)] * norm * fourier -
          corrections[i]
        let next = result[i] + term
        corrections[i] = next - result[i] - term
        result[i] = next
      }
    }
  }
  Spectrum::new(self.degree, result)
}

///|
pub fn DriscollHealyGrid::project(
  self : DriscollHealyGrid,
  field : (Direction) -> Double,
) -> Result[Spectrum, HarmonicError] {
  let values = match self.sample(field) {
    Ok(v) => v
    Err(e) => return Err(e)
  }
  self.analyze(values)
}

///|
/// Synthesize a bandlimited spectrum on DH equiangular samples.
pub fn Spectrum::synthesize_dh(
  self : Spectrum,
  grid : DriscollHealyGrid,
) -> Result[Array[Double], HarmonicError] {
  if self.degree > grid.degree {
    return Err(InvalidInput("grid band limit is below spectrum degree"))
  }
  let result = Array::make(grid.size(), 0.0)
  let nphi = grid.longitudes()
  let root_two = 2.0.sqrt()
  for ring = 0; ring < grid.rings(); ring = ring + 1 {
    let theta = @math.PI * ring.to_double() / (2 * grid.n).to_double()
    let q = normalized_legendre(self.degree, @math.cos(theta)).unwrap()
    let cosine = Array::make(self.degree + 1, 0.0)
    let sine = Array::make(self.degree + 1, 0.0)
    for l = 0; l <= self.degree; l = l + 1 {
      cosine[0] = cosine[0] +
        self.coefficients[index_of(l, 0)] * q[index_of(l, 0)]
      for m = 1; m <= l; m = m + 1 {
        let factor = root_two * q[index_of(l, m)]
        cosine[m] = cosine[m] + factor * self.coefficients[index_of(l, m)]
        sine[m] = sine[m] + factor * self.coefficients[index_of(l, -m)]
      }
    }
    for longitude = 0; longitude < nphi; longitude = longitude + 1 {
      let phi = @math.PI * longitude.to_double() / grid.n.to_double()
      let mut value = cosine[0]
      for m = 1; m <= self.degree; m = m + 1 {
        let angle = m.to_double() * phi
        value = value +
          cosine[m] * @math.cos(angle) +
          sine[m] * @math.sin(angle)
      }
      if !finite(value) {
        return Err(InvalidInput("grid synthesis overflow"))
      }
      result[ring * nphi + longitude] = value
    }
  }
  Ok(result)
}
