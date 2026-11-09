///|
/// A streaming weighted projection accumulator, with compensated sums per coefficient.
/// Weights represent solid-angle quadrature, not statistical confidence weights.
/// Arbitrary samples do not imply exact orthogonality; use fit_weighted for regression.
pub struct ProjectionAccumulator {
  degree : Int32
  sums : Array[Double]
  corrections : Array[Double]
  mut weight : Double
  mut weight_correction : Double
  mut samples : Int32
}

///|
pub fn ProjectionAccumulator::new(
  degree : Int32,
) -> Result[ProjectionAccumulator, HarmonicError] {
  let count = match coefficient_count(degree) {
    Ok(n) => n
    Err(e) => return Err(e)
  }
  Ok({
    degree,
    sums: Array::make(count, 0.0),
    corrections: Array::make(count, 0.0),
    weight: 0.0,
    weight_correction: 0.0,
    samples: 0,
  })
}

///|
/// Add one finite sample with nonnegative solid angle. On error, state remains unchanged.
pub fn ProjectionAccumulator::push(
  self : ProjectionAccumulator,
  direction : Direction,
  value : Double,
  weight : Double,
) -> Result[Unit, HarmonicError] {
  if !finite(value) || !finite(weight) || weight < 0.0 {
    return Err(InvalidInput("invalid weighted projection sample"))
  }
  if self.samples == 2_147_483_647 {
    return Err(InvalidInput("projection sample counter overflow"))
  }
  let factors = real_basis(self.degree, direction).unwrap()
  let next_sums = self.sums.copy()
  let next_corrections = self.corrections.copy()
  for i = 0; i < factors.length(); i = i + 1 {
    let term = value * weight * factors[i] - self.corrections[i]
    let next = self.sums[i] + term
    if !finite(next) {
      return Err(InvalidInput("projection coefficient overflow"))
    }
    next_corrections[i] = next - self.sums[i] - term
    next_sums[i] = next
  }
  let term = weight - self.weight_correction
  let next_weight = self.weight + term
  if !finite(next_weight) {
    return Err(InvalidInput("projection weight overflow"))
  }
  for i = 0; i < factors.length(); i = i + 1 {
    self.sums[i] = next_sums[i]
    self.corrections[i] = next_corrections[i]
  }
  self.weight_correction = next_weight - self.weight - term
  self.weight = next_weight
  self.samples = self.samples + 1
  Ok(())
}

///|
pub fn ProjectionAccumulator::samples(self : ProjectionAccumulator) -> Int32 {
  self.samples
}

///|
pub fn ProjectionAccumulator::weight(self : ProjectionAccumulator) -> Double {
  self.weight
}

///|
/// Snapshot coefficients. No division by total weight: weights already encode solid angle.
pub fn ProjectionAccumulator::spectrum(
  self : ProjectionAccumulator,
) -> Result[Spectrum, HarmonicError] {
  if self.samples == 0 || self.weight == 0.0 {
    return Err(InvalidInput("projection has no positive measure"))
  }
  Spectrum::new(self.degree, self.sums)
}

///|
/// Combine an independent partial integral. No shared mutable arrays are introduced.
/// Both accumulator states are unchanged if the combined values overflow.
pub fn ProjectionAccumulator::merge(
  self : ProjectionAccumulator,
  other : ProjectionAccumulator,
) -> Result[Unit, HarmonicError] {
  if self.degree != other.degree {
    return Err(InvalidInput("projection band limits differ"))
  }
  if other.samples > 2_147_483_647 - self.samples {
    return Err(InvalidInput("projection sample counter overflow"))
  }
  let sums = self.sums.copy()
  let corrections = self.corrections.copy()
  for i = 0; i < sums.length(); i = i + 1 {
    let term = other.sums[i] - self.corrections[i]
    let next = self.sums[i] + term
    if !finite(next) {
      return Err(InvalidInput("merged projection overflow"))
    }
    corrections[i] = next - self.sums[i] - term
    sums[i] = next
  }
  let term = other.weight - self.weight_correction
  let next_weight = self.weight + term
  if !finite(next_weight) {
    return Err(InvalidInput("merged weight overflow"))
  }
  for i = 0; i < sums.length(); i = i + 1 {
    self.sums[i] = sums[i]
    self.corrections[i] = corrections[i]
  }
  self.weight_correction = next_weight - self.weight - term
  self.weight = next_weight
  self.samples = self.samples + other.samples
  Ok(())
}

///|
/// Integrate a sampled field against another sampled field on a sphere grid.
pub fn SphereGrid::inner_product(
  self : SphereGrid,
  a : Array[Double],
  b : Array[Double],
) -> Result[Double, HarmonicError] {
  if a.length() != self.size() {
    return Err(ShapeMismatch(self.size(), a.length()))
  }
  if b.length() != self.size() {
    return Err(ShapeMismatch(self.size(), b.length()))
  }
  let products = Array::make(self.size(), 0.0)
  for i = 0; i < products.length(); i = i + 1 {
    products[i] = a[i] * b[i]
  }
  self.integrate(products)
}
