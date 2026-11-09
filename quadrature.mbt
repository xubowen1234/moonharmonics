///|
/// Nodes and weights for integrating a scalar function on [-1,1].
pub struct Quadrature {
  nodes : Array[Double]
  weights : Array[Double]
}

///|
pub fn Quadrature::size(self : Quadrature) -> Int32 {
  self.nodes.length()
}

///|
pub fn Quadrature::nodes(self : Quadrature) -> Array[Double] {
  self.nodes.copy()
}

///|
pub fn Quadrature::weights(self : Quadrature) -> Array[Double] {
  self.weights.copy()
}

///|
/// Gauss-Legendre quadrature, 1..256 nodes. Newton iteration is bounded and checked.
/// Nodes are sorted from negative to positive. Polynomial exactness is degree 2n-1.
pub fn gauss_legendre(count : Int32) -> Result[Quadrature, HarmonicError] {
  if count < 1 || count > 256 {
    return Err(InvalidInput("quadrature count must be in 1..256"))
  }
  let nodes = Array::make(count, 0.0)
  let weights = Array::make(count, 0.0)
  let nd = count.to_double()
  for i = 0; i < (count + 1) / 2; i = i + 1 {
    let mut z = @math.cos(@math.PI * (i.to_double() + 0.75) / (nd + 0.5))
    let mut converged = false
    for iteration = 0; iteration < 64; iteration = iteration + 1 {
      let (pn, previous) = legendre_pair(count, z)
      let derivative = nd * (z * pn - previous) / (z * z - 1.0)
      let next = z - pn / derivative
      if (next - z).abs() <= 4.0e-16 {
        z = next
        converged = true
        break
      }
      z = next
    }
    if !converged {
      return Err(InvalidInput("Gauss-Legendre root did not converge"))
    }
    let (pn, previous) = legendre_pair(count, z)
    let derivative = nd * (z * pn - previous) / (z * z - 1.0)
    let weight = 2.0 / ((1.0 - z * z) * derivative * derivative)
    nodes[i] = -z
    nodes[count - 1 - i] = z
    weights[i] = weight
    weights[count - 1 - i] = weight
  }
  Ok({ nodes, weights, })
}

///|
fn legendre_pair(degree : Int32, x : Double) -> (Double, Double) {
  if degree == 0 {
    return (1.0, 0.0)
  }
  let mut previous = 1.0
  let mut current = x
  for l = 2; l <= degree; l = l + 1 {
    let ld = l.to_double()
    let next = ((2.0 * ld - 1.0) * x * current - (ld - 1.0) * previous) / ld
    previous = current
    current = next
  }
  (current, previous)
}

///|
/// Integrate data evaluated at this quadrature's ordered nodes; reject overflow.
pub fn Quadrature::integrate(
  self : Quadrature,
  values : Array[Double],
) -> Result[Double, HarmonicError] {
  if values.length() != self.size() {
    return Err(ShapeMismatch(self.size(), values.length()))
  }
  let mut sum = 0.0
  let mut correction = 0.0
  for i = 0; i < values.length(); i = i + 1 {
    if !finite(values[i]) {
      return Err(InvalidInput("non-finite quadrature sample"))
    }
    let term = values[i] * self.weights[i] - correction
    let next = sum + term
    correction = next - sum - term
    sum = next
  }
  if !finite(sum) {
    Err(InvalidInput("quadrature overflow"))
  } else {
    Ok(sum)
  }
}
