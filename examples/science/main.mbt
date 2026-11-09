// Sampled scalar sphere -> spectrum -> heat smoothing and zero-mean Poisson potential.
///|
fn main {
  let grid = @h.SphereGrid::new(8).unwrap()
  let source = grid
    .project(fn(d) {
      let (x, y, z) = d.cartesian()
      z + x * y
    })
    .unwrap()
    .without_mean()
  let potential = source.poisson(0.0).unwrap()
  let residual = potential.laplacian().unwrap().distance(source).unwrap()
  if residual > 1.0e-12 {
    abort("Poisson residual exceeded tolerance")
  }
  let d = @h.Direction::new(@math.PI / 3.0, 0.4).unwrap()
  let (x, y, z) = d.cartesian()
  let exact = -z / 2.0 - x * y / 6.0
  let actual = potential.evaluate(d).unwrap()
  if (actual - exact).abs() > 1.0e-12 {
    abort("Poisson potential differs from analytic polynomial")
  }
  let smoothed = source.heat(0.2).unwrap()
  if smoothed.energy().unwrap() >= source.energy().unwrap() {
    abort("heat energy did not decrease")
  }
  println("Poisson residual: \{residual}")
  println("potential: \{actual}; heat energy: \{smoothed.energy().unwrap()}")
}
