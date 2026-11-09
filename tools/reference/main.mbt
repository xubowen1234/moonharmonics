///|
fn numbers(values : Array[Double]) -> Json {
  Json::array(values.map(fn(v) { Json::number(v) }))
}

///|
fn emit(kind : String, value : Json) -> Unit {
  println(
    Json::object({ "kind": Json::string(kind), "value": value }).stringify(),
  )
}

///|
fn main {
  for theta in [0.0, 1.0e-7, 0.4, 1.2, 2.8, @math.PI - 1.0e-7, @math.PI] {
    for phi in [0.0, 0.7, 5.2] {
      let basis = @h.real_basis(64, @h.Direction::new(theta, phi).unwrap()).unwrap()
      emit(
        "basis",
        Json::array([Json::number(theta), Json::number(phi), numbers(basis)]),
      )
    }
  }
  let q = @h.gauss_legendre(256).unwrap()
  emit("quadrature", Json::array([numbers(q.nodes()), numbers(q.weights())]))
  let grid = @h.SphereGrid::new(16).unwrap()
  let coefficients = Array::make(289, 0.0)
  for i = 0; i < coefficients.length(); i = i + 1 {
    coefficients[i] = @math.sin((i + 1).to_double()) / (i + 1).to_double()
  }
  let spectrum = @h.Spectrum::new(16, coefficients).unwrap()
  emit("synthesis", numbers(spectrum.synthesize(grid).unwrap()))
  let projection = grid
    .project(fn(d) {
      let (x, y, z) = d.cartesian()
      @math.exp(0.5 * z) + x * y + 0.2 * x
    })
    .unwrap()
  emit("analysis", numbers(projection.values()))
  for theta in [0.2, 1.1, 2.8] {
    let d = @h.Direction::new(theta, 0.7).unwrap()
    let (x, y, z) = spectrum.gradient(d).unwrap().components()
    emit("gradient", numbers([theta, 0.7, x, y, z]))
  }
  emit("lambertian", numbers(@h.ZonalKernel::lambertian(64).unwrap().gains()))
  let fit_grid = @h.SphereGrid::new(5).unwrap()
  let directions : Array[@h.Direction] = []
  let values : Array[Double] = []
  let weights : Array[Double] = []
  for i = 0; i < fit_grid.size(); i = i + 1 {
    let d = fit_grid.direction(i).unwrap()
    let (x, y, z) = d.cartesian()
    directions.push(d)
    values.push(2.0 + x - 3.0 * y + z + 0.01 * @math.sin(3.0 * i.to_double()))
    weights.push((i % 7 + 1).to_double())
  }
  let fit = @h.fit_weighted(2, directions, values, weights).unwrap()
  emit(
    "fit",
    Json::array([
      numbers(fit.spectrum().values()),
      Json::number(fit.weighted_rms()),
    ]),
  )
}
