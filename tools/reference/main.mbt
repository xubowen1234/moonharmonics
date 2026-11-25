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
  let dh_grid = @h.DriscollHealyGrid::new(16).unwrap()
  let dh_synthesis = spectrum.synthesize_dh(dh_grid).unwrap()
  emit("dh_synthesis", numbers(dh_synthesis))
  let dh_projection = dh_grid
    .project(fn(d) {
      let (x, y, z) = d.cartesian()
      @math.exp(0.5 * z) + x * y + 0.2 * x
    })
    .unwrap()
  emit("dh_analysis", numbers(dh_projection.values()))
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

  let high_values = Array::make(4225, 0.0)
  for i = 0; i < high_values.length(); i = i + 1 {
    high_values[i] = @math.sin((i + 1).to_double()) / (i + 1).to_double()
  }
  let high_spectrum = @h.Spectrum::new(64, high_values).unwrap()
  let source = @h.Direction::new(0.8, 1.3).unwrap()
  let encoded = @h.encode_direction(64, source, 1.0).unwrap()
  for
    (a, b, g) in [
      (0.4, 0.8, -1.2),
      (-0.7, -1.9, 2.3),
      (2.0, @math.PI, 0.5),
      (0.0, 1.0e-10, 0.0),
    ] {
    let d = source.rotate_zyz(a, b, g).unwrap()
    emit(
      "so3_encoding",
      Json::array([
        numbers([a, b, g, source.theta(), source.phi()]),
        numbers(encoded.rotate_zyz(a, b, g).unwrap().values()),
        numbers([d.theta(), d.phi()]),
      ]),
    )
  }
  for (a, b, g) in [(0.4, 0.8, -1.2), (-0.7, -1.9, 2.3), (2.0, @math.PI, 0.5)] {
    let rotated = high_spectrum.rotate_zyz(a, b, g).unwrap()
    let samples : Array[Json] = []
    for theta in [0.0, 1.0e-5, 0.4, 1.2, 2.8, @math.PI - 1.0e-5, @math.PI] {
      for phi in [0.0, 0.7, 5.2] {
        let d = @h.Direction::new(theta, phi).unwrap()
        samples.push(numbers([theta, phi, rotated.evaluate(d).unwrap()]))
      }
    }
    emit(
      "so3_field",
      Json::array([
        numbers([a, b, g]),
        numbers(rotated.values()),
        Json::array(samples),
      ]),
    )
  }
  let convention_values = Array::make(81, 0.0)
  let external_values = Array::make(81, 0.0)
  for i = 0; i < 81; i = i + 1 {
    convention_values[i] = @math.sin((i + 1).to_double()) / (i + 1).to_double()
    external_values[i] = @math.cos((i + 1).to_double()) / (i + 1).to_double()
  }
  let canonical = @h.Spectrum::new(8, convention_values).unwrap()
  let normalizations : Array[@h.Normalization] = [
    @h.Orthonormal,
    @h.FourPi,
    @h.Schmidt,
  ]
  let phases : Array[@h.HarmonicPhase] = [@h.CondonShortley, @h.PhaseFree]
  let orders : Array[@h.CoefficientOrder] = [
    @h.Signed,
    @h.Interleaved,
    @h.CosineThenSine,
  ]
  for ni = 0; ni < 3; ni = ni + 1 {
    for pi = 0; pi < 2; pi = pi + 1 {
      for oi = 0; oi < 3; oi = oi + 1 {
        let n = normalizations[ni]
        let p = phases[pi]
        let o = orders[oi]
        emit(
          "convention",
          Json::array([
            numbers([ni.to_double(), pi.to_double(), oi.to_double()]),
            numbers(canonical.to_convention(n, p, o).unwrap()),
            numbers(
              @h.Spectrum::from_convention(8, external_values, n, p, o)
              .unwrap()
              .values(),
            ),
          ]),
        )
      }
    }
  }
}
