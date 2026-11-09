///|
test "separable projection recovers all coefficients at representative band limits" {
  for l in [0, 1, 4, 12, 32, 64] {
    let values = Array::make((l + 1) * (l + 1), 0.0)
    for i = 0; i < values.length(); i = i + 1 {
      values[i] = @math.sin((i + 1).to_double()) / (i + 1).to_double()
    }
    let spectrum = Spectrum::new(l, values).unwrap()
    let grid = SphereGrid::new(l).unwrap()
    let recovered = grid
      .analyze(spectrum.synthesize(grid).unwrap())
      .unwrap()
      .values()
    for i = 0; i < values.length(); i = i + 1 {
      close(recovered[i], values[i], 2.0e-13)
    }
  }
}

///|
test "projection has correct known constant and dipole coefficients" {
  let s = SphereGrid::new(5)
    .unwrap()
    .project(fn(d) {
      let (_, _, z) = d.cartesian()
      2.0 + z
    })
    .unwrap()
  close(s.coefficient(0, 0).unwrap(), 2.0 * (4.0 * @math.PI).sqrt(), 1.0e-13)
  close(s.coefficient(1, 0).unwrap(), (4.0 * @math.PI / 3.0).sqrt(), 1.0e-13)
  close(s.coefficient(3, 0).unwrap(), 0.0, 1.0e-13)
}
