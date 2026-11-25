///|
test "DH equiangular grid area and bandlimited roundtrip" {
  for l in [0, 1, 4, 12, 32, 64] {
    let coefficients = Array::make((l + 1) * (l + 1), 0.0)
    for i = 0; i < coefficients.length(); i = i + 1 {
      coefficients[i] = @math.sin((i + 3).to_double()) / (i + 1).to_double()
    }
    let spectrum = Spectrum::new(l, coefficients).unwrap()
    let grid = DriscollHealyGrid::new(l).unwrap()
    assert_eq(grid.rings(), 2 * (l + 1))
    assert_eq(grid.longitudes(), 2 * (l + 1))
    close(
      grid.integrate(Array::make(grid.size(), 1.0)).unwrap(),
      4.0 * @math.PI,
      2.0e-13,
    )
    let recovered = grid
      .analyze(spectrum.synthesize_dh(grid).unwrap())
      .unwrap()
      .values()
    for i = 0; i < coefficients.length(); i = i + 1 {
      close(recovered[i], coefficients[i], 3.0e-12)
    }
  }
}

///|
test "DH grid rejects invalid degree shape and indices" {
  let grid = DriscollHealyGrid::new(4).unwrap()
  assert_true(rejected(grid.direction(-1)))
  assert_true(rejected(grid.direction(grid.size())))
  assert_true(rejected(grid.integrate([1.0])))
  assert_true(rejected(DriscollHealyGrid::new(-1)))
}
