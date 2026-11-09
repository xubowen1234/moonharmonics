///|
test "finite arithmetic errors do not leak infinities through conversion or integration" {
  let huge = Spectrum::zeros(64)
    .unwrap()
    .with_coefficient(64, 0, 1.0e308)
    .unwrap()
  assert_true(rejected(huge.to_sn3d()))
  assert_true(rejected(huge.laplacian()))
  let nan = 0.0 / 0.0
  assert_true(rejected(normalized_legendre(1, nan)))
  assert_true(rejected(Direction::new(nan, 0.0)))
  assert_true(rejected(Spectrum::new(0, [nan])))
  let grid = SphereGrid::new(2).unwrap()
  assert_true(rejected(grid.sample(fn(_) { nan })))
  assert_true(rejected(grid.analyze(Array::make(grid.size(), nan))))
  let linear = Spectrum::zeros(1)
    .unwrap()
    .with_coefficient(1, 1, -(4.0 * @math.PI / 3.0).sqrt())
    .unwrap()
  let (gx, gy, gz) = linear
    .gradient(Direction::new(1.0e-310, 0.7).unwrap())
    .unwrap()
    .components()
  close(gx, 1.0, 1.0e-14)
  close(gy, 0.0, 1.0e-14)
  close(gz, 0.0, 1.0e-14)
}

///|
test "antipodal parity holds for every signed mode through degree 64" {
  let a = real_basis(64, Direction::new(0.47, 0.63).unwrap()).unwrap()
  let b = real_basis(
    64,
    Direction::new(@math.PI - 0.47, 0.63 + @math.PI).unwrap(),
  ).unwrap()
  for l = 0; l <= 64; l = l + 1 {
    let sign = if l % 2 == 0 { 1.0 } else { -1.0 }
    for m = -l; m <= l; m = m + 1 {
      let i = index_of(l, m)
      close(b[i], sign * a[i], 2.0e-13)
    }
  }
}
