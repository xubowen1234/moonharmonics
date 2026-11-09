///|
test "pivoted QR recovers independent analytic field under unequal weights" {
  let grid = SphereGrid::new(4).unwrap()
  let dirs : Array[Direction] = []
  let values : Array[Double] = []
  let weights : Array[Double] = []
  for i = 0; i < grid.size(); i = i + 1 {
    let d = grid.direction(i).unwrap()
    let (x, y, z) = d.cartesian()
    dirs.push(d)
    values.push(2.0 + x - 3.0 * y + z)
    weights.push((i % 7 + 1).to_double())
  }
  let fit = fit_weighted(2, dirs, values, weights).unwrap()
  let k = (4.0 * @math.PI / 3.0).sqrt()
  close(
    fit.spectrum().coefficient(0, 0).unwrap(),
    2.0 * (4.0 * @math.PI).sqrt(),
    1.0e-13,
  )
  close(fit.spectrum().coefficient(1, -1).unwrap(), 3.0 * k, 1.0e-13)
  close(fit.spectrum().coefficient(1, 1).unwrap(), -k, 1.0e-13)
  close(fit.weighted_rms(), 0.0, 1.0e-13)
  assert_eq(fit.rank(), 9)
  let scaled = weights.map(fn(w) { w * 1.0e250 })
  close(
    fit
    .spectrum()
    .distance(fit_weighted(2, dirs, values, scaled).unwrap().spectrum())
    .unwrap(),
    0.0,
    1.0e-12,
  )
  let identical = Array::make(9, dirs[0])
  assert_true(
    rejected(
      fit_weighted(2, identical, Array::make(9, 1.0), Array::make(9, 1.0)),
    ),
  )
  assert_true(
    rejected(fit_weighted(2, dirs, values, Array::make(grid.size(), 0.0))),
  )
  weights[0] = -1.0
  assert_true(rejected(fit_weighted(2, dirs, values, weights)))
}
