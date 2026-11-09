///|
test "sphere grids preserve area and Cartesian moments" {
  let g = SphereGrid::new(8).unwrap()
  close(
    g.integrate(Array::make(g.size(), 1.0)).unwrap(),
    4.0 * @math.PI,
    1.0e-13,
  )
  let z2 = g
    .sample(fn(d) {
      let (_, _, z) = d.cartesian()
      z * z
    })
    .unwrap()
  close(g.integrate(z2).unwrap(), 4.0 * @math.PI / 3.0, 1.0e-13)
  assert_true(rejected(g.direction(-1)))
  assert_true(rejected(SphereGrid::oversampled(5, 5, 11)))
  assert_true(rejected(SphereGrid::oversampled(5, 6, 10)))
  assert_true(rejected(g.integrate([1.0])))
}
