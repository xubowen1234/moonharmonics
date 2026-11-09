///|
test "separable grid synthesis equals analytic linear field" {
  let k = (4.0 * @math.PI / 3.0).sqrt()
  let s = Spectrum::new(1, [0.0, -2.0 * k, 3.0 * k, -k]).unwrap()
  let g = SphereGrid::new(3).unwrap()
  let values = s.synthesize(g).unwrap()
  for i = 0; i < values.length(); i = i + 1 {
    let d = g.direction(i).unwrap()
    let (x, y, z) = d.cartesian()
    close(values[i], x + 2.0 * y + 3.0 * z, 2.0e-14)
    close(s.evaluate(d).unwrap(), values[i], 2.0e-14)
  }
  assert_true(rejected(s.synthesize(SphereGrid::new(0).unwrap())))
}
