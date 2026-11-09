///|
test "heat has semigroup law, preserves means and contracts energy" {
  let s = Spectrum::new(2, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]).unwrap()
  let a = s.heat(0.2).unwrap().heat(0.3).unwrap()
  close(a.distance(s.heat(0.5).unwrap()).unwrap(), 0.0, 1.0e-14)
  close(a.mean(), s.mean(), 0.0)
  assert_true(a.energy().unwrap() < s.energy().unwrap())
  close(s.heat(1.0e308).unwrap().coefficient(0, 0).unwrap(), 1.0, 0.0)
  close(s.cosine_taper(0, 2).unwrap().coefficient(1, 0).unwrap(), 1.5, 1.0e-14)
  close(s.bandpass(1, 1).unwrap().coefficient(2, 0).unwrap(), 0.0, 0.0)
  assert_true(rejected(s.heat(-1.0)))
  assert_true(rejected(s.cosine_taper(1, 1)))
}
