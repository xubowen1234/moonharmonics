///|
test "spectra own storage and validate shape and updates" {
  let input = [1.0, 2.0, 3.0, 4.0]
  let s = Spectrum::new(1, input).unwrap()
  input[0] = 99.0
  let output = s.values()
  output[0] = 88.0
  close(s.coefficient(0, 0).unwrap(), 1.0, 0.0)
  close(
    s.with_coefficient(1, -1, 8.0).unwrap().coefficient(1, -1).unwrap(),
    8.0,
    0.0,
  )
  close(s.coefficient(1, -1).unwrap(), 2.0, 0.0)
  assert_eq(s.resize(0).unwrap().size(), 1)
  close(s.resize(3).unwrap().coefficient(3, 0).unwrap(), 0.0, 0.0)
  assert_true(rejected(Spectrum::new(2, input)))
  assert_true(rejected(s.add(Spectrum::zeros(2).unwrap())))
  assert_true(rejected(s.scale(1.0 / 0.0)))
  close(s.add(s).unwrap().coefficient(1, 1).unwrap(), 8.0, 0.0)
}
