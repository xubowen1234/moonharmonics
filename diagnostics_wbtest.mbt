///|
test "Parseval and energy statistics use the orthonormal convention" {
  let s = Spectrum::new(1, [(4.0 * @math.PI).sqrt(), 2.0, 3.0, 4.0]).unwrap()
  close(s.mean(), 1.0, 1.0e-15)
  close(s.energy().unwrap(), 4.0 * @math.PI + 29.0, 1.0e-14)
  close(s.gradient_energy().unwrap(), 58.0, 1.0e-14)
  let g = SphereGrid::new(1).unwrap()
  close(s.parseval_defect(g, s.synthesize(g).unwrap()).unwrap(), 0.0, 1.0e-14)
  close(s.distance(s).unwrap(), 0.0, 0.0)
  let huge = Spectrum::new(0, [1.0e200]).unwrap()
  close(
    huge.distance(Spectrum::zeros(0).unwrap()).unwrap() / 1.0e200,
    1.0,
    1.0e-15,
  )
  assert_true(rejected(huge.energy()))
}
