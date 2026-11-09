///|
test "Poisson sign, gauge and phi derivative are independent analytic contracts" {
  let s = Spectrum::zeros(4)
    .unwrap()
    .with_coefficient(1, 0, 1.0)
    .unwrap()
    .with_coefficient(4, 3, 2.0)
    .unwrap()
  close(
    s.poisson(0.0).unwrap().laplacian().unwrap().distance(s).unwrap(),
    0.0,
    1.0e-14,
  )
  close(s.laplacian().unwrap().coefficient(1, 0).unwrap(), -2.0, 0.0)
  close(s.differentiate_phi().unwrap().coefficient(4, -3).unwrap(), -6.0, 0.0)
  assert_true(rejected(s.with_coefficient(0, 0, 1.0).unwrap().poisson(1.0e-14)))
  close(s.screened_poisson(0.5).unwrap().coefficient(1, 0).unwrap(), 0.5, 0.0)
}
