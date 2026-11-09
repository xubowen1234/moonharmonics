///|
test "normalized recurrence has analytic signs and polar zeros" {
  let x = 0.3
  let p = normalized_legendre(64, x).unwrap()
  close(p[index_of(1, 0)], (3.0 / (4.0 * @math.PI)).sqrt() * x, 1.0e-15)
  close(
    p[index_of(1, 1)],
    -(3.0 / (8.0 * @math.PI)).sqrt() * (1.0 - x * x).sqrt(),
    1.0e-15,
  )
  close(
    p[index_of(2, 0)],
    (5.0 / (16.0 * @math.PI)).sqrt() * (3.0 * x * x - 1.0),
    1.0e-15,
  )
  let north = normalized_legendre(64, 1.0).unwrap()
  for l = 0; l <= 64; l = l + 1 {
    close(
      north[index_of(l, 0)],
      ((2.0 * l.to_double() + 1.0) / (4.0 * @math.PI)).sqrt(),
      2.0e-13,
    )
    for m = 1; m <= l; m = m + 1 {
      close(north[index_of(l, m)], 0.0, 0.0)
    }
  }
  assert_true(rejected(normalized_legendre(65, 0.0)))
  assert_true(rejected(normalized_legendre(2, 1.01)))
}
