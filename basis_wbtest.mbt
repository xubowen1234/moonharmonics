///|
test "real basis agrees with Cartesian degree one and addition theorem" {
  let d = Direction::new(1.1, 2.3).unwrap()
  let b = real_basis(64, d).unwrap()
  let (x, y, z) = d.cartesian()
  let k = (3.0 / (4.0 * @math.PI)).sqrt()
  close(b[1], -k * y, 1.0e-14)
  close(b[2], k * z, 1.0e-14)
  close(b[3], -k * x, 1.0e-14)
  for l = 0; l <= 64; l = l + 1 {
    close(
      basis_band_norm(l, d).unwrap(),
      (2.0 * l.to_double() + 1.0) / (4.0 * @math.PI),
      1.0e-12,
    )
  }
  close(real_harmonic(1, -1, d).unwrap(), b[1], 1.0e-15)
  assert_true(rejected(real_harmonic(1, 2, d)))
}
