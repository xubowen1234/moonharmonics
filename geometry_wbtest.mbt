///|
fn close(a : Double, b : Double, tol : Double) -> Unit raise {
  assert_true(finite(a) && finite(b) && (a - b).abs() <= tol)
}

///|
fn[T] rejected(r : Result[T, HarmonicError]) -> Bool {
  match r {
    Err(_) => true
    Ok(_) => false
  }
}

///|
test "geometry normalizes longitudes and robust Cartesian magnitudes" {
  let d = Direction::new(@math.PI / 2.0, -@math.PI / 2.0).unwrap()
  close(d.phi(), 1.5 * @math.PI, 1.0e-14)
  let (x, y, z) = d.cartesian()
  close(x, 0.0, 1.0e-14)
  close(y, -1.0, 1.0e-14)
  close(z, 0.0, 1.0e-14)
  close(
    Direction::from_cartesian(1.0e300, 0.0, 0.0).unwrap().theta(),
    @math.PI / 2.0,
    1.0e-14,
  )
  close(
    d.distance(Direction::from_cartesian(0.0, 1.0, 0.0).unwrap()),
    @math.PI,
    1.0e-14,
  )
  assert_true(rejected(Direction::new(-0.1, 0.0)))
  assert_true(rejected(Direction::new(0.0, 1.0e20)))
  assert_true(rejected(Direction::from_cartesian(0.0, 0.0, 0.0)))
}

///|
test "Cartesian directions retain near-pole inclinations across scales" {
  // atan(eps) ~ eps: independent small-angle expectation, not a round-trip.
  for eps in [1.0e-8, 1.0e-12, 1.0e-100, 1.0e-200] {
    for scale in [1.0, 1.0e100, 1.0e-100] {
      let d = Direction::from_cartesian(eps * scale, 0.0, scale).unwrap()
      assert_true(d.theta() > 0.0)
      close(d.theta() / eps, 1.0, 1.0e-14)
      let basis = real_basis(1, d).unwrap()
      close(basis[3] / eps, -(3.0 / (4.0 * @math.PI)).sqrt(), 1.0e-14)
    }
  }
  for sign in [-1.0, 1.0] {
    let d = Direction::from_cartesian(3.0e-7, -4.0e-7, sign).unwrap()
    let expected = if sign > 0.0 { 5.0e-7 } else { @math.PI - 5.0e-7 }
    close(d.theta(), expected, 1.0e-15)
    let (x, y, _) = d.cartesian()
    close(x, 3.0e-7, 1.0e-15)
    close(y, -4.0e-7, 1.0e-15)
  }
}
