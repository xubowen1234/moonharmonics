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
