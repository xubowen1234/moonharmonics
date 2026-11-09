///|
test "active z rotation agrees with independent shifted evaluation" {
  let s = Spectrum::zeros(5)
    .unwrap()
    .with_coefficient(5, 3, 2.0)
    .unwrap()
    .with_coefficient(2, -1, 3.0)
    .unwrap()
  let r = s.rotate_z(0.7).unwrap()
  let d = Direction::new(0.9, 1.4).unwrap()
  close(
    r.evaluate(d).unwrap(),
    s.evaluate(Direction::new(0.9, 0.7).unwrap()).unwrap(),
    2.0e-14,
  )
  close(r.rotate_z(-0.7).unwrap().distance(s).unwrap(), 0.0, 1.0e-14)
  close(r.energy().unwrap(), s.energy().unwrap(), 1.0e-13)
  let c = Spectrum::zeros(1).unwrap().with_coefficient(1, 1, 1.0).unwrap()
  close(
    c.rotate_z(@math.PI / 2.0).unwrap().coefficient(1, -1).unwrap(),
    1.0,
    1.0e-14,
  )
  let encoded = encode_direction(3, d, 2.0).unwrap()
  close(
    encoded.coefficient(0, 0).unwrap(),
    2.0 / (4.0 * @math.PI).sqrt(),
    1.0e-15,
  )
}
