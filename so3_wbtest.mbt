///|
test "SO3 rotates first-order Cartesian fields with active right-handed signs" {
  let z = Spectrum::zeros(1).unwrap().with_coefficient(1, 0, 1.0).unwrap()
  let x = z.rotate_y(@math.PI / 2.0).unwrap()
  close(x.coefficient(1, 1).unwrap(), -1.0, 3.0e-15)
  close(x.coefficient(1, 0).unwrap(), 0.0, 3.0e-15)
  let d = Direction::from_cartesian(0.0, 0.0, 1.0).unwrap()
  let (a, b, c) = d.rotate_zyz(0.0, @math.PI / 2.0, 0.0).unwrap().cartesian()
  close(a, 1.0, 2.0e-15)
  close(b, 0.0, 2.0e-15)
  close(c, 0.0, 2.0e-15)
}

///|
test "general rotation agrees with inverse Cartesian evaluation including poles" {
  let coefficients = Array::make(81, 0.0)
  for i = 0; i < coefficients.length(); i = i + 1 {
    coefficients[i] = @math.sin(i.to_double()) / (i + 1).to_double()
  }
  let s = Spectrum::new(8, coefficients).unwrap()
  for (a, b, g) in [(0.4, 0.8, -1.2), (2.0, @math.PI, 0.5), (-0.7, -1.9, 2.3)] {
    let rotated = s.rotate_zyz(a, b, g).unwrap()
    for theta in [0.0, 0.4, 1.5, @math.PI] {
      let d = Direction::new(theta, 0.9).unwrap()
      let original = d.rotate_zyz(-g, -b, -a).unwrap()
      close(
        rotated.evaluate(d).unwrap(),
        s.evaluate(original).unwrap(),
        5.0e-13,
      )
    }
    close(
      rotated.rotate_zyz(-g, -b, -a).unwrap().distance(s).unwrap(),
      0.0,
      1.0e-13,
    )
    close(rotated.energy().unwrap(), s.energy().unwrap(), 2.0e-13)
  }
}

///|
test "degree64 rotations preserve every band and compose without factorial instability" {
  let values = Array::make(4225, 0.0)
  for m = -64; m <= 64; m = m + 1 {
    values[index_of(64, m)] = @math.sin((m + 65).to_double())
  }
  let s = Spectrum::new(64, values).unwrap()
  let r = s.rotate_y(1.8).unwrap()
  close(r.rotate_y(-1.8).unwrap().distance(s).unwrap(), 0.0, 8.0e-12)
  close(r.energy().unwrap(), s.energy().unwrap(), 3.0e-11)
  close(
    s.rotate_y(0.7).unwrap().rotate_y(1.1).unwrap().distance(r).unwrap(),
    0.0,
    8.0e-12,
  )
  for l = 0; l < 64; l = l + 1 {
    for m = -l; m <= l; m = m + 1 {
      close(r.coefficient(l, m).unwrap(), 0.0, 0.0)
    }
  }
  let d = Direction::new(0.9, 1.1).unwrap()
  close(
    r.evaluate(d).unwrap(),
    s.evaluate(d.rotate_zyz(0.0, -1.8, 0.0).unwrap()).unwrap(),
    2.0e-10,
  )
}

///|
test "SO3 validates angles, preserves constants and scales large finite inputs" {
  let s = Spectrum::new(0, [2.0]).unwrap()
  close(s.rotate_zyz(0.4, 0.8, 0.3).unwrap().distance(s).unwrap(), 0.0, 0.0)
  assert_true(rejected(s.rotate_zyz(0.0, 1.0e13, 0.0)))
  assert_true(rejected(s.rotate_zyz(0.0 / 0.0, 0.0, 0.0)))
  assert_true(
    rejected(Direction::new(0.0, 0.0).unwrap().rotate_zyz(0.0, 0.0, 1.0 / 0.0)),
  )
  let big = Spectrum::new(1, [0.0, 0.0, 1.0e300, 0.0])
    .unwrap()
    .rotate_y(0.7)
    .unwrap()
  close(big.coefficient(1, 1).unwrap() / 1.0e300, -@math.sin(0.7), 3.0e-15)
  close(big.rotate_y(0.0).unwrap().distance(big).unwrap(), 0.0, 0.0)
}
