///|
test "versioned JSON is reversible and rejects ambiguous inputs" {
  let s = Spectrum::new(1, [1.0, 2.0, 3.0, 4.0]).unwrap()
  close(
    Spectrum::from_json_string(s.to_json_string()).unwrap().distance(s).unwrap(),
    0.0,
    0.0,
  )
  assert_true(
    rejected(Spectrum::from_json_string("{\"degree\":1,\"degree\":0}")),
  )
  assert_true(
    rejected(
      Spectrum::from_json_string(
        s.to_json_string().replace(old="orthonormal-cs-l*l+l+m", new="SN3D"),
      ),
    ),
  )
  assert_true(rejected(Spectrum::from_json_string("[]")))
  assert_true(
    rejected(
      Spectrum::from_json_string(
        s.to_json_string().replace(old="\"degree\":1", new="\"degree\":1.5"),
      ),
    ),
  )
  close(
    Spectrum::from_sn3d(1, s.to_sn3d().unwrap()).unwrap().distance(s).unwrap(),
    0.0,
    2.0e-15,
  )
  let mono = Spectrum::from_sn3d(0, [1.0]).unwrap()
  close(mono.evaluate(Direction::new(0.2, 0.3).unwrap()).unwrap(), 1.0, 1.0e-15)
}
