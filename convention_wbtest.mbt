///|
test "all18 convention combinations roundtrip degree64 and preserve expansion values" {
  let values = Array::make(4225, 0.0)
  for i = 0; i < values.length(); i = i + 1 {
    values[i] = @math.cos((i + 1).to_double())
  }
  let spectrum = Spectrum::new(64, values).unwrap()
  for normalization in [Orthonormal, FourPi, Schmidt] {
    for phase in [CondonShortley, PhaseFree] {
      for order in [Signed, Interleaved, CosineThenSine] {
        let external = spectrum
          .to_convention(normalization, phase, order)
          .unwrap()
        let imported = Spectrum::from_convention(
          64, external, normalization, phase, order,
        ).unwrap()
        close(imported.distance(spectrum).unwrap(), 0.0, 2.0e-14)
        external[0] = 42.0
        close(spectrum.coefficient(0, 0).unwrap(), @math.cos(1.0), 0.0)
      }
    }
  }
}

///|
test "known first and second bands validate external order, phase and scaling" {
  // Tagged canonical coefficients: l0=1; l1=[2,3,4]; l2=[5,6,7,8,9].
  let s = Spectrum::new(2, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]).unwrap()
  let interleaved = s
    .to_convention(Orthonormal, CondonShortley, Interleaved)
    .unwrap()
  let block = s.to_convention(Orthonormal, PhaseFree, CosineThenSine).unwrap()
  let expected_a = [1.0, 3.0, 4.0, 2.0, 7.0, 8.0, 6.0, 9.0, 5.0]
  let expected_b = [1.0, 3.0, -4.0, -2.0, 7.0, -8.0, 9.0, -6.0, 5.0]
  for i = 0; i < 9; i = i + 1 {
    close(interleaved[i], expected_a[i], 0.0)
    close(block[i], expected_b[i], 0.0)
  }
  let mono = Spectrum::from_convention(0, [1.0], FourPi, PhaseFree, Signed).unwrap()
  close(mono.evaluate(Direction::new(0.2, 0.4).unwrap()).unwrap(), 1.0, 2.0e-15)
  let acn_sn3d = Spectrum::from_sn3d(1, [1.0, 2.0, 3.0, 4.0]).unwrap()
  let n3d = acn_sn3d.to_convention(FourPi, PhaseFree, Signed).unwrap()
  close(n3d[0], 1.0, 2.0e-15)
  close(n3d[1], 2.0 / 3.0.sqrt(), 2.0e-15)
  close(n3d[2], 3.0 / 3.0.sqrt(), 2.0e-15)
  close(n3d[3], 4.0 / 3.0.sqrt(), 2.0e-15)
}

///|
test "convention conversion rejects invalid inputs and both overflow directions" {
  assert_true(
    rejected(Spectrum::from_convention(1, [0.0], Schmidt, PhaseFree, Signed)),
  )
  assert_true(
    rejected(
      Spectrum::from_convention(-1, [], Orthonormal, CondonShortley, Signed),
    ),
  )
  assert_true(
    rejected(
      Spectrum::from_convention(
        0,
        [1.0 / 0.0],
        Orthonormal,
        CondonShortley,
        Signed,
      ),
    ),
  )
  assert_true(
    rejected(Spectrum::from_convention(0, [1.0e308], FourPi, PhaseFree, Signed)),
  )
  let huge = Spectrum::zeros(64)
    .unwrap()
    .with_coefficient(64, 0, 1.0e308)
    .unwrap()
  assert_true(rejected(huge.to_convention(Schmidt, CondonShortley, Signed)))
}
