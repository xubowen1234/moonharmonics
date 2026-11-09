///|
test "Cartesian tangent gradient agrees with projected linear vector including poles" {
  let k = (4.0 * @math.PI / 3.0).sqrt()
  let f = Spectrum::new(1, [0.0, -2.0 * k, 3.0 * k, -k]).unwrap()
  for theta in [0.0, 1.0e-7, 0.4, 1.8, @math.PI - 1.0e-7, @math.PI] {
    let d = Direction::new(theta, 0.7).unwrap()
    let (x, y, z) = d.cartesian()
    let dot = x + 2.0 * y + 3.0 * z
    let (gx, gy, gz) = f.gradient(d).unwrap().components()
    close(gx, 1.0 - dot * x, 2.0e-8)
    close(gy, 2.0 - dot * y, 2.0e-8)
    close(gz, 3.0 - dot * z, 2.0e-8)
    close(gx * x + gy * y + gz * z, 0.0, 1.0e-13)
  }
}

///|
test "gradient norm integrates to spectral H1 energy" {
  let f = Spectrum::zeros(8).unwrap().with_coefficient(8, -5, 1.0).unwrap()
  let g = SphereGrid::new(10).unwrap()
  let squared = g
    .sample(fn(d) {
      let n = f.gradient(d).unwrap().norm()
      n * n
    })
    .unwrap()
  close(g.integrate(squared).unwrap(), f.gradient_energy().unwrap(), 3.0e-12)
}
