// Directional radiance -> harmonic projection -> cosine convolution -> irradiance.

///|
fn main {
  let grid = @h.SphereGrid::new(6).unwrap()
  let radiance = grid
    .project(fn(d) {
      let (_, _, z) = d.cartesian()
      2.0 + z
    })
    .unwrap()
  let irradiance = radiance
    .convolve(@h.ZonalKernel::lambertian(6).unwrap())
    .unwrap()
  let north = irradiance.evaluate(@h.Direction::new(0.0, 0.0).unwrap()).unwrap()
  let exact = 2.0 * @math.PI + 2.0 * @math.PI / 3.0
  if (north - exact).abs() > 1.0e-12 {
    abort("irradiance did not match analytic hemisphere integral")
  }
  println("north irradiance: \{north}")
  let reflected = irradiance.scale(0.7 / @math.PI).unwrap()
  println(
    "reflected radiance (albedo 0.7): \{reflected.evaluate(@h.Direction::new(0.0,0.0).unwrap()).unwrap()}",
  )
}
