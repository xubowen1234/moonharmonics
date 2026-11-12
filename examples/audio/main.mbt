// Offline directional scene -> arbitrary active SO3 rotation -> N3D exchange.
// This demonstrates a finite angular expansion, not an audio device or HRTF renderer.

///|
fn main {
  let front = @h.Direction::new(@math.PI / 2.0, 0.0).unwrap()
  let side = @h.Direction::new(@math.PI / 2.0, @math.PI / 2.0).unwrap()
  let scene = @h.encode_direction(3, front, 1.0)
    .unwrap()
    .add(@h.encode_direction(3, side, 0.5).unwrap())
    .unwrap()
  // Demonstrate a complete external SN3D/ACN -> canonical -> rotation pipeline.
  let input = scene.to_sn3d().unwrap()
  let imported = @h.Spectrum::from_convention(
    3,
    input,
    @h.Schmidt,
    @h.PhaseFree,
    @h.Signed,
  ).unwrap()
  let rotated = imported.rotate_zyz(0.4, 0.8, -1.2).unwrap()
  let destination = front.rotate_zyz(0.4, 0.8, -1.2).unwrap()
  let before = scene.evaluate(front).unwrap()
  let after = rotated.evaluate(destination).unwrap()
  if (before - after).abs() > 1.0e-12 {
    abort("active SO3 directional rotation mismatch")
  }
  // FourPi + PhaseFree is N3D; the order is explicitly changed for this example.
  let interoperable = rotated
    .to_convention(@h.FourPi, @h.PhaseFree, @h.Interleaved)
    .unwrap()
  let restored = @h.Spectrum::from_convention(
    3,
    interoperable,
    @h.FourPi,
    @h.PhaseFree,
    @h.Interleaved,
  ).unwrap()
  if restored.distance(rotated).unwrap() > 1.0e-12 {
    abort("N3D convention round trip mismatch")
  }
  println("SO3 rotated directional value: \{after}")
  println("orthonormal coefficients: \{rotated.to_json_string()}")
}
