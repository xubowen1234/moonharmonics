// Offline directional scene coefficients -> active z rotation -> evaluated directional field.
// This demonstrates a finite angular expansion, not an audio device or HRTF renderer.

///|
fn main {
  let front = @h.Direction::new(@math.PI / 2.0, 0.0).unwrap()
  let side = @h.Direction::new(@math.PI / 2.0, @math.PI / 2.0).unwrap()
  let scene = @h.encode_direction(3, front, 1.0)
    .unwrap()
    .add(@h.encode_direction(3, side, 0.5).unwrap())
    .unwrap()
  let rotated = scene.rotate_z(@math.PI / 2.0).unwrap()
  let before = scene.evaluate(front).unwrap()
  let after = rotated.evaluate(side).unwrap()
  if (before - after).abs() > 1.0e-12 {
    abort("active directional rotation mismatch")
  }
  let interoperable = rotated.to_sn3d().unwrap()
  let restored = @h.Spectrum::from_sn3d(3, interoperable).unwrap()
  if restored.distance(rotated).unwrap() > 1.0e-12 {
    abort("SN3D convention round trip mismatch")
  }
  println("rotated directional value: \{after}")
  println("orthonormal coefficients: \{rotated.to_json_string()}")
}
