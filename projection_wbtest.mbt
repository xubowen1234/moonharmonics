///|
test "streaming partial sphere integrals agree with batch and reject overflow atomically" {
  let grid = SphereGrid::new(5).unwrap()
  let a = ProjectionAccumulator::new(5).unwrap()
  let b = ProjectionAccumulator::new(5).unwrap()
  let samples = grid
    .sample(fn(d) {
      let (x, _, z) = d.cartesian()
      x * z + 1.0
    })
    .unwrap()
  for i = 0; i < grid.size(); i = i + 1 {
    let target = if i % 2 == 0 { a } else { b }
    target
    .push(grid.direction(i).unwrap(), samples[i], grid.weight(i).unwrap())
    .unwrap()
  }
  a.merge(b).unwrap()
  assert_eq(a.samples(), grid.size())
  close(a.weight(), 4.0 * @math.PI, 1.0e-13)
  close(
    a.spectrum().unwrap().distance(grid.analyze(samples).unwrap()).unwrap(),
    0.0,
    1.0e-13,
  )
  let before = a.spectrum().unwrap()
  assert_true(
    rejected(a.push(Direction::new(0.0, 0.0).unwrap(), 1.0e308, 1.0e308)),
  )
  close(a.spectrum().unwrap().distance(before).unwrap(), 0.0, 0.0)
  let bad = ProjectionAccumulator::new(2).unwrap()
  assert_true(rejected(a.merge(bad)))
  close(
    grid.inner_product(samples, samples).unwrap(),
    before.energy().unwrap(),
    1.0e-13,
  )
}
