///|
test "Gauss-Legendre integrates every monomial through exactness limit" {
  for n = 1; n <= 16; n = n + 1 {
    let q = gauss_legendre(n).unwrap()
    let nodes = q.nodes()
    for power = 0; power < 2 * n; power = power + 1 {
      let values = nodes.map(fn(x) { @math.pow(x, power.to_double()) })
      let exact = if power % 2 == 0 {
        2.0 / (power + 1).to_double()
      } else {
        0.0
      }
      close(q.integrate(values).unwrap(), exact, 2.0e-14)
    }
  }
  let q = gauss_legendre(256).unwrap()
  close(q.integrate(Array::make(256, 1.0)).unwrap(), 2.0, 2.0e-14)
  assert_true(rejected(gauss_legendre(0)))
}
