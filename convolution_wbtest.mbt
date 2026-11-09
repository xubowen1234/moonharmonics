///|
test "Lambertian multipliers integrate a hemisphere and retain high even bands" {
  let kernel = ZonalKernel::lambertian(64).unwrap()
  let gains = kernel.gains()
  close(gains[0], @math.PI, 1.0e-14)
  close(gains[1], 2.0 * @math.PI / 3.0, 1.0e-14)
  close(gains[2], @math.PI / 4.0, 1.0e-14)
  close(gains[4], -@math.PI / 24.0, 1.0e-14)
  let q = gauss_legendre(40).unwrap()
  for l = 0; l <= 64; l = l + 1 {
    let mut sum = 0.0
    for i = 0; i < q.size(); i = i + 1 {
      let t = (q.nodes[i] + 1.0) / 2.0
      let (p, _) = legendre_pair(l, t)
      sum = sum + @math.PI * q.weights[i] * t * p
    }
    close(gains[l], sum, 2.0e-14)
  }
  let constant = SphereGrid::new(2).unwrap().project(fn(_) { 1.0 }).unwrap()
  close(
    constant.convolve(ZonalKernel::lambertian(2).unwrap()).unwrap().mean(),
    @math.PI,
    1.0e-14,
  )
  let smooth = ZonalKernel::project(3, 16, fn(t) { t }).unwrap().gains()
  close(smooth[1], 4.0 * @math.PI / 3.0, 1.0e-14)
  close(smooth[0], 0.0, 1.0e-14)
}
