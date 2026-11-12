///|
fn rotation_angle(angle : Double) -> Result[Double, HarmonicError] {
  if !finite(angle) || angle.abs() > 1.0e12 {
    return Err(
      InvalidInput("rotation angle must be finite and bounded by 1.0e12"),
    )
  }
  let a = Direction::new(0.0, angle).unwrap().phi
  Ok(if a > @math.PI { a - 2.0 * @math.PI } else { a })
}

///|
/// Active right-handed Rz(alpha) Ry(beta) Rz(gamma), applied right to left.
/// Angles are radians. The inverse is (-gamma,-beta,-alpha).
pub fn Direction::rotate_zyz(
  self : Direction,
  alpha : Double,
  beta : Double,
  gamma : Double,
) -> Result[Direction, HarmonicError] {
  let a = match rotation_angle(alpha) {
    Ok(a) => a
    Err(e) => return Err(e)
  }
  let b = match rotation_angle(beta) {
    Ok(b) => b
    Err(e) => return Err(e)
  }
  let g = match rotation_angle(gamma) {
    Ok(g) => g
    Err(e) => return Err(e)
  }
  let (x, y, z) = self.cartesian()
  let u = @math.cos(g) * x - @math.sin(g) * y
  let v = @math.sin(g) * x + @math.cos(g) * y
  let w = @math.cos(b) * u + @math.sin(b) * z
  Direction::from_cartesian(
    @math.cos(a) * w - @math.sin(a) * v,
    @math.sin(a) * w + @math.cos(a) * v,
    -@math.sin(b) * u + @math.cos(b) * z,
  )
}

///|
/// Active y rotation through exp(beta G_l), evaluated by bounded Taylor steps.
/// Each step has ||h G_l||_infinity < 1/2. Eighteen terms avoid factorial
/// cancellation in high-degree Wigner sums. Work O(L^3 |beta|), storage O(L^2).
pub fn Spectrum::rotate_y(
  self : Spectrum,
  angle : Double,
) -> Result[Spectrum, HarmonicError] {
  let beta = match rotation_angle(angle) {
    Ok(a) => a
    Err(e) => return Err(e)
  }
  if beta == 0.0 {
    return Ok(self)
  }
  let values = self.values()
  let root_two = 2.0.sqrt()
  for l = 1; l <= self.degree; l = l + 1 {
    let n = 2 * l + 1
    let mut scale = 0.0
    for m = -l; m <= l; m = m + 1 {
      scale = @cmp.maximum(scale, values[index_of(l, m)].abs())
    }
    if scale == 0.0 {
      continue
    }
    // Complex expansion a_m=(c_m-i*c_-m)/sqrt(2); a_-m=(-1)^m conj(a_m).
    let re = Array::make(n, 0.0)
    let im = Array::make(n, 0.0)
    re[l] = values[index_of(l, 0)] / scale
    for m = 1; m <= l; m = m + 1 {
      let sign = if m % 2 == 0 { 1.0 } else { -1.0 }
      re[l + m] = values[index_of(l, m)] / scale / root_two
      im[l + m] = -values[index_of(l, -m)] / scale / root_two
      re[l - m] = sign * re[l + m]
      im[l - m] = -sign * im[l + m]
    }
    // G[m,m+1]=sqrt((l-m)(l+m+1))/2, G[m,m-1]=-G[m-1,m].
    let edge = Array::make(n - 1, 0.0)
    for i = 0; i < n - 1; i = i + 1 {
      let m = i - l
      edge[i] = 0.5 * ((l - m) * (l + m + 1)).to_double().sqrt()
    }
    let steps = @cmp.maximum(
      1,
      (2.0 * (l + 1).to_double() * beta.abs()).ceil().to_int(),
    )
    let h = beta / steps.to_double()
    let term_re = Array::make(n, 0.0)
    let term_im = Array::make(n, 0.0)
    let next_re = Array::make(n, 0.0)
    let next_im = Array::make(n, 0.0)
    for step = 0; step < steps; step = step + 1 {
      for i = 0; i < n; i = i + 1 {
        term_re[i] = re[i]
        term_im[i] = im[i]
      }
      for k = 1; k <= 18; k = k + 1 {
        let factor = h / k.to_double()
        for i = 0; i < n; i = i + 1 {
          let upper_re = if i + 1 < n { edge[i] * term_re[i + 1] } else { 0.0 }
          let lower_re = if i > 0 { edge[i - 1] * term_re[i - 1] } else { 0.0 }
          let upper_im = if i + 1 < n { edge[i] * term_im[i + 1] } else { 0.0 }
          let lower_im = if i > 0 { edge[i - 1] * term_im[i - 1] } else { 0.0 }
          next_re[i] = factor * (upper_re - lower_re)
          next_im[i] = factor * (upper_im - lower_im)
        }
        for i = 0; i < n; i = i + 1 {
          term_re[i] = next_re[i]
          term_im[i] = next_im[i]
          re[i] = re[i] + term_re[i]
          im[i] = im[i] + term_im[i]
        }
      }
    }
    values[index_of(l, 0)] = re[l] * scale
    for m = 1; m <= l; m = m + 1 {
      values[index_of(l, m)] = re[l + m] * root_two * scale
      values[index_of(l, -m)] = -im[l + m] * root_two * scale
    }
  }
  Spectrum::new(self.degree, values)
}

///|
/// General active SO(3) rotation Rz(alpha) Ry(beta) Rz(gamma).
/// The resulting field is g(d)=f(R^-1 d); no resampling or order truncation.
pub fn Spectrum::rotate_zyz(
  self : Spectrum,
  alpha : Double,
  beta : Double,
  gamma : Double,
) -> Result[Spectrum, HarmonicError] {
  // Validate every angle before any arithmetic, including zero spectra.
  match rotation_angle(alpha) {
    Ok(_) => ()
    Err(e) => return Err(e)
  }
  match rotation_angle(beta) {
    Ok(_) => ()
    Err(e) => return Err(e)
  }
  match rotation_angle(gamma) {
    Ok(_) => ()
    Err(e) => return Err(e)
  }
  let first = match self.rotate_z(gamma) {
    Ok(s) => s
    Err(e) => return Err(e)
  }
  let second = match first.rotate_y(beta) {
    Ok(s) => s
    Err(e) => return Err(e)
  }
  second.rotate_z(alpha)
}
