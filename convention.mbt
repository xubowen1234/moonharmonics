///|
/// Basis normalization: Orthonormal has unit integral-square per harmonic;
/// FourPi has integral-square 4pi (phase-free is ambisonic N3D);
/// Schmidt has integral-square 4pi/(2l+1) (phase-free is ambisonic SN3D).
pub(all) enum Normalization {
  Orthonormal
  FourPi
  Schmidt
} derive(Eq, Debug)

///|
/// Whether odd |m| carries the Condon-Shortley minus sign.
pub(all) enum HarmonicPhase {
  CondonShortley
  PhaseFree
} derive(Eq, Debug)

///|
/// Within each degree, Signed is -l,...,0,...,+l (ACN index l*l+l+m).
/// Interleaved is 0,+1,-1,+2,-2,...; CosineThenSine is 0,+1,...,+l,-1,...,-l.
/// These are layouts of complete real bands, not FuMa channel weighting/order.
pub(all) enum CoefficientOrder {
  Signed
  Interleaved
  CosineThenSine
} derive(Eq, Debug)

///|
fn convention_index(l : Int32, m : Int32, order : CoefficientOrder) -> Int32 {
  let offset = match order {
    Signed => l + m
    Interleaved => if m > 0 { 2 * m - 1 } else { -2 * m }
    CosineThenSine => if m >= 0 { m } else { l - m }
  }
  l * l + offset
}

///|
fn normalization_factor(l : Int32, normalization : Normalization) -> Double {
  match normalization {
    Orthonormal => 1.0
    FourPi => (4.0 * @math.PI).sqrt()
    Schmidt => (4.0 * @math.PI / (2 * l + 1).to_double()).sqrt()
  }
}

///|
/// Import expansion coefficients (not basis samples) into the canonical spectrum.
/// External basis B_lm = phase * normalization_factor(l) * R_lm.
/// Every array must contain exactly (degree+1)^2 finite values; input is copied.
pub fn Spectrum::from_convention(
  degree : Int32,
  coefficients : Array[Double],
  normalization : Normalization,
  phase : HarmonicPhase,
  order : CoefficientOrder,
) -> Result[Spectrum, HarmonicError] {
  let external = match Spectrum::new(degree, coefficients) {
    Ok(s) => s
    Err(e) => return Err(e)
  }
  let values = external.values()
  for l = 0; l <= degree; l = l + 1 {
    let factor = normalization_factor(l, normalization)
    for m = -l; m <= l; m = m + 1 {
      let sign = if phase == PhaseFree && m % 2 != 0 { -1.0 } else { 1.0 }
      values[index_of(l, m)] = external.coefficients[convention_index(
          l, m, order,
        )] *
        factor *
        sign
    }
  }
  Spectrum::new(degree, values)
}

///|
/// Export expansion coefficients; division reverses the external basis scale.
/// Normalization, phase and order are independent. Overflow returns an error.
pub fn Spectrum::to_convention(
  self : Spectrum,
  normalization : Normalization,
  phase : HarmonicPhase,
  order : CoefficientOrder,
) -> Result[Array[Double], HarmonicError] {
  let values = self.values()
  for l = 0; l <= self.degree; l = l + 1 {
    let factor = normalization_factor(l, normalization)
    for m = -l; m <= l; m = m + 1 {
      let sign = if phase == PhaseFree && m % 2 != 0 { -1.0 } else { 1.0 }
      let value = self.coefficients[index_of(l, m)] / factor * sign
      if !finite(value) {
        return Err(InvalidInput("convention conversion overflow"))
      }
      values[convention_index(l, m, order)] = value
    }
  }
  Ok(values)
}

///|
pub extend Normalization with Eq::{equal, not_equal}

///|
pub extend HarmonicPhase with Eq::{equal, not_equal}

///|
pub extend CoefficientOrder with Eq::{equal, not_equal}

///|
pub extend Normalization with @debug.Debug::{to_repr}

///|
pub extend HarmonicPhase with @debug.Debug::{to_repr}

///|
pub extend CoefficientOrder with @debug.Debug::{to_repr}
