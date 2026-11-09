///|
/// Result of a full-rank weighted least-squares solve using column-pivoted Householder QR.
pub struct FitResult {
  spectrum : Spectrum
  weighted_rms : Double
  rank : Int32
  diagonal_ratio : Double
}

///|
pub fn FitResult::spectrum(self : FitResult) -> Spectrum {
  self.spectrum
}

///|
pub fn FitResult::weighted_rms(self : FitResult) -> Double {
  self.weighted_rms
}

///|
pub fn FitResult::rank(self : FitResult) -> Int32 {
  self.rank
}

///|
/// Smallest/largest absolute QR diagonal, a diagnostic, not a matrix condition number.
pub fn FitResult::diagonal_ratio(self : FitResult) -> Double {
  self.diagonal_ratio
}

///|
/// Minimize sum_i weights[i]*(f(direction[i])-values[i])^2.
/// Dense solver is deliberately limited to degree<=12 and at most 4096 rows.
/// Zero weights are permitted; negative/nonfinite weights and rank deficiency are errors.
/// rank_tolerance is relative to the first pivot and must lie in [1.0e-15,1.0e-2].
pub fn fit_weighted(
  degree : Int32,
  directions : Array[Direction],
  values : Array[Double],
  weights : Array[Double],
  rank_tolerance~ : Double = 1.0e-12,
) -> Result[FitResult, HarmonicError] {
  if degree < 0 || degree > 12 {
    return Err(InvalidInput("dense fitting supports degrees 0..12"))
  }
  let cols = (degree + 1) * (degree + 1)
  let rows = directions.length()
  if rows < cols || rows > 4096 {
    return Err(InvalidInput("fit needs (degree+1)^2..4096 rows"))
  }
  if values.length() != rows {
    return Err(ShapeMismatch(rows, values.length()))
  }
  if weights.length() != rows {
    return Err(ShapeMismatch(rows, weights.length()))
  }
  if !finite(rank_tolerance) ||
    rank_tolerance < 1.0e-15 ||
    rank_tolerance > 0.01 {
    return Err(InvalidInput("invalid relative QR rank tolerance"))
  }
  let mut max_weight = 0.0
  let mut value_scale = 0.0
  for i = 0; i < rows; i = i + 1 {
    if !finite(values[i]) || !finite(weights[i]) || weights[i] < 0.0 {
      return Err(InvalidInput("fit data or weights invalid"))
    }
    max_weight = @cmp.maximum(max_weight, weights[i])
    value_scale = @cmp.maximum(value_scale, values[i].abs())
  }
  if max_weight == 0.0 {
    return Err(SingularSystem(0))
  }
  if value_scale == 0.0 {
    value_scale = 1.0
  }
  let matrix : Array[Array[Double]] = []
  let rhs = Array::make(rows, 0.0)
  for i = 0; i < rows; i = i + 1 {
    let row = real_basis(degree, directions[i]).unwrap()
    let w = (weights[i] / max_weight).sqrt()
    for j = 0; j < cols; j = j + 1 {
      row[j] = row[j] * w
    }
    matrix.push(row)
    rhs[i] = values[i] / value_scale * w
  }
  let permutation = Array::make(cols, 0)
  for j = 0; j < cols; j = j + 1 {
    permutation[j] = j
  }
  let mut first_pivot = 0.0
  let mut min_pivot = 0.0
  for k = 0; k < cols; k = k + 1 {
    // Recompute trailing column norms instead of subtracting nearly equal squared norms.
    let mut pivot = k
    let mut best = -1.0
    for j = k; j < cols; j = j + 1 {
      let mut squared = 0.0
      for i = k; i < rows; i = i + 1 {
        squared = squared + matrix[i][j] * matrix[i][j]
      }
      if squared > best {
        best = squared
        pivot = j
      }
    }
    let norm = best.sqrt()
    if k == 0 {
      first_pivot = norm
      min_pivot = norm
    }
    if norm <= rank_tolerance * first_pivot || norm == 0.0 {
      return Err(SingularSystem(k))
    }
    min_pivot = @cmp.minimum(min_pivot, norm)
    if pivot != k {
      for i = 0; i < rows; i = i + 1 {
        let temp = matrix[i][k]
        matrix[i][k] = matrix[i][pivot]
        matrix[i][pivot] = temp
      }
      let temp = permutation[k]
      permutation[k] = permutation[pivot]
      permutation[pivot] = temp
    }
    let alpha = if matrix[k][k] >= 0.0 { -norm } else { norm }
    let vector = Array::make(rows - k, 0.0)
    for i = k; i < rows; i = i + 1 {
      vector[i - k] = matrix[i][k]
    }
    vector[0] = vector[0] - alpha
    let mut squared = 0.0
    for v in vector {
      squared = squared + v * v
    }
    let beta = 2.0 / squared
    for j = k + 1; j < cols; j = j + 1 {
      let mut dot = 0.0
      for i = k; i < rows; i = i + 1 {
        dot = dot + vector[i - k] * matrix[i][j]
      }
      for i = k; i < rows; i = i + 1 {
        matrix[i][j] = matrix[i][j] - beta * dot * vector[i - k]
      }
    }
    let mut dot = 0.0
    for i = k; i < rows; i = i + 1 {
      dot = dot + vector[i - k] * rhs[i]
    }
    for i = k; i < rows; i = i + 1 {
      rhs[i] = rhs[i] - beta * dot * vector[i - k]
    }
    matrix[k][k] = alpha
    for i = k + 1; i < rows; i = i + 1 {
      matrix[i][k] = 0.0
    }
  }
  let solution = Array::make(cols, 0.0)
  for k = cols - 1; k >= 0; k = k - 1 {
    let mut sum = rhs[k]
    for j = k + 1; j < cols; j = j + 1 {
      sum = sum - matrix[k][j] * solution[j]
    }
    solution[k] = sum / matrix[k][k]
  }
  let coefficients = Array::make(cols, 0.0)
  for j = 0; j < cols; j = j + 1 {
    coefficients[permutation[j]] = solution[j] * value_scale
  }
  let spectrum = match Spectrum::new(degree, coefficients) {
    Ok(s) => s
    Err(e) => return Err(e)
  }
  // Residuals are scaled by input magnitude and maximum weight to limit overflow.
  let mut squared_error = 0.0
  let mut total_weight = 0.0
  for i = 0; i < rows; i = i + 1 {
    let prediction = match spectrum.evaluate(directions[i]) {
      Ok(v) => v
      Err(e) => return Err(e)
    }
    let error = prediction / value_scale - values[i] / value_scale
    let w = weights[i] / max_weight
    squared_error = squared_error + w * error * error
    total_weight = total_weight + w
  }
  let weighted_rms = (squared_error / total_weight).sqrt() * value_scale
  if !finite(weighted_rms) {
    return Err(InvalidInput("fit residual overflow"))
  }
  Ok({
    spectrum,
    weighted_rms,
    rank: cols,
    diagonal_ratio: min_pivot / first_pivot,
  })
}
