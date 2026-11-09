///|
test "signed order packing is bijective for every supported mode" {
  for l = 0; l <= MAX_DEGREE; l = l + 1 {
    for m = -l; m <= l; m = m + 1 {
      let i = coefficient_index(l, m).unwrap()
      assert_eq(degree_order(i).unwrap(), (l, m))
    }
  }
  assert_eq(coefficient_count(64).unwrap(), 4225)
  assert_eq(coefficient_index(-1, 0), Err(InvalidDegree(-1)))
  assert_eq(coefficient_index(1, 2), Err(InvalidOrder(1, 2)))
  assert_true(
    match degree_order(4225) {
      Err(_) => true
      Ok(_) => false
    },
  )
}
