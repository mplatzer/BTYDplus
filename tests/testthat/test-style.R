if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    skip("skip lint checks for the time being")
    lintr::expect_lint_free()
  })
}
