if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    skip()
    lintr::expect_lint_free()
  })
}
