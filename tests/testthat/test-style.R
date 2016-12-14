if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    skip_on_cran()
    lintr::expect_lint_free()
  })
}
