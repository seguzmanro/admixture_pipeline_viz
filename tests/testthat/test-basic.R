library(testthat)
library(admixspatial)

# Test data processing functions
test_that("read_input_data validates input files", {
  # This would require actual test files
  # For now, just test that the function exists and has correct structure
  expect_true(is.function(read_input_data))
})

test_that("genetic_palette generates correct number of colors", {
  colors <- genetic_palette(5)
  expect_equal(length(colors), 5)
  expect_true(all(grepl("^#", colors)))  # All should be hex colors
})

test_that("theme_genomics returns ggplot theme", {
  theme <- theme_genomics()
  expect_true(inherits(theme, "theme"))
})

# Test main functions exist
test_that("main visualization functions exist", {
  expect_true(is.function(visualize_admixture))
  expect_true(is.function(quick_visualize))
  expect_true(is.function(execute_pipeline))
})