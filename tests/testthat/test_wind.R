### Test file for wind functions

# Load libraries
library(dplyr)

### Test for function convert_to_cardinal --------------------------------------
test_that("convert_to_cardinal works correctly", {

  # Basic cardinal points
  expect_equal(convert_to_cardinal(0), "N")
  expect_equal(convert_to_cardinal(45), "NE")
  expect_equal(convert_to_cardinal(90), "E")
  expect_equal(convert_to_cardinal(135), "SE")
  expect_equal(convert_to_cardinal(180), "S")
  expect_equal(convert_to_cardinal(225), "SW")
  expect_equal(convert_to_cardinal(270), "W")
  expect_equal(convert_to_cardinal(315), "NW")
  expect_equal(convert_to_cardinal(360), "N")  # Edge case

  # Boundary values
  expect_equal(convert_to_cardinal(22.5), "NE")
  expect_equal(convert_to_cardinal(67.5), "E")
  expect_equal(convert_to_cardinal(112.5), "SE")
  expect_equal(convert_to_cardinal(157.5), "S")
  expect_equal(convert_to_cardinal(202.5), "SW")
  expect_equal(convert_to_cardinal(247.5), "W")
  expect_equal(convert_to_cardinal(292.5), "NW")
  expect_equal(convert_to_cardinal(337.5), "N")

  # Mid-range values
  expect_equal(convert_to_cardinal(30), "NE")
  expect_equal(convert_to_cardinal(75), "E")
  expect_equal(convert_to_cardinal(120), "SE")
  expect_equal(convert_to_cardinal(165), "S")
  expect_equal(convert_to_cardinal(210), "SW")
  expect_equal(convert_to_cardinal(255), "W")
  expect_equal(convert_to_cardinal(300), "NW")
  expect_equal(convert_to_cardinal(345), "N")

  # Edge cases
  expect_equal(convert_to_cardinal(360), "N")

  # NA handling
  expect_true(is.na(convert_to_cardinal(NA)))
})


### Test for function cardinal_to_degree ---------------------------------------
test_that("cardinal_to_degree works correctly", {

  # Basic cardinal points
  expect_equal(cardinal_to_degree("N"), 0)
  expect_equal(cardinal_to_degree("NE"), 45)
  expect_equal(cardinal_to_degree("E"), 90)
  expect_equal(cardinal_to_degree("SE"), 135)
  expect_equal(cardinal_to_degree("S"), 180)
  expect_equal(cardinal_to_degree("SW"), 225)
  expect_equal(cardinal_to_degree("W"), 270)
  expect_equal(cardinal_to_degree("NW"), 315)

  # Handling NA input
  expect_true(is.na(cardinal_to_degree(NA)))

  # Invalid input should return NA
  expect_true(is.na(cardinal_to_degree("INVALID")))
  expect_true(is.na(cardinal_to_degree("")))
  expect_true(is.na(cardinal_to_degree("north")))  # Case-sensitive check
})

### Test for function get_mode_dir ---------------------------------------------
test_that("get_mode_dir works correctly", {

  # Basic case: clear mode
  expect_equal(get_mode_dir(c("N", "N", "E", "N", "S")), "N")

  # Tie-breaking: returns first most frequent value
  expect_equal(get_mode_dir(c("N", "E", "E", "N")), "N")
  # Both "N" and "E" appear twice, "N" comes first

  # All unique values: should return the first element
  expect_equal(get_mode_dir(c("N", "E", "S", "W")), "N")

  # Single element input: should return the same value
  expect_equal(get_mode_dir(c("N")), "N")

  # Handling NA values: should ignore NA
  expect_equal(get_mode_dir(c("N", "N", "E", "N", "S", NA)), "N")
  # "N" is still the most frequent

  # Only NA values: should return NA
  expect_true(is.na(get_mode_dir(c(NA, NA, NA))))
})


### Test for function compute_wind_episode -------------------------------------
test_that("compute_wind_episode works correctly", {
  # Sample wind data
  wind_df <- data.frame(
    datetime = as.POSIXct(c("2023-01-01 12:00:00", "2023-01-01 13:00:00",
                            "2023-01-01 14:00:00"), tz = "UTC"),
    FF = c(10, 12, 14),
    DD = c(90, 180, 270),
    cardDir = c("E", "S", "W")
  )

  # Sample episode data
  episode <- data.frame(
    location1 = c(1.2, 2.3, 0.5),
    row.names = c("2023-01-01 12:00:00", "2023-01-01 13:00:00",
                  "2023-01-01 14:00:00")
  )

  mode_card <- get_mode_dir(wind_df$cardDir)
  mode_deg <- cardinal_to_degree(mode_card)
  mode_deg_mean <- mean(cardinal_to_degree(wind_df$cardDir[wind_df$cardDir ==
                                                             mode_card]))

  # Test with valid input
  result <- compute_wind_episode(episode, "location1", 1, wind_df, 2)
  expect_equal(result$FF, 14)  # The wind speed at t0 = delta + 1
  # Mean direction for most frequent cardDir
  expect_equal(result$DD, mode_deg_mean)
  expect_equal(result$DD_deg, mode_deg)  # W -> 270 degrees
  expect_equal(result$cardDir, mode_card) # Most frequent direction
  expect_equal(result$DD_t0, 270)  # Direction at t0
  expect_equal(result$cardDir_t0, "W")

  # Test when no timestamps match
  empty_episode <- data.frame(location1 = c(3.0),
                            row.names = "2025-01-01 00:00:00")  # Future date
  result_empty <- compute_wind_episode(empty_episode, "location1", 1,
                                        wind_df, 2)
  expect_true(all(is.na(result_empty)))

  # Test with a single timestamp
  single_episode <- data.frame(location1 = c(1.5),
                               row.names = "2023-01-01 12:00:00")
  result_single <- compute_wind_episode(single_episode, "location1", 1,
                                        wind_df, 0)
  expect_equal(result_single$FF, 10)  # Wind speed at t0
  expect_equal(result_single$DD, 90)  # Mean wind direction
  expect_equal(result_single$DD_deg, 90)  # S -> 180 degrees
  expect_equal(result_single$cardDir, "E") # Only direction
  expect_equal(result_single$DD_t0, 90)  # Direction at t0

  # Test with NA values with a timestamp not in the wind data
  episode_with_na <- data.frame(location1 = c(NA),
                                row.names = "2023-01-01 15:00:00")
  result_na <- compute_wind_episode(episode_with_na, "location1", 1, wind_df, 0)
  expect_true(all(is.na(result_na)))  # Should return all NAs

})
