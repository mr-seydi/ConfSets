# Starter script for the ConfSets package hex sticker.
# Run from the package root after installing hexSticker.
# Adjust the artwork placement, colors, typography, and sticker layout below.

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

# Keep the SVG as the editable source artwork in inst/sticker/.
# The PNG copy is used here so hexSticker can place it reliably as a subplot.
curve_image <- "inst/sticker/aurora_curve.png"

hexSticker::sticker(
  subplot = curve_image,
  package = "ConfSets",
  p_size = 42,
  p_color = "#afbcbf",
  p_x = 1,
  p_y = 0.65,
  s_x = 1,
  s_y = 1.2,
  s_width = 0.35,
  # aurora_curve.png is 1600 x 960 pixels (width / height).
  asp = 1600 / 960,
  h_fill = "#253b52",
  h_color = "#4c656f",
  h_size = 1.5,
  filename = "man/figures/ConfSets-hex.png",
  dpi = 600
)

# pkgdown discovers the package logo from this conventional filename.
file.copy("man/figures/ConfSets-hex.png", "man/figures/logo.png", overwrite = TRUE)
