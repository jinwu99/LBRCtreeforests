rm(list=ls())
library(pdftools)
library(magick)

input_pdf <- "combined_trees_3x1.pdf"  # Your input PDF file
output_png <- "combined_trees_3x1.png"  # Your output PNG file

input_pdf <- "LC Brier score.pdf"  # Your input PDF file
output_png <- "LC Brier score.png"  # Your output PNG file
# Convert PDF to PNG (at 300 dpi for high quality)
pdf_images <- pdf_convert(input_pdf, format = "png", dpi = 600)
# Rename and save the first page as a PNG file
file.rename(pdf_images[1], output_png)

############################

setwd("C:/Users/jinwu/Desktop/LBRCtree_forests/analysis/results_real_data")

files <- c(
  "LTRC-CIT.pdf", "LBRC-CIT-C.pdf", "LBRC-CIT-F.pdf"
)
imgs <- lapply(files, image_read_pdf, density = 1000)
imgs <- do.call(c, imgs)
# 3) Build each row by appending three images horizontally
rows <- lapply(seq(1, 3, by = 3), function(i) {
  image_append(imgs[i:(i+2)], stack = F)
})
# 4) Stack the three rows vertically
grid3 <- image_append(do.call(c, rows), stack = F)
image_write(grid3, path = "combined_trees_1x3.pdf", format = "pdf")

####################

setwd("C:/Users/jinwu/Desktop/LBRCtree_forests/analysis/results_real_data")

library(magick)

# 1. Read in your three logos
logo3 <- image_read_svg("Ewha_logo.svg")
logo2 <- image_read_svg("Yonsei_logo.svg")
logo1 <- image_read("NCC_logo.png")

# 2. Normalize them to the same square size S (use the smallest width to avoid up‐sampling)
infos  <- lapply(list(logo1, logo2, logo3), image_info)
widths <- sapply(infos, function(i) i$width)
S      <- min(widths)

logo1_r <- image_scale(logo1, paste0(S, "x", S))
logo2_r <- image_scale(logo2, paste0(S, "x", S))
logo3_r <- image_scale(logo3, paste0(S, "x", S))

# 3. Compute equilateral triangle geometry
d        <- S                    # choose center‐to‐center distance = S
H_tri    <- sqrt(3)/2 * d        # vertical offset between top and bottom centers
W_canvas <- 2 * S                # canvas width
H_canvas <- ceiling(S + H_tri)   # canvas height

# 4. Make a blank canvas
canvas <- image_blank(width = W_canvas, height = H_canvas, color = "white")

# 5. Compute placement offsets (top‐left corner of each logo)
#    Offsets are in "+x+y" form; origin is top‐left of the canvas.
top_off       <- paste0("+", round(S/2), "+", 0)
bottom_left   <- paste0("+", 0,           "+", round(H_tri))
bottom_right  <- paste0("+", S,           "+", round(H_tri))

# 6. Composite each logo onto the canvas
canvas <- image_composite(canvas, logo1_r, offset = top_off)
canvas <- image_composite(canvas, logo2_r, offset = bottom_left)
canvas <- image_composite(canvas, logo3_r, offset = bottom_right)

# 7. Save the result
image_write(canvas, path = "triangle_logos.png", format = "png")


#################################################

library(magick)

# 1. Read in each PDF (first page only) at a good resolution
img1 <- image_read_pdf("combined_trees_3x1.pdf", density = 300)[1]
img2 <- image_read_pdf("LC Brier score.pdf", density = 300)[1]

# 2. Get their pixel dimensions
info1 <- image_info(img1)
info2 <- image_info(img2)

# 3. Compute the new height for img2: 2/3 of img1’s height
new_height <- round(info1$height * 1/2)

# 4. Resize img2 to that height (keeping aspect ratio)
img2_resized <- image_scale(img2, geometry = paste0("x", new_height))

# 5. Pad img2 on top/bottom so its canvas height equals img1’s height,
#    centering the content vertically, with a white background
info2r <- image_info(img2_resized)
img2_padded <- image_extent(
  img2_resized,
  geometry = paste0(info2r$width, "x", info1$height),
  gravity  = "center",
  color    = "white"
)

# 6. Append them horizontally
combined <- image_append(c(img1, img2_padded), stack = FALSE)

# 7. Write out to a single-page PDF
image_write(combined, path = "three_trees_Brier.pdf", format = "pdf")

###################

library(magick)

# 1. Read in your three inputs:
img1 <- image_read_pdf("combined_trees_3x1.pdf", density = 300)[1]
img2 <- image_read_pdf("LC Brier score.pdf",   density = 300)[1]
img3 <- image_read("Poster_table.png",         density = 300)  # PNG

# 2. Figure out img1’s height:
info1 <- image_info(img1)
H <- info1$height

# 3. Decide how to split that height between img2 (top) and img3 (bottom).
#    Here we give each half of the total height; you can tweak f_top.
f_top    <- 0.7
H2       <- round(H * f_top)    # height for img2
H3       <- H - H2              # remainder for img3

# 4. Resize each to its allotted height (preserving aspect):
img2_r <- image_scale(img2, paste0("x", H2))
img3_r <- image_scale(img3, paste0("x", H3))

# 5. Find the maximum width of those two so the column is uniform:
w2 <- image_info(img2_r)$width
w3 <- image_info(img3_r)$width
W_right <- max(w2, w3)

# 6. Pad each to exactly (W_right × its height), centering them in each block:
img2_block <- image_extent(
  img2_r,
  geometry = paste0(W_right, "x", H2),
  gravity  = "center",
  color    = "white"
)
img3_block <- image_extent(
  img3_r,
  geometry = paste0(W_right, "x", H3),
  gravity  = "center",
  color    = "white"
)

# 7. Stack those two blocks vertically to form the right column:
right_column <- image_append(c(img2_block, img3_block), stack = TRUE)

# 8. Finally, append that column to the right of img1:
final    <- image_append(c(img1, right_column), stack = FALSE)

# 9. Write out a single‐page PDF:
image_write(final, path = "three_trees_Brier_table.pdf", format = "pdf")
