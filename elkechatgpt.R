library(httr)
library(magick)

# Function to download and format PNG with detailed debugging
download_and_format_png <- function(cid, output_dir, width = NULL, height = NULL, background_color = "white") {
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/PNG")

  response <- GET(url)

  if (status_code(response) == 200) {
    file_path <- file.path(output_dir, paste0(cid, ".png"))
    writeBin(content(response, "raw"), file_path)

    img <- image_read(file_path)

    # Resize the image if width and height are provided
    if (!is.null(width) && !is.null(height)) {
      img <- image_resize(img, paste0(width, "x", height))
    }

    # Set background color to white
    img <- image_background(img, color = background_color, flatten = TRUE)
    print(c("img_back",img))
    plot(img)
    img <- image_transparent(img, color = background_color, fuzz = 15)
    print(c("img_opacity",img))
    plot(img)
    img <- image_border(img, color = "black", geometry = "3x3", operator = "copy")
    print(c("img_border",img))
    plot(img)

    # Save the final formatted image
    image_write(img, path = file_path, format = "png")

    print(paste("Downloaded and formatted PNG for CID", cid))
  } else {
    print(paste("Failed to download PNG for CID", cid, ", status code:", status_code(response)))
  }
}

# Example usage
cids <- c(443939, 65348, 62770, 36462, 23663976)

cids <- readRDS("data/elkepcids.RDS")
output_dir <- "www/elkeimg"
width <- 200
height <- 200
background_color <- "white"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (cid in cids) {
  download_and_format_png(cid, output_dir, width, height, background_color)
}
