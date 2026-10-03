# Start here — flat GumOptic app

This ZIP places `app.R`, `model/`, and `www/` side by side to avoid the missing model/logo path error. Extract the complete ZIP. In RStudio, open this folder, open `app.R`, and click **Run App**. Do not copy or run app.R alone.

Required R packages (install once):
```r
install.packages(c("shiny", "magick", "jsonlite"))
```

Expected files beside app.R:
- `model/metadata.json`
- `model/forest_nodes.csv`
- `model/centroids.csv`
- `model/scaler.csv`
- `www/gumoptic_logo.svg`

The trained model is bundled. Readers only upload an unknown image. If a red error remains, send its exact text. This screening model is not a guaranteed species confirmation.
