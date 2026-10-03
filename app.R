# GumOptic Consensus: 13-species unknown gum photo screener
# Reader-facing app. The portable, already-trained forest and centroid tables are in model/.
# Required R packages: shiny, magick, jsonlite
library(shiny)
library(magick)
library(jsonlite)

load_gumoptic_model <- function(model_dir = NULL) {
  if (is.null(model_dir)) {
    candidates <- c("model", "app/model", "../app/model", "GumOptic_Final_13_species/app/model")
    found <- candidates[file.exists(file.path(candidates, "metadata.json"))]
    if (!length(found)) stop(paste0("Cannot find bundled app/model/metadata.json. Extract the complete ZIP (including app/model/) and launch the app from its app folder. Current working directory: ", getwd()))
    model_dir <- found[1]
  }
  meta <- jsonlite::fromJSON(file.path(model_dir, "metadata.json"), simplifyVector = TRUE)
  nodes <- read.csv(file.path(model_dir, "forest_nodes.csv"), check.names = FALSE)
  scaler <- read.csv(file.path(model_dir, "scaler.csv"), stringsAsFactors = FALSE)
  cent <- read.csv(file.path(model_dir, "centroids.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  features <- as.character(unlist(meta$features))
  classes <- as.character(unlist(meta$classes))
  if (!identical(as.character(meta$feature_version), "gumoptic23-v1")) stop("Incompatible model feature version.")
  if (!all(features %in% scaler$Feature) || !all(features %in% names(cent))) stop("Model feature files are incomplete.")
  scaler <- scaler[match(features, scaler$Feature), ]
  cent <- cent[match(classes, cent$Species), , drop = FALSE]
  nodes <- nodes[order(nodes$tree, nodes$node), , drop = FALSE]
  trees <- lapply(seq_len(as.integer(meta$n_trees)), function(i) { d <- nodes[nodes$tree == i, , drop = FALSE]; d[order(d$node), , drop = FALSE] })
  list(meta = meta, features = features, classes = classes, trees = trees,
       centroids = as.matrix(cent[, features, drop = FALSE]),
       med = scaler$Median, iqr = scaler$IQR,
       temp = as.numeric(meta$centroid_temp),
       reject = as.numeric(meta$oof_distance_q95),
       prob_cols = paste0("p", seq_along(classes)))
}

# Must match build_model.py's 23 optical descriptors and image resize convention.
extract_gum_features <- function(path) {
  im <- magick::image_read(path)[1]
  im <- tryCatch(magick::image_auto_orient(im), error = function(e) im)
  im <- magick::image_resize(im, "96x54!")
  im <- magick::image_convert(im, colorspace = "sRGB")
  a <- magick::image_data(im, channels = "rgb")
  R <- as.numeric(a[1, , ]) / 255
  G <- as.numeric(a[2, , ]) / 255
  B <- as.numeric(a[3, , ]) / 255
  mx <- pmax(R, G, B); mn <- pmin(R, G, B); delta <- mx - mn
  S <- ifelse(mx > 1e-8, delta / (mx + 1e-8), 0); V <- mx
  H <- rep(0, length(mx)); nz <- delta > 1e-8
  ir <- nz & mx == R; ig <- nz & mx == G; ib <- nz & mx == B
  H[ir] <- ((G[ir] - B[ir]) / (delta[ir] + 1e-8)) %% 6
  H[ig] <- (B[ig] - R[ig]) / (delta[ig] + 1e-8) + 2
  H[ib] <- (R[ib] - G[ib]) / (delta[ib] + 1e-8) + 4
  H <- H / 6
  roi <- S > .12 & V > .05 & V < .99
  warm <- roi & R >= .82 * G & G >= .72 * B
  if (sum(warm) >= 20) roi <- warm
  if (sum(roi) < 15) roi <- rep(TRUE, length(R))
  stat4 <- function(x, nm) setNames(c(mean(x[roi]), sd(x[roi]),
      unname(quantile(x[roi], .10)), unname(quantile(x[roi], .90))),
      paste0(nm, c("_mean", "_sd", "_q10", "_q90")))
  gray <- .299 * R + .587 * G + .114 * B
  counts <- hist(gray[roi], breaks = seq(0, 1, length.out = 17),
                 include.lowest = TRUE, plot = FALSE)$counts
  p <- counts[counts > 0] / sum(counts)
  entropy <- -sum(p * log2(p))
  grad <- abs(diff(gray)); rg <- R - G; yb <- (R + G) / 2 - B
  c(stat4(R, "R"), stat4(G, "G"), stat4(B, "B"),
    hue_mean = mean(H[roi]), hue_sd = sd(H[roi]),
    sat_mean = mean(S[roi]), brightness_mean = mean(V[roi]),
    gray_mean = mean(gray[roi]), gray_sd = sd(gray[roi]),
    gray_entropy = entropy, edge_mean = mean(grad), edge_frac = mean(grad > .08),
    colorfulness = sqrt(sd(rg)^2 + sd(yb)^2) + .3 * sqrt(mean(rg)^2 + mean(yb)^2),
    roi_fraction = mean(roi))
}

predict_gumoptic <- function(f, m) {
  x <- as.numeric(f[m$features]); names(x) <- m$features
  if (length(x) != length(m$features) || any(!is.finite(x))) stop("Feature extraction returned missing/non-finite values; try a clear RGB JPG/PNG image.")
  k <- length(m$classes); rf <- numeric(k)
  for (t in seq_along(m$trees)) {
    tr <- m$trees[[t]]; node_id <- 0L
    for (step in seq_len(1000L)) {
      row <- tr[node_id + 1L, , drop = FALSE]
      fi <- as.integer(row$feature[1])
      if (fi < 0L) {
        rf <- rf + as.numeric(row[1, m$prob_cols, drop = TRUE])
        break
      }
      node_id <- if (x[fi + 1L] < row$threshold[1]) as.integer(row$left[1]) else as.integer(row$right[1])
    }
  }
  rf <- rf / length(m$trees); rf <- rf / sum(rf)
  z <- (x - m$med) / m$iqr
  d <- sqrt(rowSums((m$centroids - matrix(z, nrow = k, ncol = length(z), byrow = TRUE))^2))
  sim <- exp(-d / max(m$temp, 1e-8)); sim <- sim / sum(sim)
  fused <- .70 * rf + .30 * sim
  data.frame(Species = m$classes, RF_vote = rf, Centroid_similarity = sim,
             Fusion_score = fused, Centroid_distance = d, check.names = FALSE)
}

model <- tryCatch(load_gumoptic_model(), error = function(e) e)
ui <- fluidPage(
  tags$head(tags$link(rel = "icon", type = "image/svg+xml", href = "gumoptic_logo.svg"),
    tags$style(HTML("body{background:#f5faf9;color:#25333b}.brand{display:flex;align-items:center;gap:16px;padding:12px 18px;margin:10px 0 18px;border-radius:18px;background:linear-gradient(110deg,#fff,#e6f7f0);box-shadow:0 3px 14px #1232}.brand img{width:88px;height:88px}.brand h1{font-weight:850;font-size:34px;margin:0;background:linear-gradient(90deg,var(--c1),var(--c2),var(--c3));-webkit-background-clip:text;background-clip:text;color:transparent;animation:titleHue 8s ease-in-out infinite alternate}.brand p{margin:5px 0;color:#526a66}.score{margin:9px 0}.scoretop{display:flex;justify-content:space-between;font-weight:700}.track{height:18px;background:#e5eeea;border-radius:10px;overflow:hidden}.fill{height:100%;border-radius:10px;transition:width .65s ease}.notice{padding:12px;border-radius:12px;background:#fff7e7;border-left:5px solid #e7a522;margin:12px 0}.metric{padding:12px;border-radius:12px;background:#fff;border-left:5px solid #179b78;margin:10px 0}@keyframes titleHue{to{filter:hue-rotate(45deg)}}@media(prefers-reduced-motion:reduce){.brand h1{animation:none}}"))),
  uiOutput("brand"),
  sidebarLayout(
    sidebarPanel(width = 4,
      h4("Identify an unknown gum image"),
      p("The bundled 13-species reference model loads automatically. Readers upload one image; they do not need the training archives."),
      fileInput("unknown", "Unknown sample (JPG, PNG, or TIFF)", accept = c(".jpg", ".jpeg", ".png", ".tif", ".tiff")),
      downloadButton("download_features", "Download optical features CSV"),
      uiOutput("model_status"),
      actionButton("cycle", "Cycle title/bar colours", class = "btn-success"),
      p(class = "text-muted", "Screening aid only. The model is closed-set and does not prove botanical identity.")),
    mainPanel(width = 8,
      tabsetPanel(
        tabPanel("Identify", h4("Candidate ranking"), uiOutput("decision"),
          uiOutput("score_bars"), tableOutput("score_table"),
          h4("Optical feature report"), tableOutput("feature_table")),
        tabPanel("About / methods", h4("Bundled reference model"),
          p("This package contains a portable 400-tree random-forest model trained on 710 images across 13 supplied species archives, plus a robust-scaled centroid reference. The app combines the RF vote and centroid similarity as 70:30. It needs no training ZIPs at prediction time."),
          p("The model development used 23 color, brightness, grayscale-texture, edge, colorfulness, and ROI descriptors from resized images. The 5-fold image-level cross-validation was not specimen-independent because tree/specimen identifiers were not available. The reported scores are not external-validation estimates and can be optimistic if photos are repeated views of the same biological material."),
          p("RGB means are displayed on a 0-255 scale. Lighting, camera, background, gum thickness, and ROI segmentation can affect descriptors. Use consistent imaging and independently verify conclusions."),
          p("The fusion score is a ranking score, not a calibrated probability. The alert is a distance-based review heuristic, not a reliable detector for every species outside the training set."),
          p("For local use: install.packages(c('shiny','magick','jsonlite')); then run shiny::runApp('app'). GitHub stores the code/model files but does not itself host a running Shiny application."))
      ))
  )
)
server <- function(input, output, session) {
  palette_id <- reactiveVal(0L)
  palettes <- list(c("#168650", "#1d63b5", "#f1a426"), c("#7b2cbf", "#e63946", "#f4a261"), c("#005f73", "#0a9396", "#ee9b00"), c("#b5179e", "#4361ee", "#4cc9f0"))
  colors <- reactive(palettes[[palette_id() + 1L]])
  observeEvent(input$cycle, { palette_id((palette_id() + 1L) %% length(palettes)) })
  output$brand <- renderUI({ cc <- colors(); div(class = "brand", style = paste0("--c1:", cc[1], ";--c2:", cc[2], ";--c3:", cc[3]),
      tags$img(src = "gumoptic_logo.svg", alt = "GumOptic Consensus logo"),
      div(h1("GumOptic Consensus"), p("13-species gum image identification screening"))) })
  output$model_status <- renderUI({
    if (inherits(model, "error")) div(class = "notice", strong("Model loading problem"), p(conditionMessage(model)))
    else div(class = "metric", strong("Bundled model ready"), p(paste(length(model$classes), "species;", model$meta$n_images, "training photos;", model$meta$n_trees, "trees.")))
  })
  feats <- reactive({ req(input$unknown); extract_gum_features(input$unknown$datapath) })
  # Capture prediction errors and display them in the page instead of leaving the name area blank.
  preds <- reactive({
    req(input$unknown)
    if (inherits(model, "error")) return(list(error = paste("Bundled model did not load:", conditionMessage(model))))
    tryCatch(list(scores = predict_gumoptic(feats(), model)), error = function(e) list(error = conditionMessage(e)))
  })
  output$decision <- renderUI({
    req(input$unknown); result <- preds()
    if (!is.null(result$error)) return(div(class = "alert alert-danger", strong("Could not classify this image"), p(result$error), p("Confirm the full ZIP was extracted and the model files remain in app/model/.")))
    s <- result$scores[order(-result$scores$Fusion_score), , drop = FALSE]
    gap <- s$Fusion_score[1] - s$Fusion_score[2]
    alert <- s$Centroid_distance[1] > model$reject || gap < .03
    div(class = if (alert) "notice" else "metric",
      h3(if (alert) paste("Review candidate:", s$Species[1]) else paste("Leading species candidate:", s$Species[1])),
      p(sprintf("Fusion ranking score: %.3f", s$Fusion_score[1])),
      p(if (alert) "Low support / close candidates: verify independently. The name above is still the top-ranked class." else "This is a closed-set screening result, not a confirmed identification."))
  })
  output$score_bars <- renderUI({
    req(input$unknown); result <- preds()
    if (!is.null(result$error)) return(NULL)
    s <- result$scores[order(-result$scores$Fusion_score), , drop = FALSE]
    lapply(seq_len(nrow(s)), function(i) { hue <- (i * 31 + palette_id() * 45) %% 360; col <- sprintf("hsl(%d,68%%,46%%)", hue)
      div(class = "score", div(class = "scoretop", span(s$Species[i]), span(sprintf("%.1f score", 100*s$Fusion_score[i]))),
        div(class = "track", div(class = "fill", style = sprintf("width:%.2f%%;background:%s", 100*s$Fusion_score[i], col)))) })
  })
  output$score_table <- renderTable({ req(input$unknown); result <- preds(); if (!is.null(result$error)) return(NULL); s <- result$scores; s[order(-s$Fusion_score), ] }, digits = 4)
  output$feature_table <- renderTable({ req(input$unknown); f <- feats(); vals <- as.numeric(f); names(vals) <- names(f); rgb <- names(vals) %in% c("R_mean", "G_mean", "B_mean"); vals[rgb] <- vals[rgb] * 255
    data.frame(Feature = names(vals), Value = round(vals, 4), Unit = ifelse(rgb, "8-bit mean (0-255)", ifelse(names(vals) == "gray_entropy", "bits", "descriptor units"))) }, digits = 4)
  output$download_features <- downloadHandler(filename = function() "gumoptic_unknown_optical_features.csv", content = function(file) {
    f <- as.data.frame(as.list(feats()), check.names = FALSE); f$R_mean_8bit <- f$R_mean * 255; f$G_mean_8bit <- f$G_mean * 255; f$B_mean_8bit <- f$B_mean * 255; write.csv(f, file, row.names = FALSE)
  })
}
shinyApp(ui, server)
