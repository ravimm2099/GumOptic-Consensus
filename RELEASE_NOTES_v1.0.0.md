# GumOptic-Consensus v1.0.0

Publication release accompanying the GumOptic-Consensus manuscript.

## Included
- R Shiny reader-facing application.
- Portable 400-tree Random Forest.
- Robust scaler and 13 species centroids.
- Model metadata and feature schema `gumoptic23-v1`.
- Python model-building/export script.
- Five-fold image-level fused confusion matrix and species recall summary.
- User/troubleshooting supplementary document.

## Model summary
- Reference photographs: 710.
- Species: 13.
- Optical descriptors: 23.
- Final Random Forest: 400 trees.
- Consensus weights: RF 0.70; centroid similarity 0.30.
- RF image-level CV: accuracy 75.63%, balanced accuracy 74.25%, macro-F1 74.46%.
- Consensus image-level CV: accuracy 75.21%, balanced accuracy 73.85%, macro-F1 73.96%.

## Scope
Validation is stratified five-fold at the image level. Specimen/tree identifiers and an independent external test set were unavailable. Performance therefore describes internal image-level discrimination and should not be interpreted as independent field-identification accuracy. The app is a closed-set screening aid.
