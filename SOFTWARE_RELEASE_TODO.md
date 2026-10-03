# Complete before public release

- Add all author names/affiliations and a software license.
- Record R, `shiny`, `magick`, and `jsonlite` versions used for the published app.
- Record Python, NumPy, pandas, Pillow, and scikit-learn versions used to rebuild the model.
- Add a release tag and a persistent archive DOI; verify that the tagged repository includes the exact app and model used in the paper.
- State whether the 710 source images can be redistributed; do not include them without permission.
- Consider adding the complete feature matrix and fold-wise validation predictions if permitted.
- Smoke-test `shiny::runApp(".")` from a clean R installation after extraction.
