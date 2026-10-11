## Test environments

* Local Windows 11 x64, R 4.6.1
* Win-builder Windows Server 2022, R-release 4.6.1
* Win-builder Windows Server 2022, R-devel
  (2026-10-09 r90655 ucrt)
* GitHub Actions: package-check workflow completed successfully

## R CMD check results

0 errors | 0 warnings

The local Windows check and both win-builder checks each reported
one CRAN incoming feasibility NOTE.

This is a new submission.

Both win-builder checks also flagged "sparsification" as a possibly
misspelled word in DESCRIPTION. This is an intentional technical
term describing row-wise post-estimation sparsification paths.

Examples, tests, vignettes and PDF and HTML manuals passed.
