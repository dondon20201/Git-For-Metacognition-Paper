# What Generalizes in Metacognition?

Code to reproduce the analyses in "What Generalizes in Metacognition? Confidence Level, Informativeness, and a Measurement-Reliability Account of Apparent Dissociations" by Brandon Chow.

## How to run

1. Install R 4.6.1.
2. Open R in this folder and run `renv::restore()` to install the exact package versions.
3. From this folder, run the scripts in `analysis/` in order: `Rscript analysis/00_fetch_verify.R`, then 01, 02 and so on up to 17.

Script 00 downloads the Lund et al. (2025) data and script 11 downloads the Confidence Database (Rahnev et al., 2020). Both check every file's MD5 before use. Results are written to `results/` and figures to `figures/`.
