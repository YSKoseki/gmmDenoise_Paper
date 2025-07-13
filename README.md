# gmmDenoise_Paper

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/540440691.svg)](https://zenodo.org/badge/latestdoi/540440691)
<!-- badges: end -->

## Overview
This repository includes several directories containing the R code used in the analyses presented in Koseki et al. (2025). Below is a brief overview of each directory:

- `Example_1/`
Contains scripts for analyzing stream fish community eDNA data from Nakagawa et al. (2018). The analysis serves as an illustrative example of how `gmmDenoise` can be used to derive accurate intraspecific diversity estimates and population genetic inferences.

- `Example_2/`
Contains scripts for another example analysis using `gmmDenoise`. This example uses estuarine fish community eDNA data from Ahn et al. (2020).

- `Example_comp_fig/`
Contains a script to create a set of histograms illustrating the effect of the stringency-controlling parameter of `UNOISE3` on the read count distributions of amplicon sequence variants (ASVs) obtained in `Example_1` and `Example_2`.

- `PCR_Artifact_Sim/`
Includes scripts for simulating eDNA metabarcoding processes, including PCR, to analyze the read count distributions of true sequences and artifacts derived from them.

- `gmmDenoise_Test_1/`
Contains scripts used to evaluate the performance of `gmmDenoise` with mock eDNA data of _Plecoglossus altivelis_ from Tsuji, Miya et al. (2020).

- `gmmDenoise_Test_2/`
Contains scripts for another performance test of `gmmDenoise` using stream eDNA data of _P. altivelis_ from Tsuji, Maruyama et al. (2020).

Each script in these directories includes a brief header and inline comments explaining its purpose and usage. The core `gmmDenoise` package is available at the following repository: [https://github.com/YSKoseki/gmmDenoise](https://github.com/YSKoseki/gmmDenoise).

## References
- Ahn, H., Kume, M., Terashima, Y., Ye, F., Kameyama, S., Miya, M., Yamashita, Y., & Kasai, A. (2020). Evaluation of fish biodiversity in estuaries using environmental DNA metabarcoding. PLOS ONE, 15(10), e0231127.
- Koseki, Y., Takeshima, H., Yoneda, R., Katayanagi, K., Ito, G., & Yamanaka, H. (2025). _gmmDenoise_: a new method and _R_ package for high-confidence sequence variant filtering in environmental DNA amplicon analysis. Authorea.
- Nakagawa, H., Yamamoto, S., Sato, Y., Sado, T., Minamoto, T., & Miya, M. (2018). Comparing local‐ and regional‐scale estimations of the diversity of stream fish using eDNA metabarcoding and conventional observation methods. Freshwater Biology, 63(6), 569–580.
- Tsuji, S., Maruyama, A., Miya, M., Ushio, M., Sato, H., Minamoto, T., & Yamanaka, H. (2020). Environmental DNA analysis shows high potential as a tool for estimating intraspecific genetic diversity in a wild fish population. Molecular Ecology Resources, 20(5), 1248–1258.
- Tsuji, S., Miya, M., Ushio, M., Sato, H., Minamoto, T., & Yamanaka, H. (2020). Evaluating intraspecific genetic diversity using environmental DNA and denoising approach: A case study using tank water. Environmental DNA, 2(1), 42–52.

## Funding
This work was supported by JSPS KAKENHI Grant Numbers JP21K12329 and JP22K14908.
