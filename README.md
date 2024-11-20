This repository contains public code for the manuscript:

Zou, Heng-Xing, Xinyi Yan, and Volker HW Rudolf. 2024. "Time-dependent Interaction Modification Generated from Plant-soil Feedback." *Ecology Letters*. [https://onlinelibrary.wiley.com/doi/10.1111/ele.14432]([https://doi.org/10.1101/2023.11.03.565336](https://onlinelibrary.wiley.com/doi/10.1111/ele.14432))

For any questions, please [contact Heng-Xing Zou](hxzou.ecology@gmail.edu)

# Abstract

Pairwise interactions between species can be modified by other community members, leading to emergent dynamics contingent on community composition. Despite the prevalence of such higher-order interactions, little is known about how they are linked to the timing and order of species’ arrival. We generate population dynamics from a mechanistic plant-soil feedback model, then apply a general theoretical framework to show that the modification of a pairwise interaction by a third plant depends on its germination phenology. These time-dependent interaction modifications emerge from concurrent changes in plant and microbe populations and are strengthened by higher overlap between plants’ associated microbiomes. The interaction between this overlap and the specificity of microbiomes further determines plant coexistence. Our framework is widely applicable to mechanisms in other systems from which similar time-dependent interaction modifications can emerge, highlighting the need to integrate temporal shifts of species interactions to predict the emergent dynamics of natural communities.

# Code Files

All simulations were run in R version 4.2.1.

- `Functions_PSF.R` contains all functions used for contructing the plant-soil feedback model.
- `GenerateData.R` contains code for generating data using the mechanistic plant-soil feedback model.
- `GenerateData_MultipleParams.R` contains code for generating data using randomly selected parameters; this is required for the sensitivity analysis (see Appendix I).
- `ExamineModels.R` contains code for comparing the performance of several dynamic models and generate Figure S1.
- `FeedbackRatios.R` contains code for calculating the ratio of intra- vs. interspecific effect from microbiomes using empirical data (Yan et al. 2022 [_PNAS_](https://www.pnas.org/doi/abs/10.1073/pnas.2122088119)).
- `ModelFitting.R` contains all code for fitting the simulated data to phenomenological models of higher-order interactions.
- `Analysis_Tidy.Rmd` contains code for creating all figures in the main text and Appendix.

# Data

The `Data` directory contains six subdirectories that are categorized below. 

The name of raw data (i.e., data simulated from the plant-soil feedback model) has three components:
- `ij` or `ijk`: data generated from two- or three-plant models;
- `al`, `nc`, or `co`: representing different degrees of overlap between microbiomes; `al` represents low overlap, `nc` represents no overlap, and `co` represents high overlap (see also Figure 1);
- `ts` or `rs`: `ts` represents data generated from a time series, `rs` represents data generated by a response surface design for one year (see Methods).

All data are `.csv` files.

## Population Dynamics

- `Standard` is the default directory that contains population dynamics used to produce all results in the main text.
- `m_ii_small` contains population dynamics with a set of parameters that specifies a smaller intraspecific microbial effect.
- `m_ij_small` contains population dynamics with a set of parameters that specifies a larger intraspecific microbial effect.

`Group` denotes the type of population (plants, seeds, or microbiomes), and `init1`, `init2`, and `init3` denotes the initial population of plants, seeds, or microbiomes  $i$, $j$, and $k$. `pi`, `pj`, and `pk` denotes germination phenology of plants $i$, $j$, and $k$; `pij` is the phenological difference calculated by $p_j-p_i$ (see Methods and Figure 1). For datasets ending in `rs` (response surface), `Time` ranges from 1 to 12 because we only data of one year to estimate the interaction coefficients (see Methods). 

## Interaction Coefficients

- `Repetition` contains *fitted interaction coefficients* and their standard error reported by the model from simulations with random parameters. This set of simulations is only performed with three-plant models.
- `RealisticProp` contains *fitted interaction coefficients* and their standard error reported by the model from simulations with random parameters drawn from distribution of ratios of intra- vs. interspecific effect from microbiomes using empirical data (see below). This set of simulations is only performed with three-plant models.

Coefficients starting with `alpha` are pairwise, coefficients starting with `beta` are higher-order (see Methods). `pi`, `pj`, and `pk` denotes germination phenology of plants $i$, $j$, and $k$; `pij` is the phenological difference calculated by $p_j-p_i$ (see Methods and Figure 1). `rep` denotes the repetition (from 1 to 100). 

## Empirical Data

- `Metaanalysis_Yan` contains data used to calculate the ratio of intra- vs. interspecific effect from microbiomes using empirical data.

See [Yan et al. 2022](https://doi.org/10.5281/zenodo.6513066) for the original data and formatting.


