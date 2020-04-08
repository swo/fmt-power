# Power calculations for detecting differences in efficacy of fecal microbiota transplantation donors

## Getting started

Dependencies:

- \*nix environment (for the downloading and extracting commands to work)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (not essential, but removes the requirement to run scripts in order by hand)
- R (developed against 3.6.0), including packages:
    - [tidyverse](https://tidyverse.tidyverse.org/)
    - [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
    - [memoise](https://cran.r-project.org/web/packages/memoise/index.html) (again, not strictly essential, but made development much faster)

## File structure

- `README.md`: This file
- `Snakefile`: Workflow file
- `config.yaml`: Helper file for the Snakefile
- `cache/`: Results of simulations are cached here
- `data/`: Microbiome data is downloaded here
- `download-data.R`: Script to download microbiome data
- `simulate-XXX.R`: Scripts to run each of the simulations
- `utils.R`: Helper file for the simulation scripts
- `find-min-effect-sizes.R`: Final analysis script
- `results/`: Results of the simulations

The abbreviation `gb` (for "good-bad") is short-hand for the model with two
donor efficacies and dichotomous patient outcomes.

## Author

Scott Olesen <solesen@openbiome.org>
