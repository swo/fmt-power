# Power calculations for detecting differences in efficacy of fecal microbiota transplantation donors

## Dependencies

- \*nix environment (for the downloading and extracting commands to work)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (not essential, but removes the requirement to run scripts in order by hand)
- R (developed against 3.6.0), including packages:
    - [tidyverse](https://tidyverse.tidyverse.org/)
    - [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
    - [memoise](https://cran.r-project.org/web/packages/memoise/index.html) (again, not strictly essential, but made development much faster)

## To do

## Outline

### Model donors as Gaussian LO

- Fig 1a: Distributions of donor efficacies
- 1b: Show drawing from that distribution, getting yes/no, counting up
- 1c: What is the critical σ to detect with FFH @ 80% power?

### Model donors as good/bad

- Fig 2: Show distribution
- 2b: Same as 1b?
- 2c: Same table

### Good/bad donors with microbiome

- 3a: Same as 2a
- 3b: Assign good/bad to random 16S's, simulate, test for separation
- 3c: Same as 2c; with diarrhea, critical Δp
- 3d: Same as 3c; with weak m.b. signal (IBD?)

### Biomarker

- 3a: Donors have some mean values, drawn from donor envelope
- 3b: Patients have values drawn from envelopes around each donor
- 3c: Table of values; looking for critical σ_donor / σ_patient

Look in lit for σ_patient across IBD patients. Look at σ_donor as saying a
strong drug's effect?
