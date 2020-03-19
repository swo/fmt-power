TARGETS := fig/sigma.pdf fig/gb.pdf fig/anova.pdf fig/schubert.pdf fig/sigma-spread.pdf results/min-effect-sizes.tsv
SIM_RESULTS := results/sigma.tsv results/gb.tsv results/anova.tsv results/schubert.tsv
SIGMA_SCRIPT := simulate-sigma.R
GB_SCRIPT := simulate-gb.R
ANOVA_SCRIPT := simulate-anova.R
SPREAD_SCRIPT := plot-sigma-spread.R
SCHUBERT_SCRIPT := simulate-schubert.R
MINEF_SCRIPT := find-min-effect-sizes.R

all: $(TARGETS)

clean:
	rm -rf fig/* cache/*

fig/sigma.pdf results/sigma.tsv: $(SIGMA_SCRIPT) utils.R
	./$(SIGMA_SCRIPT)

fig/gb.pdf results/gb.tsv: $(GB_SCRIPT) utils.R
	./$(GB_SCRIPT)

fig/anova.pdf results/anova.tsv: $(ANOVA_SCRIPT) utils.R
	./$(ANOVA_SCRIPT)

fig/schubert.pdf results/schubert.tsv: $(SCHUBERT_SCRIPT) utils.R
	./$(SCHUBERT_SCRIPT)

fig/sigma-spread.pdf: $(SPREAD_SCRIPT)
	./$(SPREAD_SCRIPT)

results/min-effect-sizes.tsv: $(SIM_RESULTS) $(MINEF_SCRIPT)
	./$(MINEF_SCRIPT)
