TARGETS := fig/sigma.pdf fig/gb.pdf fig/anova.pdf fig/schubert.pdf fig/sigma-spread.pdf
SIGMA_SCRIPT := simulate-sigma.R
GB_SCRIPT := simulate-gb.R
ANOVA_SCRIPT := simulate-anova.R
SPREAD_SCRIPT := plot-sigma-spread.R
SCHUBERT_SCRIPT := simulate-schubert.R

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
