TARGETS := fig/sigma.pdf fig/gb.pdf fig/anova.pdf fig/sigma-spread.pdf
SIGMA_SCRIPT := simulate-sigma.R
GB_SCRIPT := simulate-gb.R
ANOVA_SCRIPT := simulate-anova.R
SPREAD_SCRIPT := plot-sigma-spread.R

all: $(TARGETS)

clean:
	rm -rf fig/* cache/*

fig/sigma.pdf: $(SIGMA_SCRIPT) utils.R
	./$(SIGMA_SCRIPT)

fig/gb.pdf: $(GB_SCRIPT) utils.R
	./$(GB_SCRIPT)

fig/anova.pdf: $(ANOVA_SCRIPT) utils.R
	./$(ANOVA_SCRIPT)

fig/sigma-spread.pdf: $(SPREAD_SCRIPT)
	./$(SPREAD_SCRIPT)
