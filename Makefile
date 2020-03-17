TARGETS := fig/anova.pdf
ANOVA_SCRIPT := simulate-anova.R

all: $(TARGETS)

clean:
	rm -rf fig/*

fig/anova.pdf: $(ANOVA_SCRIPT)
	./$(ANOVA_SCRIPT)
