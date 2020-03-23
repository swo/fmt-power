FIGS = ["sigma", "gb", "anova", "16s", "sigma-spread"]
SIMS_BASIC = ["sigma", "gb", "anova"]
SIMS_16S = ["ob_goodrich", "cdi_schubert"]

configfile: "config.yaml"

wildcard_constraints:
    study="[^\\.]+"

rule all:
    input:
        expand("fig/{fig}.pdf", fig=FIGS),
        "results/min-effect-sizes.tsv"

rule clean:
    shell: "rm -rf fig/* cache/* data/*"

rule x16s_plot:
    output: "fig/16s.pdf"
    input: "plot-16s.R", expand("results/{sim}.tsv", sim=SIMS_16S)
    shell: "./plot-16s.R"

rule sigma_plot:
    output: "fig/sigma-spread.R"
    input: "plot-sigma-spread.R"
    shell: "./plot-sigma-spread.R"

rule table:
    output: "results/min-effect-sizes.tsv"
    input:
        files=expand("results/{sim}.tsv", sim=SIMS_BASIC + SIMS_16S),
        script="find-min-effect-sizes.R"
    shell: "./find-min-effect-sizes.R {input.files}"

rule analyze:
    output: "fig/{x}.pdf", "results/{x}.tsv"
    input: "simulate-{x}.R", "utils.R"
    shell: "./simulate-{wildcards.x}.R"

rule analyze_16s:
    output: "results/16s.tsv"
    input:
        expand("data/{study}_results/{study}.file}", study=SIMS_16S, faile = ["otu_table.100.denovo", "metadata.txt"]),
        "simulate-16s.R", "utils.R"
    shell: "./simulate-16s.R"

rule extract:
    output: "data/{study}_results/{file}"
    input: "data/{study}_results.tar.gz"
    shell:
        "tar xvf {input}"
        " -C raw/ {wildcards.study}_results/{wildcards.file}"

rule download:
    output: "data/{study}_results.tar.gz"
    params:
        md5=lambda wc: config["studies"][wc["study"]]["md5"],
        url=lambda wc, output: config["base_url"] + "/" + os.path.basename(output[0])
    run: "./download-data.R {params.url} {params.md5} {output[0]}"
