SIMS = ["sigma", "gb", "anova", "16s"]
STUDIES = ["ob_goodrich", "cdi_schubert"]

configfile: "config.yaml"

wildcard_constraints:
    study="[^\\.]+"

rule all:
    input:
        expand("fig/{sim}.pdf", sim=SIMS),
        "results/min-effect-sizes.tsv"

rule clean:
    shell: "rm -rf fig/* cache/* data/* results/*"

rule plot_16s:
    output: "fig/16s.pdf"
    input: "plot-16s.R", "results/16s.tsv"
    shell: "./plot-16s.R"

rule table:
    output: "results/min-effect-sizes.tsv"
    input:
        files=expand("results/{sim}.tsv", sim=SIMS),
        script="find-min-effect-sizes.R"
    shell: "./find-min-effect-sizes.R {input.files}"

rule simulate:
    wildcard_constraints: x="(gb|sigma|anova)"
    output: "fig/{x}.pdf", "results/{x}.tsv"
    input: "simulate-{x}.R", "utils.R"
    shell: "./simulate-{wildcards.x}.R"

rule simulate_16s:
    output: "results/16s.tsv"
    input:
        expand("data/{study}_results/{study}.{file}", study=STUDIES, file=["otu_table.100.denovo", "metadata.txt"]),
        "simulate-16s.R", "utils.R"
    shell: "./simulate-16s.R"

rule extract:
    output: "data/{study}_results/{file}"
    input: "data/{study}_results.tar.gz"
    shell:
        "tar xvf {input}"
        " -C data/ {wildcards.study}_results/{wildcards.file}"

rule download:
    output: "data/{study}_results.tar.gz"
    params:
        md5=lambda wc: config["studies"][wc["study"]]["md5"],
        url=lambda wc, output: config["base_url"] + "/" + os.path.basename(output[0])
    shell: "./download-data.R {params.url} {params.md5} {output[0]}"
