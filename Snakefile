SIMS = ["mannwhitney", "contingency", "anova", "16s"]
STUDIES = ["ob_goodrich", "cdi_schubert"]

configfile: "config.yaml"

rule all:
    input: "results/min-effect-sizes.tsv"

rule clean:
    shell: "rm -rf cache/* data/* results/*"

rule table:
    output: "results/min-effect-sizes.tsv"
    input:
        files=expand("results/{sim}.tsv", sim=SIMS),
        script="find-min-effect-sizes.R"
    shell: "./{input.script} {input.files}"

rule simulate_16s:
    output: "results/16s.tsv"
    input:
        expand("data/{study}_results/{study}.{file}", study=STUDIES, file=["otu_table.100.denovo", "metadata.txt"]),
        "utils.R", script="simulate-16s.R"
    shell: "./{input.script}"

rule extract:
    wildcard_constraints: study="[^\\.]+"
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

rule simulate:
    output: "results/{x}.tsv"
    input: "utils.R", script="simulate-{x}.R"
    shell: "./{input.script}"

