import yaml

from snakemake.utils import min_version

min_version("8.10")


configfile: "config/default.yaml"


# Add all your includes here.
include: "rules/automatic.smk"
include: "rules/dummy.smk"


rule all:
    message:
        "Run project."
    input:
        "results/combined_text.md",
    default_target: True
