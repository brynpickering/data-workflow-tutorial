"""Rules to used to download automatic resource files."""


rule dummy_download:
    message:
        "Download the clio README file."
    params:
        url=config["resources"]["automatic"]["dummy_readme"],
    output:
        readme="resources/automatic/dummy_readme.md",
    conda:
        "../envs/shell.yaml"
    shell:
        'curl -sSLo {output.readme} "{params.url}"'
