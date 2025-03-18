rule dummy_add_text:
    message:
        "Dummy rule combining user inputs and automatic downloads."
    params:
        config_text=config["dummy_text"],
    input:
        user_file="resources/user/user_message.md",
        readme="resources/automatic/dummy_readme.md",
    output:
        combined="results/combined_text.md",
    conda:
        "../envs/shell.yaml"
    log:
        "logs/dummy_add_text.log"  # relative to calling file
    script:
        "../scripts/dummy_script.py"


rule multiply_csv:
    message: "multiply my CSV file by a number"
    input:
        csv_to_update = "resources/user/data.csv"
    params:
        multiplier = config["parameters"]["multiplier"],
        adder = config["parameters"]["adder"],
        apply_to_columns = config["parameters"]["apply_to_columns"],
    output:
        final_csv = "results/final_data.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/data_manipulator.py"
