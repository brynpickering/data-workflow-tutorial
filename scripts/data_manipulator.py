import pandas as pd

if __name__ == "__main__":
    breakpoint()
    df = pd.read_csv(snakemake.input.csv_to_update)
    df_new = df.copy()

    df_new[snakemake.params.apply_to_columns] *= snakemake.params.multiplier
    df_new[snakemake.params.apply_to_columns] += snakemake.params.adder

    df_new.to_csv(snakemake.output.final_csv)
