import logging
from typing import TYPE_CHECKING, Any

import pandas as pd
import utils

if TYPE_CHECKING:
    snakemake: Any

LOGGER = logging.getLogger(__name__)
utils.set_logger(LOGGER, snakemake.log[0])


def multiply_and_add(
    df: pd.DataFrame, multiplier: float, adder: float, apply_to_columns: list[str]
) -> pd.DataFrame:
    df_new = df.copy()

    df_new[apply_to_columns] *= multiplier
    df_new[apply_to_columns] += adder
    return df_new


if __name__ == "__main__":
    breakpoint()
    df = pd.read_csv(snakemake.input.csv_to_update)

    df_new = multiply_and_add(
        df,
        snakemake.params.multiplier,
        snakemake.params.adder,
        snakemake.params.apply_to_columns,
    )

    df_new.to_csv(snakemake.output.final_csv)
