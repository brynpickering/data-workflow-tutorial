import logging
from typing import TYPE_CHECKING, Any

import pandas as pd
import utils
import pandera as pa
if TYPE_CHECKING:
    snakemake: Any

LOGGER = logging.getLogger(__name__)
utils.set_logger(LOGGER, snakemake.log[0])


def multiply_and_add(
    df: pd.DataFrame, multiplier: float, adder: float, apply_to_columns: list[str]
) -> pa.DataFrameSchema:
    """Apply arithmetic to specific columns of a data table.

    Args:
        df (pd.DataFrame): Table containing pertinent columns.
        multiplier (float): Multiplier to apply to column data.
        adder (float): Number to add to column data _after_ applying the multiplier.
        apply_to_columns (list[str]): Columns on which to apply the arithmetic.

    Returns:
        pd.DataFrame: `df` with arithmetic applied to specific columns.
    """
    df_new = df.copy()

    df_new[apply_to_columns] *= multiplier
    df_new[apply_to_columns] += adder

    # Apply data schema
    validate_df(df_new)

    return df_new

def validate_df(df: pd.DataFrame):
    """Schema to validate the manipulated data table.

    Args:
        df (pd.DataFrame): Table to check

    Returns:
        pa.DataFrameSchema: Schema
    """
    schema = pa.DataFrameSchema({
        "population": pa.Column(int, pa.Check.ge(0)),
        "demand_heat_kwh": pa.Column(float, pa.Check.ge(0)),
        "demand_electricity_kwh": pa.Column(float, pa.Check.ge(0))
    },
    strict=True, coerce=True)

    schema.validate(df)

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input.csv_to_update, index_col=0)

    df_new = multiply_and_add(
        df,
        snakemake.params.multiplier,
        snakemake.params.adder,
        snakemake.params.apply_to_columns,
    )

    df_new.to_csv(snakemake.output.final_csv)
