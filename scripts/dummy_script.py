"""A simple script to serve as an example.

Should be deleted in real workflows.
"""

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

import utils

if TYPE_CHECKING:
    snakemake: Any

LOGGER = logging.getLogger(__name__)
utils.set_logger(LOGGER, snakemake.log[0])

config = snakemake.params.config_text
readme = Path(snakemake.input.readme).read_text()
user = Path(snakemake.input.user_file).read_text()

output_text = "\n\n".join([readme, user, config])

LOGGER.info("Hello world")

Path(snakemake.output.combined).write_text(output_text)
