"""Utility Functions."""

import logging
import sys


def set_logger(logger: logging.Logger, logfile: str | None = None):
    """Set up the synthpop loggers with a specific format.

    Args:
        logger (logging.Logger): Logger to set up (e.g. linked to a specific source code file).
        logfile (str | None, optional):
            If given, add a logging file handler, to pipe logs to that file.
            Logs will still be piped to STDOUT if this is given.
            Defaults to None.
    """
    formatter = logging.Formatter(
        "[%(asctime)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    console = logging.StreamHandler(stream=sys.stdout)
    console.setFormatter(formatter)
    logger.addHandler(console)
    logger.setLevel(logging.INFO)

    if logfile is not None:
        file_handler = logging.FileHandler(logfile)
        file_handler.setLevel(5)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
