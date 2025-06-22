# -*- coding: utf-8 -*-
import sys
import logging
from importlib.metadata import version, PackageNotFoundError

if sys.version_info < (3, 8):
    logging.warning(
        "chaintools_bio will require Python 3.8 or higher in the next release. "
        "Your current Python version is {sys.version_info.major}.{sys.version_info.minor}."
    )
    # raise ImportError(
    #     "chaintools requires Python 3.8 or higher. "
    #     f"Current version: {sys.version}"
    # )

try:
    __version__ = version("chaintools_bio")
except PackageNotFoundError:
    # package is not installed
    pass
