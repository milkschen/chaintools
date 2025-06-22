import sys
import logging
from importlib.metadata import version, PackageNotFoundError

if sys.version_info < (3, 8):
    logging.warning(
        "chaintools will require Python 3.8 or higher in the next release. "
        "Current version: %s",
        sys.version,
    )
    # raise RuntimeError(
    #     "chaintools requires Python 3.8 or higher. "
    #     f"Current version: {sys.version}"
    # )

try:
    __version__ = version("chaintools")
except PackageNotFoundError:
    # package is not installed
    pass
