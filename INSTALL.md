# Installation

## Using pip (recommended)

Install chaintools directly from PyPI:

```bash
pip install chaintools
```

## Using uv (fast alternative)

If you have [uv](https://docs.astral.sh/uv/) installed:

```bash
uv pip install chaintools
```

## From source

Clone the repository and install in development mode:

```bash
git clone https://github.com/milkschen/chaintools.git
cd chaintools
pip install -e .
```

Or with uv:

```bash
git clone https://github.com/milkschen/chaintools.git
cd chaintools
uv pip install -e .
```

Or use `uv sync` to install dependencies from the lock file:

```bash
git clone https://github.com/milkschen/chaintools.git
cd chaintools
uv sync
```

## Requirements

- Python ≥3.8
- pysam ≥0.15.3
- intervaltree ≥3.1.0
- typer ≥0.9.0
- pandas ≥0.24.1

Dependencies will be automatically installed with any of the above methods.
