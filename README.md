# MLM
**M**ulti **L**evel **M**igration is a pipeline for ancestral status inference of nodes from phylogeographic analysis, proudly introduced by Naclist.

## Table of Contents
1. [Installation](#installation)
2. [Dependencies](#dependencies)
3. [Usage](#usage)
    - [MLM.py Usage](#mlmpy-usage)
    - [Treetime Migration Usage](#treetime-migration-usage)
4. [License](#license)

---

## Installation

To use MLM.py, you need to have Python 3.x installed. The script relies on several third-party packages, which can be installed using `pip`.

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/MLM.py.git
    cd MLM.py
    ```

2. Create a virtual environment (optional but recommended):
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows, use venv\Scripts\activate
    ```

3. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

---

## Dependencies

The following Python packages are required to run `MLM.py`:

- `biopython`: For working with phylogenetic trees (Newick format).
- `pandas`: For handling confidence matrices in tabular format.
- `treetime`: For applying the Treetime migration model.
- `numpy` (optional, for numerical operations).

Install the dependencies with:

```bash
pip install biopython pandas treetime numpy

