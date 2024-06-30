# Near Linear Cophenetic Distance
This repository contains a rust crate to compute the cophenetic distance between two rooted phylogenetic trees in near-linear time.

## Installation

To install you must first have cargo and rustup installed:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After installing the above command you can run the following to install nlcd:
```bash
cargo install --git=https://github.com/sriram98v/near-linear-cophenetic-distance
```

Alternatively, you can install seq_class by cloning this repository and building it locally:
```bash
git clone https://github.com/sriram98v/near-linear-cophenetic-distance
cd near-linear-cophenetic-distance
cargo install --path=./
```

## Usage
### Finding the Cophenetic distance between a pair of trees
To compute the cophenetic distance between a pair of trees, please create a single file with the extension ```.tre``` containing the two trees in Newick format (line-separated). The run the following command to compute the cophenetic distance with depth as the path function:
```bash
nlcd dist -i <PATH TO .TRE FILE> -p <NORM>
```

Please refer the help page for details on how to use other path functions using:
```bash
nlcd -h
```

### Reproduce empirical results
In order to reproduce the results as seen in the article, run the following command
```bash
nlcd repr-emp -n 1000 -x 100 -p 10 -t 1 -o ./emp-study
```

In order to plot the results, create a local python virtual environment and install the dependencies using the following commands:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The command above will create a file containing the distributions as seen in the Results section of the article. In order to reproduce the plots, please run the python script provided in the ```scripts``` directory as follows (run from the base of the repository):
```bash
./scripts/plot-distribs.py
```

### Reproduce scalability analysis
**Note: Reproduction of these results may not looked identical to that in the aritcle**
In order to reproduce the scalability analysis as seen in the article, run the following command
```bash
nlcd repr-sca -s 200 -e 10000 -x 200 -i 20 -p 5 -o ./sca-study
```

In order to plot the results, create a local python virtual environment and install the dependencies using the following commands:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The command above will create a file containing the distributions as seen in the Results section of the article. In order to reproduce the plots, please run the python script provided in the ```scripts``` directory as follows (run from the base of the repository):
```bash
./scripts/plot-sca.py
```