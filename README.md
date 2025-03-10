[![DOI](https://zenodo.org/badge/816474443.svg)](https://zenodo.org/doi/10.5281/zenodo.13386029)

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

### Examples
An example input file can be found in the ```EXAMPLES``` directory. This file contains a pair of edge-weight phylogenetic trees over the same taxa. Below are the flags that can be used:

```bash
Usage: nlcd dist --norm <NORM> --input_file <FILE_PATH> --method <METHOD> --weighted <WEIGHTED>

Options:
  -p, --norm <NORM>             nth norm
  -i, --input_file <FILE_PATH>  Input tree file in Newick format
  -m, --method <METHOD>         One of size, depth, height, or local [default: depth]
  -w, --weighted <WEIGHTED>     Use edge weights [default: false] [possible values: true, false]
  -h, --help                    Print help
```
To compute the cophenetic distance between the trees taking only topology into account, set weighted to false, which will set all edge-weights to 1. When method is set to ```depth```, the depth of the vertex will be computed as the sum of all the edge-weights in the path from the vertex to the root. Similarly, When method is set to ```height```, the height of a vertex will be computed as the sum of all the edge-weights in the path from the vertex to the closest leaf. When method is set to ```size```, the size of the vertex will be computed as the number of leaves in the subtree rooted at that vertex.

For example, the following command computes the distance between the pair of example trees by comparing the depth of the vertices under the first norm.
```bash
nlcd dist -i ./EXAMPLES/trees.nwk -p 1 -m depth -w false
```

Similarly, the following command computes the distance between the pair of example trees by comparing the weighted distances of the vertices from the root vertex under the first norm.
```bash
nlcd dist -i ./EXAMPLES/trees.nwk -p 1 -m depth -w true
```

In the case that the contributing function is externally computed, you can set the weight of the edge between a vertex and it's parent to be the value of the contributing function for that vertex and then run the folloing command.
```bash
nlcd dist -i ./EXAMPLES/trees.nwk -p 1 -m local -w true
```

Note that when mode is set to ```local``` the weighted flag is treated ```true``` regardless of user input, and the weight of the edge leading to the root vertex is always zero. 