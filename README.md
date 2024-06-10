# tenOpt
Optimal tensor-tensor products
This repository contains code for computing optimal tensor-tensor products using the star-M product introduced in []. 

## Installation

Clone the repository using
```
git clone https://github.com/elizabethnewman/tenOpt.git
```

## Setup

To setup the paths, open MATLAB, make ```tenOpt``` the working directory, and run ```tenOptSetup.m``` from the command line.


## Required Matlab Toolboxes

* Financial Toolbox, Datafeed Toolbox
* PDE Toolbox

## Organization:

- **optimizers:** Options for optimization algorithms (gradient descent, alternating descent).
- **objectiveFunctions:** Options for problem-specific objective functions (least squares, low-rank approximation) to pass to optimizer.
- **linesearch:** Options for line search algorithms for optimizer.
- **products:** Functions for tensor-tensor (star-M, facewise) and tensor-matrix (mode-k) products.
- **examples:** Examples applying optimal tensor-tensor products to various applications on small datasets.
- **unitTests:** unit tests for other functions in repository.
- **utils:** additional tools for tensor operations (folding/unfolding, Frobenius norm, transpose).
- **tutorials:** introductory tutorials to demonstrate code functions.

## Introductory Materials

To illustrate the functions available in this repository, we have provided some [tutorials](tutorials) formatted as .mlx live scripts.

## How to cite

```
@{

}
```

## Acknowledgements
