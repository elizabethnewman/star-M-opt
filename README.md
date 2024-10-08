# star-M-opt
This repository contains code for learning optimal tensor-tensor products using the star-M product introduced in [Tensor–tensor products with invertible linear transforms](https://www.sciencedirect.com/science/article/pii/S0024379515004358). 

## Installation

Clone the repository using
```
git clone https://github.com/elizabethnewman/star-M-opt.git
```

## Setup

To setup the paths, open MATLAB, make ```star-M-opt``` the working directory, and run 
```
starMOptSetup
```
in the MATLAB console.


## Required Matlab Toolboxes

* Financial Toolbox, Datafeed Toolbox
* PDE Toolbox

## Organization:

- **optimizers:** options for optimization algorithms (gradient descent, alternating descent).
- **objectiveFunctions:** options for problem-specific objective functions (least squares, low-rank approximation) to pass to optimizer.
- **products:** functions for tensor-tensor (star-M, facewise) and tensor-matrix (mode-k) products.
- **tensorSVD:** functions to compute the tensor SVD and corresponding Jacobians.
- **examples:** examples applying optimal tensor-tensor products to various applications.
- **utils:** additional tools for tensor operations (folding/unfolding, Frobenius norm, transpose).
- **unitTests:** unit tests for other functions in repository.  To test, run ```starMOptUnitTests``` in a MATLAB console.
- **tutorials:** introductory tutorials to demonstrate code functions.

## Introductory Materials

To illustrate the functions available in this repository, we have provided some [tutorials](https://github.com/elizabethnewman/star-M-opt/tree/main/tutorials) formatted as .mlx live scripts.

## How to cite

```
@misc{newman2024optimalmatrixmimetictensoralgebras,
      title={Optimal Matrix-Mimetic Tensor Algebras via Variable Projection}, 
      author={Elizabeth Newman and Katherine Keegan},
      year={2024},
      eprint={2406.06942},
      archivePrefix={arXiv},
      primaryClass={math.NA},
      url={https://arxiv.org/abs/2406.06942}, 
}
```

