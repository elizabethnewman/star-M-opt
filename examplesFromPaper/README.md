## Organization

The examples are from the paper
```console
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


To replicate the experiment results in the paper, run the following from the MATLAB console:
* Example 3.4: Lipschitz Continuity Convergence Gaurantees
```
LipschitzConvergencePaperResults
```
* Example 3.5: Algorithm Comparison
```
convergencePaperResults
``` 
* Section 4.2: Index Tracking
```
indexTrackingPaperResults
```
*  Section 4.3: Digits
```
digitsPaperResults
```
*  Section 4.4: Reduced Order Modeling
```
romPaperResults
```

For the least squares, index tracking, and reduced order modeling examples, the corresponding folders have the following four MATLAB functions:
```
<name>ExperimentParameters.m
<name>Run.m
<name>SetupData.m
```
The experiment parameters are based on the object
```
starMOptExperimentParameters.m
```
Visualization tools are included in corresponding folders. The high resolution images in the paper were made offline using Tikz and PGFPlots from the stored results. 
