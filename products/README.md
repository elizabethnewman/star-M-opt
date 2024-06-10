# Tensor-Tensor Products

The structure of every product looks like the following:

```
function[output,Jac1,Jac2] = tensorProduct(input1,input2)
```

where ```Jac1``` is the Jacobian of the output with respect to the first input, ```Jac2``` with respect to the second input, and so on. 

Currently, each Jacobian is a struct with four fields,

```
Jac.inputSize = size(input)
Jac.outputSize = size(output)
Jac.A  = @(x) ...
Jac.AT = @(y) ...
```

The first function handle ```Jac.A``` is a mapping from ```size(input)``` to ```size(output)```.  The adjoint, ```Jac.AT``` is a mapping from ```size(output)``` to ```size(input)```.
