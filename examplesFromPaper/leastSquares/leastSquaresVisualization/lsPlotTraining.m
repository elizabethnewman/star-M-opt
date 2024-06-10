

function[fig] = lsPlotTraining(results,iter)

[A,B] = lsSetupData(results.options);

fctn = objFctnTikhonovLeastSquaresVarPro(A,B,results.options.lambda,results.options.alpha);

M = reshape(results.optInfo.x(:,iter),size(A,3),[]);
X = fctn.solve(M);

fig = lsPlotHyperplanes(A, B, X, M, 1, 0);

end

