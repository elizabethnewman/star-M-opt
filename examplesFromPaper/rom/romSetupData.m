function[A] = romSetupData(options)


A = [];

% choose parameters

cList = linspace(options.cMin,options.cMax,options.cN);

% choose time points 
tList = linspace(options.tMin,options.tMax,options.tN);

% solve for each parameter
for c = cList
    
    % create model
    model = romSetupWaveEquation();

    specifyCoefficients(model,"m",1,...
                          "d",0,...
                          "c",c,...
                          "a",0,...
                          "f",0);

    % solve
    results = solvepde(model,tList);

    % store results
    A = cat(3,A,results.NodalSolution);
end

end

