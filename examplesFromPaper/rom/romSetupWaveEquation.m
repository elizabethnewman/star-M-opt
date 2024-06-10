function[model] = romSetupWaveEquation()

% setup PDE model
model = createpde();

% create square boundary
geometryFromEdges(model,@squareg);

% boundary conditions
applyBoundaryCondition(model,"dirichlet","Edge",[2,4],"u",0);
applyBoundaryCondition(model,"neumann","Edge",[1 3],"g",0);

% initial conditions
u0  = @(location) atan(cos((pi/2) * location.x));
ut0 = @(location) 3 * sin(pi * location.x) .* exp(sin((pi/2) * location.y));
setInitialConditions(model,u0,ut0);

% create mesh
generateMesh(model,'GeometricOrder','linear','Hmax',0.1);

end