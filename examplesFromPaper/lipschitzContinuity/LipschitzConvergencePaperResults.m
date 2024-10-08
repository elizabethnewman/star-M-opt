
clear; clc;

% setup paths
run("../../starMOptSetup.m")

% setup directory
dirName    = 'LipschitzConvergenceResults/';

if ~exist(dirName,'dir'), mkdir(dirName); end
cd(dirName)

subDirName = ['./',datestr((datetime('today'))),'/'];
if ~exist(subDirName,'dir'), mkdir(subDirName); end
cd(subDirName)

% return to main directory
cd ../..

%% setup data

% create tensor
A(:,:,1) = [1,0;0,1;0,0];
A(:,:,2) = [0,0;1,0;0,1];

B(:,:,1) = [1;1;1];
B(:,:,2) = [1;1;1];

X(:,:,1) = [1;1];
X(:,:,2) = [1;1];

% Givens rotation matrix
Q = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];

% objective function
fctn = objFctnTikhonovLeastSquaresVarPro(A,B,0,0,'avgFlag',0);

%% optimize, starting at different angles

LStar = 1e-1;
for i = 0:359
    thetai = pi / 180 * i;

    opt = gradientDescent('verbose',1,'maxIter',500,'manifold','Stiefel','logInterval',1,'storeIterates',1,'verbose',0);
    opt.linesearch = constantLinesearch('alpha',LStar,'manifold','Stiefel','gamma',1e-3,'maxIter',40);
   
    % solve
    [MSol,optInfo] = opt.solve(fctn,Q(thetai));
    
    % get path and quadrants
    thetaPath = [];
    quadPath  = [];
    for j = 1:size(optInfo.x,2)
        c     = optInfo.x(1,j);
        s     = optInfo.x(2,j);
        [t,q] = findQuadrant(c,s);
        thetaPath = cat(2,thetaPath,t);
        quadPath  = cat(2,quadPath,q);
    end

    T = array2table([optInfo.values,thetaPath(:),quadPath(:)], 'VariableNames',cat(2,optInfo.header,'thetaPath','quadPath'));
    writetable(T,sprintf([dirName,subDirName,'init_%0.3d.csv'],i))
end


%% helper functions
function[t,q] = findQuadrant(c,s)
    t     = atan(abs(s / c));
    signC = sign(c);
    signS = sign(s);
    
    q = 1;
    if signC < 0
        if signS > 0
            % quad II
            t = pi - t;
            q = 2;
        else
            % quad III
            t = pi + t;
            q = 3;
        end
    else
        if signS < 0
            % quad IV
            t = 2 * pi - t;
            q = 4;
        end
            
    end
end