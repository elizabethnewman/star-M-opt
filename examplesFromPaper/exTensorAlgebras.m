
clear; clc;



% setup 3 x 3 examples
n3 = 8;

% setup symbolic vector
a = [];
for i = 1:n3
    syms(['a',num2str(i)],'real')
    a = cat(1,a,eval(['a',num2str(i)]));
end

% algebraic structure
R = @(M) M \ (diag(M * a(:)) * M);

%% Matrices

I = eye(n3);
F = dftmtx(n3);         
S = triu(ones(n3));
D = eye(n3) - diag(ones(n3-1,1),1);

RI = diag(a);
RF = circulant(a);
RS = triu(a * ones(n3,1)',1) + diag(S * a);
RD = triu(ones(n3,1) * (D * a)') + triu(ones(n3,1) * (D' * a)',1) ;


sympref('FloatingPointOutput',true);

% experiment
M = {'I',I,RI,'F',F,RF,'S',S,RS,'D',D,RD};

for i = 1:3:length(M)

    fprintf('============== %s ============== \n',M{i})
    disp('M^{-1} * diag(M * a) * M')
    disp(simplify(R(M{i + 1})))

    disp('By Hand')
    disp(simplify(M{i + 2}))

    fprintf('Equal? %d\n',isequal(simplify(R(M{i + 1})),simplify(M{i + 2})))
    fprintf('=============================== \n')
end

sympref('FloatingPointOutput','default');


%% Helper Functions
function[C] = circulant(v)
% adapted from https://www.mathworks.com/matlabcentral/fileexchange/42-circulant-m
m   = length(v);

idx = (1:m)' * ones(1,m);
idx = mod(idx - idx',m) + 1;

C = v(idx);

end



