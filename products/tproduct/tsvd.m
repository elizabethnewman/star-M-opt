function[U,S,V] = tsvd(A, k)

if ~exist('k','var'), k = []; end

[UHat,SHat,VHat] = facewiseSVD(ttransform(A), k);
U = ttransform(UHat,1);
S = ttransform(SHat,1);
V = ttransform(VHat,1);


end