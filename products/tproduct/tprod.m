function[C] = tprod(A, B)

C = ttransform(facewise(ttransform(A),ttransform(B)),1);

end