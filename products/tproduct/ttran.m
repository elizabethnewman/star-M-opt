function[A] = ttran(A)

A = permute(A,[2,1,3]);
A(:,:,2:end) = A(:,:,size(A,3):-1:2);

end