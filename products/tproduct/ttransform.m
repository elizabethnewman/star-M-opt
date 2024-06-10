function[A] = ttransform(A, invFlag)

if ~exist('invFlag', 'var') || isempty(invFlag)
    invFlag = false;
end

if invFlag
    A = ifft(A,[],3);
else
    A = fft(A,[],3);
end

end