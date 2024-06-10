% helper function for nice printouts
function[significand,exponent] = convert2base(num,base)

    exponent    = floor(log2(num) / log2(base));
    significand = base^(log2(num) / log2(base) - exponent);

end