function[I] = facewiseEye(szA)
    I = eye(szA(1:2)) .* ones([1,1,szA(3:end)]);
end



