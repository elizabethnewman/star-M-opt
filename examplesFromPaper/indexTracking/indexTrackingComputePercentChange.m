function[A] = indexTrackingComputePercentChange(A,dailyFlag)


if ~exist('dailyFlag','var') || isempty(dailyFlag), dailyFlag = false; end

colons = repmat({':'}, ndims(A) - 1);

if dailyFlag
    % daily percent change
    A = (A(2:end,colons{:}) - A(1:end-1,colons{:})) ./ A(1:end-1,colons{:});
else
    % percent change from start date
    A = (A - A(1,colons{:})) ./ (A(1,colons{:}));
end

end

