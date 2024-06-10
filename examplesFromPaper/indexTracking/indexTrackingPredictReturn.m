function[results] = indexTrackingPredictReturn(A,X)

priceChangeFromPurchase = indexTrackingComputePercentChange(A(:,:),0);
percentReturnsFromPurchase = priceChangeFromPurchase * X(:);

results.percentReturnsFromPurchase         = percentReturnsFromPurchase;
results.averagePercentReturnsFromPurchase  = mean(percentReturnsFromPurchase);
results.stdPercentReturnsFromPurchase      = std(percentReturnsFromPurchase);


priceChangePerDay    = indexTrackingComputePercentChange(A(:,:),1);
percentReturnsPerDay = priceChangePerDay * X(:);

results.percentReturnsPerDay        = percentReturnsPerDay;
results.averagePercentReturnsPerDay = mean(percentReturnsPerDay);
results.stdPercentReturnsPerDay     = std(percentReturnsPerDay);

end