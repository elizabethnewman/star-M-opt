function[predictResults,trackingResults,trackingResultsMatrixPerSector] = indexTrackingPrediction(X, M, matrixResultsPerSector, initDate, endDate)

[ATest,BTest,BTestMat] = indexTrackingSetupData(initDate,endDate,'raw');

trackingResults = cell(1,length(X));
trackingResults{1} = fronorm(ATest(:,:) * X{1} - BTestMat);
for i = 2:length(X)
    trackingResults{i} = fronorm(mprod(ATest,X{i},M{i}) - BTest);
end

trackingResultsMatrixPerSector = cell(size(matrixResultsPerSector));
for i = 1:length(trackingResultsMatrixPerSector)
    trackingResultsMatrixPerSector{i} = fronorm(ATest(:,:,i) * matrixResultsPerSector.X{i} - BTest(:,1,i));
end



predictResults = cell(1,length(X));
for i = 1:length(X)
    predictResults{i} = indexTrackingPredictReturn(ATest, X{i});
end

end

