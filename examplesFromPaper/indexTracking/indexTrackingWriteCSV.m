function[T] = indexTrackingWriteCSV(results, startDate, endDate, diffFlag, idx, showFutureDate, A, B)

if ~exist('diffFlag','var') || isempty(diffFlag), diffFlag = false; end

dataFlag = 0;
if (~exist('A','var') || isempty(A)) || (~exist('B','var') || isempty(B))
    dataFlag = 1;
    [A,B,~,~] = indexTrackingSetupData(startDate,endDate,'first');
end

flag = 0;
if exist('showFutureDate','var') && ~isempty(showFutureDate)
    flag = 1;
    
    if dataFlag
        % form new B as well
        [~,B,~,~] = indexTrackingSetupData(startDate,showFutureDate,'first');
    end
end

startDate = datetime(startDate, 'InputFormat', 'dd-MMM-yyyy');
endDate   = datetime(endDate, 'InputFormat', 'dd-MMM-yyyy');

if flag
    % x = linspace(startDate, showFutureDate, size(B,1));
    x = busdays(startDate, showFutureDate);
    x = x(1:size(B,1));
else
    % x = linspace(startDate, endDate, size(B,1));
    x = busdays(startDate, endDate);
    x = x(1:size(B,1));
end

b         = B(:,1,idx);

headers = cat(2,{'index'}, results.Mnames(2:end));
T = [b(:)];


for i = 2:length(results.M)
    tmp = mprod(A,results.X{i},results.M{i});
    n = numel(x) - size(tmp,1);
    z = NaN * ones(n,1);
    T = cat(2,T,[tmp(:,1,idx); z]);
end

T = array2table(T,'VariableNames',headers);
T = [array2table(x,'VariableNames',{'date'}),...
    array2table(datenum(x),'VariableNames',{'datenum'}), ...
    T];



end