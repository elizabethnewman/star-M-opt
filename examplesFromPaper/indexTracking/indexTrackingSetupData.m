function[A,B,BMat,sectors,SP500] = indexTrackingSetupData(initDate,endDate,percentChangeFlag, saveFlag, fileName, dirName)

if ~exist('percentChangeFlag','var'), percentChangeFlag = 'raw'; end
if ~exist('saveFlag','var') || isempty(saveFlag), saveFlag = false; end
if ~exist('fileName','var') || isempty(fileName), fileName = 'tmp'; end
if ~exist('dirName','var') || isempty(dirName),   dirName  = 'indexTrackingResults/'; end



sectors.communicationServices   = {'CSCO', 'TMUS', 'VZ','CMCSA','AMX','ORAN','DIS','T','DASH','ZM'};
sectors.consumerDiscretionary   = {'TGT', 'AMZN', 'WMT','RCL','HD','LVMHF','TM','MCD', 'NKE','SBUX'};
sectors.consumerStaples         = {'PG','KO', 'PEP','NSRGY','LRLCY','COST','PM','UL','BUD','EL'};
sectors.energy                  = {'XOM', 'CVX', 'COP','SHEL','TTE','SLB','BP','EQNR','PBR','EOG'};
sectors.financials              = {'JPM', 'BAC', 'WFC','HSBC','HDB','MS','SCHW','TD','GS','C'};
sectors.healthcare              = {'UNH', 'JNJ', 'LLY','NVO','MRK','RHHBY','PFE','TMO','ABT','DHR'};
sectors.industrials             = {'UPS', 'BA', 'ACN','CAT','RTX','HON','UNP','GE','DE','ADP'};
sectors.informationTechnology   = {'AAPL', 'MSFT', 'GOOGL','IBM','CRM','CSCO','FSLR','ACN','ENPH','AVGO'};
sectors.materials               = {'LIN', 'FCX', 'SHW','BHP','RIO','APD','SCCO','ECL','GOLD','VALE'};
sectors.realEstate              = {'SPG', 'PSA', 'O','PLD','AMT','EQIX','CCI','WELL','CSGP','DLR'};
sectors.utilities               = {'NEE', 'D', 'ED','SO','DUK','NGG','SRE','D','AEP','XEL'};


sectors.names                   = fields(sectors);


% https://www.spglobal.com/spdji/en/indices/equity/sp-500-health-care-sector/#overview

SP500.communicationServices   = {'^SP500-50'};
SP500.consumerDiscretionary   = {'^SP500-25'};
SP500.consumerStaples         = {'^SP500-30'};
SP500.energy                  = {'^GSPE'};
SP500.financials              = {'^SP500-40'};
SP500.healthcare              = {'^SP500-35'};
SP500.industrials             = {'^SP500-20'};
SP500.informationTechnology   = {'^SP500-45'};
SP500.materials               = {'^SP500-15'};
SP500.realEstate              = {'^SP500-60'};
SP500.utilities               = {'^SP500-55'};
SP500.names                   = fields(SP500);
SP500.index                   = {'^GSPC'};

%% gather data

% CHANGES (Oct. 3, 2024): new data access fields

A = [];
d = [];
Aj = [];
for i = 1:length(sectors.names)
    mySector = sectors.(sectors.names{i});
    for j = 1:length(mySector)
        yahoo_raw = getMarketDataViaYahoo(mySector{j}, initDate, endDate);
        
        if isempty(d)
            % d = yahoo_raw.Date; 
            d = yahoo_raw.chart.result.timestamp;
            idx = 1:numel(d);
        else
            [d,idx] = updateDate(d,yahoo_raw);
            Aj = Aj(idx,:);
        end
        % tmp = yahoo_raw.AdjClose;
        tmp = yahoo_raw.chart.result.indicators.adjclose.adjclose;
        if numel(tmp) > numel(idx), tmp = tmp(idx); end
        Aj        = cat(2,Aj,tmp);
    end
    A = cat(3,A,Aj);

    % reset Aj
    Aj = zeros(size(Aj,1),0);
end

% right-hand side (for all data)
yahoo_raw   = getMarketDataViaYahoo(SP500.index{1}, initDate, endDate);
[d,idx] = updateDate(d,yahoo_raw);
A       = A(idx,:,:);
% BMat    = yahoo_raw.AdjClose(idx);
BMat = yahoo_raw.chart.result.indicators.adjclose.adjclose(idx);

% right-hand side (per sector
B = zeros(size(A,1),1,0);
for i = 1:length(SP500.names)
    mySectorIndex   = SP500.(SP500.names{i});
    yahoo_raw       = getMarketDataViaYahoo(mySectorIndex{1}, initDate, endDate);
   
    % update date info
    [d,idx] = updateDate(d,yahoo_raw);
    
    % tmp = yahoo_raw.AdjClose;
    tmp = yahoo_raw.chart.result.indicators.adjclose.adjclose;
    if numel(tmp) > numel(idx), tmp = tmp(idx); end

    B = cat(3,B(idx,:,:),tmp);

    % update data
    A = A(idx,:,:);
    BMat = BMat(idx,:);
end


% saving price data (not percent change data) as separate variables 
% A0 = A;
% B0 = B;

% % compute percent change
switch percentChangeFlag
    case 'daily'
        dailyFlag = 1;
        A    = indexTrackingComputePercentChange(A,dailyFlag);
        B    = indexTrackingComputePercentChange(B,dailyFlag);
        BMat = indexTrackingComputePercentChange(BMat,dailyFlag);
    case 'first'
        dailyFlag = 0;
        A    = indexTrackingComputePercentChange(A,dailyFlag);
        B    = indexTrackingComputePercentChange(B,dailyFlag);
        BMat = indexTrackingComputePercentChange(BMat,dailyFlag);
end

if saveFlag
    save(sprintf([dirName, fileName,'--%s--%s'],initDate,endDate),"A", "B", "BMat")
end

end


function[d,idx] = updateDate(d,yahoo_raw)

% dNew = intersect(d,yahoo_raw.Date);
dNew = intersect(d, yahoo_raw.chart.result.timestamp);
idx  = datefind(dNew,d);
d    = dNew;

end

