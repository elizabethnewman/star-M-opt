function[fig,err,nrmB] = indexTrackingPlotDataOneSectorTracking(results, startDate, endDate, diffFlag, idx, showFutureDate, A, B)

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

sectorNames = fields(results.sectors);
sectorNames = sectorNames(1:end-1);

% sector colors
cmap = jet(length(sectorNames));


% choice of M colors
colors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% colors = winter(length(results.M));

% plot indexTracking (with m-product)
markers = {'-o','-+','-*','-x','-s','-d','-v','->','-h','-^'};


fig = figure(1); clf;

err  = zeros(1,length(results.M) - 1);
ymin = min(B(:,1,idx),[],'all') - 0.05;
ymax = max(B(:,1,idx),[],'all') + 0.05;

b         = B(:,1,idx);
nrmB(idx) = fronorm(b);

if ~diffFlag
    plot(x(:), b,'Color',cmap(idx,:),'LineWidth',4, 'DisplayName',['S&P 500 - ', sectorNames{idx}])
end

hold on;

for i = 2:length(results.M)
    tmp = mprod(A,results.X{i},results.M{i});

    tmpmin = min(tmp(:,1,idx),[],'all');
    tmpmax = max(tmp(:,1,idx),[],'all');

    if tmpmin < ymin, ymin = tmpmin - 0.05; end
    if tmpmax > ymax, ymax = tmpmax + 0.05; end

    if diffFlag
        plot(x(1:size(tmp,1)), abs(tmp(:,1,idx) - b(1:size(tmp,1))),markers{i-1},'MarkerSize',15,'Color',colors(i-1,:),'LineWidth',2, 'DisplayName',results.Mnames{i})
    else
        h = plot(x(1:size(tmp,1)), tmp(:,1,idx),markers{i-1},'MarkerSize',15,'Color',colors(i-1,:),'LineWidth',2, 'DisplayName',results.Mnames{i});
    end
    
    err(idx,i-1) = fronorm(tmp(:,1,idx) - b(1:size(tmp,1)));
end
    
% ii = (3 * (idx - 1) + 1):(3 * idx);
% tmp = A(:,:,i) * results.matrixPerSector.X{i};
% tmpmin = min(tmp(:,1),[],'all');
% tmpmax = max(tmp(:,1),[],'all');
% 
% if tmpmin < ymin, ymin = tmpmin - 0.05; end
% if tmpmax > ymax, ymax = tmpmax + 0.05; end
% 
% 
% % plot indexTracking (with m-product)
% if diffFlag
%     plot(x(:), abs(tmp(:) - b),'Color','w','LineWidth',4, 'DisplayName','matrix')
% else
%     plot(x(:), tmp(:),'Color','w','LineWidth',4, 'DisplayName','matrix')
% end


hold off;


grid on;

xlabel('time')

ylabel('return from start')
ylim([ymin, ymax])
datetick('x', 'mmm dd, yyyy', 'keepticks')
% xlim([x(1), x(end) + 0.01 * (x(end) - x(1))])
xlim([x(1), x(end)])
legend()

set(gca,'Color',[0.3,0.3,0.3]);
set(gca,'GridColor',[0.95, 0.95, 0.95]);
set(gca,'GridLineWidth',3);
set(gca,'FontSize',18)

hold off;

% subplot(3,4,12)
% 
% for i = 2:length(results.M)
%     plot(err(:,i-1) ./ nrmB,markers{i-1},'MarkerSize',15,'Color',colors(i-1,:),'LineWidth',2, 'DisplayName',results.Mnames{i})
%     hold on;
% end
% 
% set(gca,'FontSize',18)
% xlabel('sector')
% ylabel('rel. err.')
% legend('Location','eastoutside')
% hold off;



% % matrix case
% ymin = min(BMat,[],'all') - 0.05;
% ymax = max(BMat,[],'all') + 0.05;
% 
% fig = figure(2); clf; fig.Name = 'matrix vs. SP500';
% plot(x(:), BMat(:,1),'Color','w','LineWidth',4, 'DisplayName','S&P 500')
% 
% hold on;
% 
% tmp = A(:,:) * results.X{1};
% tmpmin = min(tmp(:,1),[],'all');
% tmpmax = max(tmp(:,1),[],'all');
% 
% if tmpmin < ymin, ymin = tmpmin - 0.05; end
% if tmpmax > ymax, ymax = tmpmax + 0.05; end
% 
% 
% % plot indexTracking (with m-product)
% plot(x(:), tmp(:),'Color','b','LineWidth',4, 'DisplayName','matrix')
% 
% grid on;
% legend()
% 
% xlabel('time')
% 
% ylabel('return from start')
% ylim([ymin, ymax])
% datetick('x', 'mmm dd, yyyy', 'keepticks')
% % xlim([x(1), x(end) + 0.01 * (x(end) - x(1))])
% xlim([x(1), x(end)])
% 
% set(gca,'Color',[0.2,0.2,0.2]);
% set(gca,'GridColor',[0.95, 0.95, 0.95]);
% set(gca,'GridLineWidth',3);
% set(gca,'FontSize',18)
% 
% hold off;


% fig = figure(3); clf; fig.Name = 'matrix vs. SP500-sector';
% plot(x(:), B(:,1,idx),'Color',cmap(idx,:),'LineWidth',4, 'DisplayName','S&P 500')
% 
% hold on;
% 
% ii = (3 * (idx - 1) + 1):(3 * idx);
% tmp = A(:,ii) * results.X{1}(ii);
% tmpmin = min(tmp(:,1),[],'all');
% tmpmax = max(tmp(:,1),[],'all');
% 
% if tmpmin < ymin, ymin = tmpmin - 0.05; end
% if tmpmax > ymax, ymax = tmpmax + 0.05; end
% 
% 
% % plot indexTracking (with m-product)
% plot(x(:), tmp(:),'Color','b','LineWidth',4, 'DisplayName','matrix')
% 
% grid on;
% 
% xlabel('time')
% 
% ylabel('return from start')
% ylim([ymin, ymax])
% datetick('x', 'mmm dd, yyyy', 'keepticks')
% % xlim([x(1), x(end) + 0.01 * (x(end) - x(1))])
% xlim([x(1), x(end)])
% 
% set(gca,'Color',[0.2,0.2,0.2]);
% set(gca,'GridColor',[0.95, 0.95, 0.95]);
% set(gca,'GridLineWidth',3);
% set(gca,'FontSize',18)
% 
% hold off;




end

