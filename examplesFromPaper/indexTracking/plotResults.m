
% reload data

[A,B,BMat] = indexTrackingSetupData(results.options.initDateTrain,...
    results.options.endDateTrain);

[ATest,BTest,BTestMat] = indexTrackingSetupData(results.options.initDateTest,...
    results.options.endDateTest);

%% 
AMat = A(:,:);
priceChangeFromPurchase = ((AMat(2:end,:) - A(1,:)) ./ A(1,:));
priceChangePerDay       = ((AMat(2:end,:) - A(1:end-1,:)) ./ A(1:end-1,:));

%% plot truth
figure(1); clf;

% plot(BMat(:), 'k', 'LineWidth', 1)
% hold on;
% semilogy(abs(A(:,:) * results.X{1} - BMat(:)))

AMat = priceChangeFromPurchase(:,:);
for i = 1:size(AMat,2)
    plot(AMat(:,i),'LineWidth', 2)
    hold on;
end
% legend()

% hold off;
% for k = 1:size(B,3)
%     for 
%     plot()
% end

hold off;

% names = {};
% 
% for i = 1:length(sectors.names)
%     names = cat(2,names,sectors.(sectors.names{i}));
% end

for i = 1:size(AMat,2)
    text(size(AMat,1) + 1, AMat(end,i), names{i})
end



%% plot results

figure(1); clf;
subplot(1, 3, 1)
for i = 1:length(results.predictResults)
    tmp = results.predictResults{i}.percentReturnsFromPurchase;
    plot(tmp, 'LineWidth', 2, 'DisplayName', results.Mnames{i})
    hold on;
end
% plot(BTestMat,'k', 'LineWidth', 2, 'DisplayName', results.Mnames{i})



hold off;
% legend()
grid on;
xlabel('time')
ylabel('return(M)')
set(gca,'FontSize',18)


% figure(2); clf;
subplot(1, 3, 2)
tmpOpt = results.predictResults{end}.percentReturnsFromPurchase;
for i = 1:(length(results.predictResults) - 1)
    tmp = results.predictResults{i}.percentReturnsFromPurchase;
    plot(tmp - tmpOpt, 'LineWidth', 2, 'DisplayName', results.Mnames{i})
    hold on;
    sum(tmp - tmpOpt < 0) / numel(tmp)
end

hold off;
%legend()
grid on;
xlabel('time')
ylabel('return(M) - return(M^*)')
set(gca,'FontSize',18)



% figure(3); clf;
subplot(1, 3, 3)
lineColors = get(gca,'ColorOrder');
for i = 1:length(results.predictResults)
    y = results.predictResults{i}.averagePercentReturnsFromPurchase;
    x = results.predictResults{i}.stdPercentReturnsFromPurchase;
    % plot(x, y, 'o', 'LineWidth', 2, 'DisplayName', results.Mnames{i})
    plot(x,y, 'o', 'LineWidth', 3, 'MarkerSize', 15, 'MarkerEdgeColor','k', 'MarkerFaceColor',lineColors(i,:),'DisplayName', results.Mnames{i})
    
    hold on;
    text(x, y, results.Mnames{i},'FontSize',18)
end
hold off;
legend();
grid on;
xlabel('standard deviation of returns')
ylabel('mean of returns')
set(gca,'FontSize',18)

title('percent from purchase')

%% plot results

figure(2); clf;
subplot(1, 3, 1)
for i = 1:length(results.predictResults)
    tmp = results.predictResults{i}.percentReturnsPerDay;
    plot(tmp, 'LineWidth', 2, 'DisplayName', results.Mnames{i})
    hold on;
end
% plot(BTestMat,'k', 'LineWidth', 2, 'DisplayName', results.Mnames{i})


hold off;
% legend()
grid on;
xlabel('time')
ylabel('return(M)')
set(gca,'FontSize',18)


% figure(2); clf;
subplot(1, 3, 2)
tmpOpt = results.predictResults{end}.percentReturnsPerDay;
for i = 1:(length(results.predictResults) - 1)
    tmp = results.predictResults{i}.percentReturnsPerDay;
    plot(tmp - tmpOpt, 'LineWidth', 2, 'DisplayName', results.Mnames{i})
    hold on;
    sum(tmp - tmpOpt < 0) / numel(tmp)
end

hold off;
%legend()
grid on;
xlabel('time')
ylabel('return(M) - return(M^*)')
set(gca,'FontSize',18)



% figure(3); clf;
subplot(1, 3, 3)
lineColors = get(gca,'ColorOrder');
for i = 1:length(results.predictResults)
    y = results.predictResults{i}.averagePercentReturnsPerDay;
    x = results.predictResults{i}.stdPercentReturnsPerDay;
    % plot(x, y, 'o', 'LineWidth', 2, 'DisplayName', results.Mnames{i})
    plot(x,y, 'o', 'LineWidth', 3, 'MarkerSize', 15, 'MarkerEdgeColor','k', 'MarkerFaceColor',lineColors(i,:),'DisplayName', results.Mnames{i})
    
    hold on;
    text(x, y, results.Mnames{i},'FontSize',18)
end
hold off;
legend();
grid on;
xlabel('standard deviation of returns')
ylabel('mean of returns')
set(gca,'FontSize',18)

title('percent per day')


%% 
clear sectors
sectors.informationTechnology   = {'AAPL', 'MSFT', 'GOOGL'};
sectors.healthcare              = {'UNH', 'JNJ', 'LLY'};
sectors.financials              = {'JPM', 'BAC', 'WFC'};
sectors.consumerDiscetionary    = {'TSLA', 'AMZN', 'WMT'};
sectors.communicationServices   = {'CSCO', 'TMUS', 'VZ'};
sectors.industrials             = {'V', 'MA', 'ACN'};
sectors.consumerStaples         = {'PG','KO', 'PEP'};
sectors.energy                  = {'XOM', 'CVX', 'COP'};
sectors.utilities               = {'NEE', 'D', 'ED'};
sectors.realEstate              = {'SPG', 'PSA', 'O'};
sectors.materials               = {'LIN', 'FCX', 'LYB'};
sectors.names                   = fields(sectors);

figure(5); clf;
for i = 1:length(results.Mnames)
    subplot(2, 3, i)
    pie(results.sectorDistribution{i})
    title(results.Mnames{i})
    colormap jet;
end
subplot(2,3,6)
pie((1 / 11) * ones(11,1), sectors.names)
colormap jet;
title('legend')

%% 

sectors.informationTechnology   = {'AAPL', 'MSFT', 'GOOGL'};
sectors.healthcare              = {'UNH', 'JNJ', 'LLY'};
sectors.financials              = {'JPM', 'BAC', 'WFC'};
sectors.consumerDiscetionary    = {'TSLA', 'AMZN', 'WMT'};
sectors.communicationServices   = {'CSCO', 'TMUS', 'VZ'};
sectors.industrials             = {'V', 'MA', 'ACN'};
sectors.consumerStaples         = {'PG','KO', 'PEP'};
sectors.energy                  = {'XOM', 'CVX', 'COP'};
sectors.utilities               = {'NEE', 'D', 'ED'};
sectors.realEstate              = {'SPG', 'PSA', 'O'};
sectors.materials               = {'LIN', 'FCX', 'LYB'};
sectors.names                   = fields(sectors);

figure(4); clf;

colors = parula(11);
x1 = 0.5;
x2 = 3.5;
for i = 1:size(A,3)

    area([x1, x1, x2, x2], [0; 1; 1; 0], 'FaceColor', colors(i,:),'FaceAlpha', 0.75, 'HandleVisibility', 'off');
    hold on;
    x1 = x2;
    x2 = x2 + 3;
end

set(gca,'ColorOrderIndex',1)
lineColors = get(gca,'ColorOrder');
for i = 1:length(results.X)
    plot(results.X{i}(:), '-o', 'LineWidth', 3, 'MarkerSize', 15, 'MarkerEdgeColor','k', 'MarkerFaceColor',lineColors(i,:),'DisplayName', results.Mnames{i})
end

xlabel('stock')
xlim([0.5, 33.5])
ylabel('weight')
ylim([0, 0.15])
legend()
set(gca, 'FontSize', 18)
set(gca, 'XTick', 2:3:33)
set(gca, 'XTickLabels', sectors.names)
hold off;



