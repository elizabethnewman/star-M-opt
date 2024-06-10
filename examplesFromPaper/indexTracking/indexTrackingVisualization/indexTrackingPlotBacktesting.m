function[fig] = indexTrackingPlotBacktesting(results, startDate, endDate)


[A,B,~,~] = indexTrackingSetupData(startDate,endDate,'first');


% choice of M colors
colors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% colors = winter(length(results.M));

% plot portfolio (with m-product)
markers = {'-o','-+','-*','-x','-s','-d','-v','->','-h','-^'};


err  = zeros(11,length(results.M) - 1);

for i = 2:length(results.M)
    tmp = mprod(A,results.X{i},results.M{i});
    err(:,i-1) = sqrt(sum((tmp - B).^2,1));
end

fig = figure(1); clf;

[~,idx] = sort(results.sectors.names);

nrmB = sqrt(squeeze(sum(B.^2,1)));
for i = 2:length(results.M)
    plot(1:11,err(idx,i-1) ./ nrmB(idx),markers{i-1},'MarkerSize',15,'Color',colors(i-1,:),'LineWidth',2, 'DisplayName',results.Mnames{i})
    hold on;
end

set(gca,'FontSize',18)
xlabel('sector')
xticks(1:11)
xticklabels(results.sectors.names)
ylabel('rel. err.')
legend('Location','eastoutside')
hold off; 


end
