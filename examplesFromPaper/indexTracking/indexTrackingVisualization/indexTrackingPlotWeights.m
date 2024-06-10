function[fig] = indexTrackingPlotWeights(results)

fig = figure;
n = length(results.sectors.names);

colors = jet(n);
x1 = 0.5;
x2 = 10.5;
for i = 1:n

    area([x1, x1, x2, x2], [0; 1; 1; 0], 'FaceColor', colors(i,:),'FaceAlpha', 0.75, 'HandleVisibility', 'off');
    hold on;
    x1 = x2;
    x2 = x2 + 10;
end

set(gca,'ColorOrderIndex',1)
lineColors = get(gca,'ColorOrder');
for i = 1:length(results.X)
    plot(results.X{i}(:), '-o', 'LineWidth', 3, 'MarkerSize', 15, 'MarkerEdgeColor','k', 'MarkerFaceColor',lineColors(i,:),'DisplayName', results.Mnames{i})
end

xlabel('stock')
xlim([0.5, 110.5])
ylabel('weight')
% ylim([0, 0.15])
ylim([0, 1])
legend()
set(gca, 'FontSize', 18)
set(gca, 'XTick', 2:10:110)
set(gca, 'XTickLabels', results.sectors.names)
hold off;

end

