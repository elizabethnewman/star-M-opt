function romPlotErrorPerParam(results,plotOptions)

if ~exist('plotOptions','var') || isempty(plotOptions)
    plotOptions = {'LineWidth',1,'MarkerSize',15};
end

names        = results.Mnames;
cList        = linspace(results.options.cMin, results.options.cMax, results.options.cN);
nrmAPerParam = results.dataInfo.nrmAPerParam;

for i = 1:length(names)
    errPerParam = results.approxResults{i}.errPerParam;
    if ~strcmp(names{i},'matrix')
        semilogy(cList(:), errPerParam(:) ./ nrmAPerParam(:), '-o', 'DisplayName', names{i}, plotOptions{:})
        hold on;
    end
    
end
legend()
xlabel('c')
ylabel('error')
hold off;
set(gca,'FontSize',18)

end

