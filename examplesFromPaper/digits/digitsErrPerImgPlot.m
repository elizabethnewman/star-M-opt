fig = figure(1); clf;
for k = 5
    load(sprintf('digitsResults/07-Jun-2023--k-%0.2d.mat',k));
    
    [~,idx] = sort(results.approxResults{end}.res);
    for i = [1,3, 5, 7]
        % plot(results.approxResults{i}.res(idx),'-','LineWidth',3,'DisplayName',results.Mnames{i})
        histogram(results.approxResults{i}.res(idx),'DisplayName',results.Mnames{i},'Normalization','count')

        hold on;
        disp([results.Mnames{i}, ': ', num2str(results.approxResults{i}.err)])
    end
    
end
legend()
xlabel('error')
ylabel('number of images')
grid on
hold off;
set(gca,'FontSize',18)

% matlab2tikz('digitsResults/07-Jun-2023--errPerImgHistogram.tex')