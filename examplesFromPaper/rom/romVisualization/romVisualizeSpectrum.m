
function[] = romVisualizeSpectrum(results)

A = results.A;

% colorbar
mymap = jet(size(A,3));

fig = figure(1); clf;
hCB = colorbar('south');
colormap(mymap);
hCB.TickLabels  = [];
hCB.Ticks       = [];
set(gca,'Visible',false)
hCB.Position    = [0.05 0.425 0.9 0.05];
% hCB.Position(4) = 0.1000;
% exportgraphics(fig,'spectra_initZ^T/colorbar.png','BackgroundColor','none')

for i = [3,4,6,7]
    figure(i); clf;
    SkHat = modeProduct(results.Sk{i},results.M{i});
    AHat  = modeProduct(results.A,results.M{i});
    
    % diag(mean(SkHat,3))
    for k = 1:size(A,3)

        sk = diag(SkHat(1:5,1:5,k));
        csk = cumsum(sk.^2) / fronorm(AHat(:,:,k))^2;

        % semilogy(diag(diag(SkHat(1:5,1:5,k))).^2 / results.dataInfo.nrmAPerParam(k)^2,'o-','LineWidth',3,'MarkerSize',10,'DisplayName',num2str(k),'Color',mycolors(k,:))
        plot(csk,'o-','LineWidth',3,'MarkerSize',10,'DisplayName',num2str(k),'Color',mymap(k,:))

        hold on;
    end
    hold off;
    colormap(mymap);
    % cbh = colorbar ; 
    % 
    % cbh.Ticks = linspace(0, 1, 11) ; 
    % cbh.TickLabels = num2cell([1,5:5:50]);

    xlabel('number of singular values stored')
    ylabel('percentage of energy')
    ylim([0,1])
    grid on;
    set(gca,'FontSize',18);
    
    % matlab2tikz(['spectra_initZ^T/spectrum_',results.Mnames{i},'.tex'])
end

figure(2); clf;
for i = [3,4,6,7]
    AHat = modeProduct(results.A,results.M{i});
    SkHat = modeProduct(results.Sk{i},results.M{i});
    tmp  = (reshape(sum(AHat.^2,[1,2]),[],1)) / fronorm(AHat)^2;
    % tmp  = (reshape(sum(SkHat(1:2,1:2,:).^2,[1,2]),[],1)) / fronorm(AHat)^2;
    % sum(tmp(:))
    plot(tmp, '-o','MarkerSize',10,'LineWidth',3,'DisplayName',results.Mnames{i})
    hold on;
end
hold off;
legend()
xlabel('frontal slice')
ylabel('||hat{X}(:,:,i)||_F^2 / ||X||_F^2')
set(gca,'FontSize',18)
% matlab2tikz('spectra_initZ^T/frontal_slices.tex')


end

%% 
% load romResults/05-Sep-2023--k-02_initC.mat


%% 
% look at the spectrum of the norms of the singular tubes
% figure(1); clf;
% for i = 3:length(results.M)
%     semilogy(diag(vecnorm(results.Sk{i},2,3)),'-o','LineWidth',3,'MarkerSize',10,'DisplayName',results.Mnames{i})
%     hold on;
% end
% hold off;
% legend()

% %% 
% load romResults/05-Sep-2023--k-02_initZ.mat
% 
% 
% mycolors = jet(50);
% for i = [3,4,6,7]
%     fig = figure(i); clf; fig.Name = results.Mnames{i};
%     SkHat = modeProduct(results.Sk{i},results.M{i});
%     AHat  = modeProduct(results.A,results.M{i});
% 
%     diag(mean(SkHat,3))
%     for k = 1:50
%         if mod(k,3) == 0, mm = '^'; end
%         if mod(k,3) == 1, mm = 'o'; end
%         if mod(k,3) == 2, mm = 'x'; end
% 
%         sk = diag(SkHat(1:5,1:5,k));
%         csk = cumsum(sk.^2) / fronorm(AHat(:,:,k))^2;
% 
%         % semilogy(diag(diag(SkHat(1:5,1:5,k))).^2 / results.dataInfo.nrmAPerParam(k)^2,'o-','LineWidth',3,'MarkerSize',10,'DisplayName',num2str(k),'Color',mycolors(k,:))
%         plot(csk,'o-','LineWidth',3,'MarkerSize',10,'DisplayName',num2str(k),'Color',mycolors(k,:))
% 
%         hold on;
%     end
%     hold off;
%     colormap(mycolors)
%     cbh = colorbar ; 
% 
%     cbh.Ticks = linspace(0, 1, 11) ; 
%     cbh.TickLabels = num2cell([1,5:5:50]);
% 
%     xlabel('singular value index')
%     ylabel('singular value')
%     % ylim([1e-8,1e1])
%     ylim([0.45,1])
%     grid on;
%     set(gca,'FontSize',18)
%     % matlab2tikz(['~/Desktop/spectrum_',results.Mnames{i},'.tex'])
% end
% 
% 
% % tmp = cellfun(@num2str,mat2cell(1:50,1,ones(1,50)),'UniformOutput',false);
% % colorbar('Ticks',1:50,'TickLabels',[tmp{:}])
% % legend()
% 
% %% 
% 
% figure(10); clf;
% for i = [6,7]
%     AHat  = modeProduct(results.A,results.M{i});
%     % disp(fronorm(AHat).^2)
%     tmp = reshape(sum(AHat.^2,[1,2]),[],1) / fronorm(AHat).^2;
%     semilogy(tmp, '-o','MarkerSize',20,'LineWidth',5,'DisplayName',results.Mnames{i})
%     hold on;
% end
% hold off;
% legend()


%% 


