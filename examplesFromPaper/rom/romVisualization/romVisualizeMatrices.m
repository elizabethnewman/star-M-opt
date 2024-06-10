
function romVisualizeMatrices(results,dirName)

% create color range
cmin = -1; cmax = 1;
for i = 3:length(results.M)
    cmin = min(cmin, min(results.M{i}(:)));
    cmax = max(cmax, max(results.M{i}(:)));
end
[cmin, cmax]
% save('Mpics/colorlims',"cmax","cmin")

% choose colormap
mymap1 = hot(256);
mymap2 = mymap1(:,[3,2,1]);
mymap = [flipud(mymap2);mymap1];
mymap = mymap(21:end-20,:);


fig = figure(1); clf;
hCB = colorbar('east');
colormap(mymap)
hCB.TickLabels = [];
hCB.Ticks = [];
set(gca,'Visible',false)
hCB.Position = [0.425 0.05 0.14 0.9];
% hCB.Position(4) = 0.1000;

% exportgraphics(fig,'Mpics/colorbar.png','BackgroundColor','none')


for i = 3:length(results.M)
    fig = figure(i - 2); clf;
    imagesc(results.M{i});
    axis off; axis square;
    clim([cmin, cmax])
    colormap(mymap);
    exportgraphics(fig,[dirName,results.Mnames{i},'--init-',results.options.MInitialize,'.png'],'BackgroundColor','none')

    % colorbar; 
end

end
