function digitsSaveApproximationPlot(results,cmap)

if ~exist('cmap','var') || isempty(cmap), cmap = 'parula'; end

fig = figure(1); clf;
hCB = colorbar('east');
colormap(cmap)
hCB.TickLabels = [];
hCB.Ticks = [];
set(gca,'Visible',false)
hCB.Position = [0.425 0.05 0.14 0.9];
% hCB.Position(4) = 0.1000;

exportgraphics(fig,['digitsResults/',cmap,'.png'],'BackgroundColor','none')


[A,labels,nrmAPerClass,idx] = digitsSetupData(results.options);

nClasses = length(unique(labels));
for i = 0:nClasses - 1

    idxLabel    = find(labels == i);
    idxPerM = [];
    for Mname = {'MOpt','I','Z^T'}
        ii = strcmp(results.Mnames,Mname{1});
        [~,idx]  = sort(results.approxResults{ii}.res(idxLabel),'descend');  
        
        if strcmp(Mname{1},'MOpt')
            idxPerM = cat(2,idxPerM,[idxLabel(idx(1)),idxLabel(idx(end))]);
        else
            idxPerM = cat(2,idxPerM,idxLabel(idx(1)));
        end
    end
    
    dirName     = ['digitsResults/',results.options.filename,'/class_',num2str(i),'/'];
    subDirNames = {'MOptBest/','MOptWorst/','MIBest/','MZBest/'};
    for j = 1:length(idxPerM)
        
        ii = idxPerM(j);
    
    
        figure(1); clf;
        imagesc(A(:,:,ii));
        axis('square')
        axis('image')
        axis('off')
        myCLim = clim;
            
        tmpName = [dirName,subDirNames{j}];
        if ~exist(tmpName, 'dir'), mkdir(tmpName); end
        
        exportgraphics(fig,[tmpName,'orig.png'],'BackgroundColor','none')
  
   
        % get clim for difference
        myCLimErr = [];
        for jj = 1:length(results.Mnames)
            tmp = abs(A(:,:,ii) - results.Ak{jj}(:,:,ii));
            myCLimErr = cat(2,myCLimErr,[min(tmp(:)); max(tmp(:))]);
        end
        myCLimErr = [min(myCLimErr(1,:)); max(myCLimErr(2,:))];
        
        
        save([tmpName,'colorbar'],'myCLim','myCLimErr')

        count = 1;
        for jj = 1:length(results.Mnames)
            figure(1); clf;
            imagesc(results.Ak{jj}(:,:,ii));
            clim(myCLim);
            axis('square')
            axis('image')
            axis('off')
            exportgraphics(fig,[tmpName,results.Mnames{jj},'_approx.png'],'BackgroundColor','none')
    
            imagesc(abs(A(:,:,ii) - results.Ak{jj}(:,:,ii)));
            axis('square')
            axis('image')
            axis('off')
            clim(myCLimErr);
            exportgraphics(fig,[tmpName,results.Mnames{jj},'_diff.png'],'BackgroundColor','none')
            
            count = count + 1;

        end
    
    end

end


end

