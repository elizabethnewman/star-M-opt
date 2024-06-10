function romSaveApproximationPlot(results,idx)

model = romSetupWaveEquation();
A     = results.A;
Ak    = results.Ak;
names = results.Mnames;

myClim = [min(A(:,:,idx),[],1); max(A(:,:,idx),[],1)];

myClimErr = [];
for i = 1:length(names)
    B = abs(A(:,:,idx) - Ak{i}(:,:,idx));
    tmp = [min(B,[],1);max(B,[],1)];
    myClimErr = cat(3,myClimErr,[-max(abs(tmp),[],1);max(abs(tmp),[],1)]);
end

myClimErr = [min(myClimErr(1,:,:),[],3); max(myClimErr(2,:,:),[],3)];

myCmap = parula(256);

n    = 512;
tmp1 = hot(512);
tmp2 = tmp1(:,[2,3,1]);
myCmapErr = [tmp2(floor(n / 2):-1:1,:);tmp1(1:floor(n / 2),:)];



dirName = ['romResults/',date,'--idx_',num2str(idx)];
if ~exist(dirName, 'dir'), mkdir(dirName); end

save([dirName,'/colorbar'], 'myClim', 'myClimErr')

cbarNames = {'colorbar_approx','colorbar_diff'};
count = 1;
for cmap = {myCmap, myCmapErr}
    fig = figure(1); clf;
    hCB = colorbar('east');
    colormap(cmap{1})
    hCB.TickLabels = [];
    hCB.Ticks = [];
    set(gca,'Visible',false)
    hCB.Position = [0.425 0.05 0.14 0.9];
    % hCB.Position(4) = 0.1000;
    
    exportgraphics(fig,[dirName,'/',cbarNames{count},'.png'],'BackgroundColor','none')
    count = count + 1;
end


saveImg(model,A(:,:,idx),myCmap,myClim,[dirName,'/orig/'],'img_');

for i = 1:length(names)
    saveImg(model,Ak{i}(:,:,idx),myCmap,myClim,[dirName,'/',names{i},'/'],'img_');
    saveImg(model,Ak{i}(:,:,idx) - A(:,:,idx),myCmapErr,myClimErr,[dirName,'/',names{i},'/'],'diff_');
end

end


function saveImg(model,B,cmap,myClim,dirName,fileName)

if ~exist(dirName, 'dir'), mkdir(dirName); end

fig = figure(1);
for i = 1:size(B,2)
    clf;
    h = pdeplot(model,"XYData",B(:,i),Colormap=cmap);
    
    if myClim(1,i) < myClim(2,i), clim(myClim(:,i)); end
    h(2).Visible = 'off';

    axis('square')
    axis('off')
    exportgraphics(fig,[dirName,fileName,num2str(i),'.png'],'BackgroundColor','none')
end

end


