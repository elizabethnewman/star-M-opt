function[fig] = lsPlotHyperplanes(A,B,X,M,hyperplaneFlag,spatialFlag)

if ~exist('M','var') || isempty(M), M = 1; end
if ~exist('hyperplaneFlag','var') || isempty(hyperplaneFlag), hyperplaneFlag = 1; end
if ~exist('spatialFlag','var') || isempty(spatialFlag), spatialFlag = 1; end

if ~spatialFlag && size(X,3) == 2
    % convert to transform domain
    A = modeProduct(A,M);
    B = modeProduct(B,M);
    X = modeProduct(X,M);
end

fig = figure(1); clf;
if size(X,3) == 2
    scatter3(A(:,1,1),A(:,2,1),B(:,1,1),100,'ko','MarkerFaceColor','r','DisplayName','$i=1$');
    hold on;
    scatter3(A(:,1,2),A(:,2,2),B(:,1,2),100,'ko','MarkerFaceColor','b','DisplayName','$i=2$');
else
    scatter3(A(:,1,1),A(:,2,1),B(:,1,1),100,'ko','MarkerFaceColor','k','DisplayName','$i=1,2$');
    hold on;
    scatter3(A(:,1,2),A(:,2,2),B(:,1,2),100,'ko','MarkerFaceColor','k','HandleVisibility','off');

end

if hyperplaneFlag
    if size(X,3) == 2
        fsurf(@(a1,a2) X(1,1,1) * a1 + X(2,1,1) * a2,[-3, 3],'FaceColor','r','FaceAlpha',0.5,'HandleVisibility','off');
        fsurf(@(a1,a2) X(1,1,2) * a1 + X(2,1,2) * a2,[-3, 3],'FaceColor','b','FaceAlpha',0.5,'HandleVisibility','off');
    else
        fsurf(@(a1,a2) X(1,1,1) * a1 + X(2,1,1) * a2,[-3, 3],'FaceColor','k','FaceAlpha',0.5,'HandleVisibility','off');
    end
end
hold off;


xlim([-3, 3])
ylim([-3, 3])
zlim([-8,8])
set(gca,'FontSize', 18)
set(gcf,'Color','w','Position', [500 100 600 600])
view(-45,45)

end


