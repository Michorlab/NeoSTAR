function h=mydotplot(arrayP,arrayR,labelX,labelY)

r = size(arrayP,1);
c = size(arrayP,2);

arrayRow = repmat(1:r,c,1);
arrayCol = repmat((1:c)',1,r);
scalefactor = 5;

h=scatter(arrayRow(:),arrayCol(:),power(-log(arrayP(:)),1.6)*scalefactor+1,arrayR(:),'filled');
colormap(redbluecmap);
colorbar;
set(gca,'xtick',1:r);
set(gca,'xticklabels',labelX,'TickLabelInterpreter','none');
set(gca,'ytick',1:c);
set(gca,'yticklabels',labelY,'TickLabelInterpreter','none');
caxis([0 2]);

hold on;
scatter(r+1,1,power(-log(0.05),1.6)*scalefactor+1,'k');
text(r+0.5,2,'p=0.05');
scatter(r+1,3,power(-log(0.01),1.6)*scalefactor+1,'k');
text(r+0.5,4,'p=0.01');
scatter(r+1,5,power(-log(0.001),1.6)*scalefactor+1,'k');
text(r+0.5,6,'p=0.001');
xlims = xlim;
xlim([xlims(1),xlims(2)+1]);

return
