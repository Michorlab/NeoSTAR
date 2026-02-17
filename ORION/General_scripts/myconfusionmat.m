function mat1 = myconfusionmat(label1,label2,flag1)
%% my own confusion matrix
%  Jerry Lin 2020/11/01


%% Initialization
if nargin < 3
    flag1 = true;
end

temp1 = tabulate(label1);
temp2 = tabulate(label2);

dim1 = size(temp1,1);
dim2 = size(temp2,1);

mat1 = zeros(dim1,dim2);

%% calculate confusion matrix
for i =1:dim1
    for j = 1:dim2
        mat1(i,j) = sum((ismember(label1,temp1(i,1))) & (ismember(label2,temp2(j,1))));
    end
end


%% Plots

if flag1
    figure('units','normalized','outerposition',[0.6 0 0.4 1]);

    subplot(3,1,1);
    imagesctext(mat1,12,'k','%2.0f');
    colorbar;
    set(gca,'xtick',1:dim2);
    set(gca,'xticklabel',temp2(:,1));
    set(gca,'ytick',1:dim1);
    set(gca,'yticklabel',temp1(:,1));
    xlabel('label 2');
    ylabel('label 1');
    title('All counts');
    ylims = ylim;
    for i = 1:dim2-1
        line([i+0.5 i+0.5],[ylims(1),ylims(2)],'color','k','LineStyle','--','LineWidth',1.5);
    end

    subplot(3,1,2);
    imagesctext(mat1./sum(mat1,1)*100,12,'k','%2.0f','%');
    colorbar;
    set(gca,'xtick',1:dim2);
    set(gca,'xticklabel',temp2(:,1));
    set(gca,'ytick',1:dim1);
    set(gca,'yticklabel',temp1(:,1));
    xlabel('label 2');
    ylabel('label 1');
    title('Normalized by column');
    ylims = ylim;
    for i = 1:dim2-1
        line([i+0.5 i+0.5],[ylims(1),ylims(2)],'color','k','LineStyle','--','LineWidth',1.5);
    end
    %caxis([0 0.6]);

    subplot(3,1,3);
    imagesctext(mat1./sum(mat1,2)*100,12,'k','%2.0f','%');
    colorbar;
    set(gca,'xtick',1:dim2);
    set(gca,'xticklabel',temp2(:,1));
    set(gca,'ytick',1:dim1);
    set(gca,'yticklabel',temp1(:,1));
    xlabel('label 2');
    ylabel('label 1');
    title('Normalized by row');
    %caxis([0 0.6]);
    xlims = xlim;
    for i = 1:dim1-1
        line([xlims(1),xlims(2)],[i+0.5 i+0.5],'color','k','LineStyle','--','LineWidth',1.5);
    end
    colormap(cool);
end
set(gcf,'color','w');
return;
