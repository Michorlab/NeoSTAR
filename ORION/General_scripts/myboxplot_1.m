function h1=myboxplot_1(data1,group1,msize,mcolor,lwidth,csize)
%% my boxplot for mean+/-std
%  Jerry Lin 2019/08/23
%  Usage: myboxplot(data1,'b'group1,msize,mcolor,csize,lwidth)
%

%% Initialization
mytable = table(data1,group1);
stat1 = grpstats(mytable,'group1',{'mean','std'});

if (nargin < 3) 
    msize = 20;
    mcolor = 'r';
    lwidth = 4;
    csize = 20;
elseif (nargin <4)
    mcolor = 'r';
    lwidth = 4;
    csize = 20;
elseif (nargin <5)
    lwidth = 4;
    csize = 20;
else
    csize = 20;
end

%% plot 
h1=errorbar(1:size(stat1,1),stat1{:,3},stat1{:,4},'s','MarkerSize',msize,'MarkerEdgeColor',mcolor,'MarkerFaceColor',mcolor,'CapSize',csize,'LineWidth',lwidth);
set(gca,'xtick',1:size(stat1,1));
set(gca,'xticklabels',stat1{:,1});
xlims = xlim;
xlim([xlims(1)-0.5,xlims(2)+0.5]);
return;
