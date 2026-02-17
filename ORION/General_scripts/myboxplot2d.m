function [h1,p]=myboxplot2d(data1,group1,grouporder,fontsize)
%% My function boxplot2
%  Jerry Lin 2020/04/17
%  Usage myboxplot2(groupdata,groups,grouporder)
%  2021/09/19 add features to array only plot

if nargin <=1
    h1 = boxplot(data1);
elseif nargin <=2
    h1=boxplot(data1,group1);
else
    h1=boxplot(data1,group1,'GroupOrder',grouporder);
end
hold on;

ylims = ylim;
step1 = 0.1*(ylims(2)-ylims(1));
ylim([ylims(1),ylims(2)+2*step1]);
line([1,2],[ylims(2),ylims(2)],'LineWidth',1);

if nargin <=3
    fontsize = 11;
end

if nargin > 1
    group2 = grp2idx(group1);
    d1 = data1(group2==1);
    d2 = data1(group2==2);
else
    d1 = data1(:,1);
    d2 = data1(:,2);
end
scatter(repmat(1,length(d1),1),d1,50,'k');
scatter(repmat(2,length(d2),1),d2,50,'k');
hold off;

[~,p3]=ttest2(d1,d2,'Vartype','unequal');
[~,p1]=ttest2(d1,d2,'Tail','left','Vartype','unequal');
[~,p2]=ttest2(d1,d2,'Tail','right','Vartype','unequal');
%p = min([p1,p2,p3]);
p = p3;

if p<0.05
    color1 ='r';
    fontsize = fontsize+2;
else
    color1 = 'k';
end

text(1.3,ylims(2)+step1,strcat('p=',num2str(p,'%0.3f')),'FontSize',fontsize,'Color',color1);

return;
