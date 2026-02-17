function [h1,p]=myboxplot2p(data1,groupcolor)
%% My function boxplot2
%  Jerry Lin 2025/02/04
%  Usage myboxplot2(groupdata,groups,grouporder)

h1 = boxplot(data1);
hold on;

ylims = ylim;
step1 = 0.1*(ylims(2)-ylims(1));
ylim([ylims(1),ylims(2)+2*step1]);
line([1,2],[ylims(2),ylims(2)],'LineWidth',1);

fontsize = 12;

d1 = data1(:,1);
d2 = data1(:,2);

colormatrix = zeros(length(groupcolor),3);
for i = 1:length(groupcolor)
    switch(groupcolor{i})
        case 'r'
            colormatrix(i,:)=[1,0,0];
        case 'g'
            colormatrix(i,:)=[0,1,0];
        case 'b'
            colormatrix(i,:)=[0,0,1];
    end
end


scatter(repmat(1,length(d1),1),d1,20,colormatrix,'filled');
scatter(repmat(2,length(d2),1),d2,20,colormatrix,'filled');
for i =1:length(d1)
    line([1 2],[d1(i),d2(i)],'Color',groupcolor{i},'LineWidth',0.5);
end
boxplot(data1);
hold off;
[~,p]=ttest2(d1,d2,'Vartype','unequal');

if p<0.05
    color1 ='r';
    fontsize = fontsize+1;
else
    color1 = 'k';
end
if p>0.001
    text(1.3,ylims(2)+step1,strcat('p=',num2str(p,'%0.3f')),'FontSize',fontsize,'Color',color1);
else
    text(1.3,ylims(2)+step1,'p<0.001','FontSize',fontsize,'Color',color1);
end
return;
