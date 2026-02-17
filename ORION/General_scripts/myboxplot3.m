function h1=myboxplot3(data1,group1,fontsize,sw1,grouporder)
%% My function boxplot3
%  Jerry Lin 2021/02/03
%  Usage myboxplot3(groupdata,groups,grouporder)
% 
%  2025/01/12   add support for number matrix only
%

if nargin <=1
   h1 = boxplot(data1);
   sw1 = false;
   fontsize = 9;
elseif nargin ==2
   h1=boxplot(data1,group1);
   fontsize=9;
   sw1 = false;
elseif nargin == 3
   h1=boxplot(data1,group1);
   sw1 = false; 
elseif nargin == 4
   h1=boxplot(data1,group1);
elseif nargin == 5
    h1=boxplot(data1,group1,'GroupOrder',grouporder);
end


if nargin > 1
    group3 = grp2idx(group1);
elseif nargin >=5
    group3 = zeros(length(group1));
    group3(ismember(group1,grouporder{1}))=1;
    group3(ismember(group1,grouporder{2}))=2;
    group3(ismember(group1,grouporder{3}))=3;
end

if sw1
    str1 = 'p=';
else
    str1 = '';
end

ylims = ylim;
step1 = 0.1*(ylims(2)-ylims(1));
ylim([ylims(1),ylims(2)+4*step1]);
line([1.1,1.9],[ylims(2),ylims(2)],'LineWidth',0.5,'Color','k');
line([2.1,2.9],[ylims(2),ylims(2)],'LineWidth',0.5,'Color','k');
line([1.1,2.9],[ylims(2)+2*step1,ylims(2)+2*step1],'LineWidth',0.5,'Color','k');

if nargin ==1
    d1 = data1(:,1);
    d2 = data1(:,2);
    d3 = data1(:,3);
else    
    d1 = data1(group3==1);
    d2 = data1(group3==2);
    d3 = data1(group3==3);
end

[~,p12]=ttest2(d1,d2,'Vartype','unequal');
[~,p23]=ttest2(d2,d3,'Vartype','unequal');
[~,p13]=ttest2(d1,d3,'Vartype','unequal');

if p12 < 0.05
    clr1 = 'r';
else
    clr1 = 'k';
end

if p12>0.001
    text(1.25,ylims(2)+step1,strcat(str1,num2str(p12,'%0.3f')),'FontSize',fontsize,'color',clr1);
else
    text(1.25,ylims(2)+step1,'p<0.001','FontSize',fontsize,'color',clr1);
end

if p23 < 0.05
    clr2 = 'r';
else
    clr2 = 'k';
end

if p23>0.001
    text(2.25,ylims(2)+step1,strcat(str1,num2str(p23,'%0.3f')),'FontSize',fontsize,'color',clr2);
else
    text(2.25,ylims(2)+step1,'p<0.001','FontSize',fontsize,'color',clr2);
end

if p13 < 0.05
    clr3 = 'r';
else
    clr3 = 'k';
end
if p13>0.001
    text(1.75,ylims(2)+3*step1,strcat(str1,num2str(p13,'%0.3f')),'FontSize',fontsize,'color',clr3);
else
    text(1.75,ylims(2)+3*step1,'p<0.001','FontSize',fontsize,'color',clr3);
end

line([1.1,1.9],[ylims(2),ylims(2)],'LineWidth',0.5,'Color',clr1);
line([2.1,2.9],[ylims(2),ylims(2)],'LineWidth',0.5,'Color',clr2);
line([1.1,2.9],[ylims(2)+2*step1,ylims(2)+2*step1],'LineWidth',0.5,'Color',clr3);

return;
