function [p,mean_r]=myttest2(data1,group1)
%% My ttest2 for array with group
%  Jerry Lin 2024/03/29
%

if nargin > 1
    group2 = grp2idx(group1);
    d1 = data1(group2==1);
    d2 = data1(group2==2);
else
    d1 = data1(:,1);
    d2 = data1(:,2);
end

[~,p]=ttest2(d1,d2,'Vartype','unequal');
mean_r = mean(d1)./mean(d2);

return;
