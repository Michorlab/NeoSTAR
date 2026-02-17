function [h1,p]=myboxplot2dj(data1,group1,grouporder,fontsize)
%% My function boxplot2dj
%  Jerry Lin 2025/09/05
%  Usage myboxplot2(groupdata,groups,grouporder,fontsize)
%  add jitter and colors


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
line([1,2],[ylims(2),ylims(2)],'LineWidth',0.5,'Color','k');

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

x1 = 0.2*rand(length(d1),1)+0.9;
x2 = 0.2*rand(length(d2),1)+1.9;
scatter(x1,d1,25,'r','filled');
scatter(x2,d2,25,'b','filled');
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

% boxes = findobj(gca, 'Tag', 'Box');
% patchX = get(boxes(1), 'XData');
% patchY = get(boxes(1), 'YData');
% line(patchX, patchY, 'Color', 'r', 'LineWidth', 1);
% patchX = get(boxes(2), 'XData');
% patchY = get(boxes(2), 'YData');
% line(patchX, patchY, 'Color', 'b', 'LineWidth', 1);

% ---- change color -----

ax = gca;

boxes      = flipud(findobj(ax, 'Tag','Box'));
medians    = flipud(findobj(ax, 'Tag','Median'));
whiskU     = flipud(findobj(ax, 'Tag','Upper Whisker'));
whiskL     = flipud(findobj(ax, 'Tag','Lower Whisker'));
capsU      = flipud(findobj(ax, 'Tag','Upper Adjacent Value'));
capsL      = flipud(findobj(ax, 'Tag','Lower Adjacent Value'));
outliers   = flipud(findobj(ax, 'Tag','Outliers'));  % one line object per group

n = numel(boxes);
colors = [1,0,0;0,0,1];

for i = 1:n
    % Box outline
    set(boxes(i),   'Color', colors(i,:), 'LineWidth', 1);

    % Median line
    set(medians(i), 'Color', colors(i,:), 'LineWidth', 1);

    % Whiskers
    set(whiskU(i),  'Color', colors(i,:), 'LineWidth', 1, 'LineStyle','-');
    set(whiskL(i),  'Color', colors(i,:), 'LineWidth', 1, 'LineStyle','-');

    % Caps
    set(capsU(i),   'Color', colors(i,:), 'LineWidth', 1);
    set(capsL(i),   'Color', colors(i,:), 'LineWidth', 1);

    % Outliers (markers)
    if i <= numel(outliers) && isvalid(outliers(i))
        set(outliers(i), 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', 'none');
        % If you used 'Symbol','' above, there may be no outlier objects.
    end
end

box on


return;
