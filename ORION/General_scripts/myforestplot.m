function h1=myforestplot(params,odds_ratios,conf_intervals)

%% Parameters and their corresponding Odds Ratios and Confidence Intervals
%params = {'Distance travelled for use', 'Year', 'Wealth index', 'Distance to deep tube well', ...
%          'Deep tube well use', 'Maternal education', 'Population aged <5 y'};
%odds_ratios = [1.002, 0.817, 0.963, 1.000, 0.513, 0.965, 1.095];
%conf_intervals = [0.998 1.006; 0.753 0.886; 0.876 1.059; 0.999 1.000; 0.365 0.722; 0.931 1.000; 1.076 1.111];

% Create the figure
h1=figure('units','normalized','outerposition',[0.5 0 0.5 1]);
hold on;

% Number of parameters
n_params = length(params);

% Plot each parameter

for i = 1:n_params
    if conf_intervals(i,1)>1 || conf_intervals(i,2)<1
        mycolor = 'r';
    else
        mycolor = 'k';
    end
    plot(conf_intervals(i,:), [i i], '-', 'Color', mycolor);
    plot(odds_ratios(i), i, 'ko', 'MarkerFaceColor', mycolor,'MarkerSize',8);
end

set(gca, 'ydir', 'reverse', 'ytick', 1:n_params, 'yticklabel', params,'TickLabelInterpreter','none');
xlabel('Odds Ratio');
plot([1 1], [0 n_params+1], 'k--');     % Add a reference line at OR = 1

%xlim([0.2 1.2]);

title('Forest Plot of Odds Ratios with 95% Confidence Intervals');

grid on;

hold off;
set(gcf,'color','w');
return;