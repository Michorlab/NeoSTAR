function [fold_difference, lower_ci, upper_ci] = mygroupfold (data1, group1,sw1)
%% Calculate fold difference and CI from two gorup
%  jerry Lin 2025/03/04


%% initialization
[idx1,~,gl] = grp2idx(group1);
d1 = data1(idx1==1);
d2 = data1(idx1==2);
disp(gl);
mean1 = mean(d1); % Mean of group 1
std1 = std(d1);   % Standard deviation of group 1
n1 = length(d1);     % Sample size of group 1

mean2 = mean(d2);  % Mean of group 2
std2 = std(d2);    % Standard deviation of group 2
n2 = length(d2);     % Sample size of group 2

% Calculate the fold difference
fold_difference = mean1 / mean2;

% Log-transform the means and standard deviations
logMean1 = log(mean1);
logMean2 = log(mean2);
logStd1 = std1 / mean1;
logStd2 = std2 / mean2;

% Calculate the variance of the log-transformed ratio
variance = (logStd1^2 / n1) + (logStd2^2 / n2);

% Calculate the 95% confidence interval of the log-transformed ratio
ci95_log = 1.96 * sqrt(variance);

% Transform back to original scale
lower_ci = exp(log(fold_difference) - ci95_log);
upper_ci = exp(log(fold_difference) + ci95_log);

% Display the results
if sw1
    disp(['Fold Difference: ', num2str(fold_difference)]);
    disp(['95% Confidence Interval: [', num2str(lower_ci), ', ', num2str(upper_ci), ']']);
end

return;
