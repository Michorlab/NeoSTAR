function foldx = odds(pcells1,pcells2)
%% Calculate odd ratio for two array
%  Jerry Lin 2024/03/28
%

q1 = mean(pcells1 .* pcells2);
q2 = mean(pcells1 .* ~pcells2);
%q3 = mean(~pcells1 .* ~pcells2);
q4 = mean(~pcells1 .* pcells2);

foldx = mean(q1/((q1+q2)*(q1+q4)));
return;
