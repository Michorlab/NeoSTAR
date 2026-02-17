function [TPR,FPR]= myTPR_FPR(listT,listP)
%% function to calculate TPR and FPR for ROC curve
%  input listT (true label) & listP (predicted label)
%  return TPR = TP/(TP+FN) & FPR = FP/(FP+TN)



%% Confusion matrix elements
TP = sum((listT == 1) & (listP == 1));
FP = sum((listT == 0) & (listP == 1));
TN = sum((listT == 0) & (listP == 0));
FN = sum((listT == 1) & (listP == 0));

%% TPR and FPR
TPR = TP / (TP + FN);  % Sensitivity or Recall
FPR = FP / (FP + TN);  % 1 - Specificity

return;
