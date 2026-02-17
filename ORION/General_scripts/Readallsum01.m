%% Read Summary data (for cell counting)
%  Jerry Lin 2016/10/07
% 
function [allsum,eachsum] = Readallsum01(myName,allcyc,myDIR)

%myDIR = uigetdir('F:\PAN0809');

%myName = '1188';
%allch = 6;

%allsum =zeros;
%totalframe = 143; %%rows*cols;
eachsum = cell(allcyc,1);

for ch = 1:allcyc;
    filename = strcat(myDIR,myName,'-cycle-',num2str(ch),'.csv');
    temp1 = import_summary01(filename,2,inf);
    allsum(:,ch)=temp1.Count;
    eachsum{ch} = temp1;
end

return;
