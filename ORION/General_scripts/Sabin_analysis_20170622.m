

slideName = {'468PRI1','468193RP3','SCP01BT1','SCP3BT1','SCP17BT1','SCP28T1','SCP29T1','SCP32T1'};

figure;
eval(strcat('temp1=sample',slideName{1},';'));
temp2 = NaN(length(temp1.X),length(slideName));

for slide = 1:length(slideName)
    myName = slideName{slide};
    eval(strcat('temp1=sample',myName,';'));
    data1 = asinh(temp1.Ki67);
    subplot(2,4,slide);
    [idx,c,gate, pluscells] = Kmeangate(data1,2,1);
    title(myName);
    temp2(:,slide) = data1;
end
    
    
figure,boxplot(temp2);
