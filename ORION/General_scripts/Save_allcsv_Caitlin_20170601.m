slideName = {'CMD1','CMD2','CMD3','CMD5','CMD6','CMD7','CMD9','CMD10','CMD11','CMD12'};
titles = {'Ctrl1','Ctrl2','PalbL1','PalbL2','PalbH1','PalbH2','AbemL1','AbemL2','AbemH1','AbemH2'};



for slide = 1:length(slideName)
    myName = slideName{slide};
    filename = strcat(titles{slide},'.csv');
    eval(strcat('temp1=sample',myName,';'));
    writetable(temp1,filename);
end