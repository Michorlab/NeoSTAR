%% Zoltan/Nicole  CycIF Analysis 2017/06/13


%% Gating data

slideName = {'BP10','CSS10','LP10','PS10','i074_10','i109_10','i151_10','i493_10','i779_10','i781_10','i943_10'};
%titles = {'D1400257A2','D1400257B7','D1400257B10','D1400617A4','D1400617B2','D1400617C4','D1401032A7','D1401032B1','D1401032B2','D1401032C1','98DTF14515','98DTF1697','98DTF1425','98DTF1433'};
titles = slideName;
markers = {'PD1','PD1_S','CD8A','CD4','MITF','STAT3','cJUN','Langerin','Langerin_S','CD3D','pSTAT1','PDL1_S','CD1c','CD1c_S','CD20','FoxP3','CD14','pSTAT3'};
gates = [500 400 1000 2000 400 2000 1000 3000 1000 1000 500 1000 4000 1000 2000 400 1500 500];
clear temp1 temp2 counts;

gatecounts = NaN(length(slideName),length(markers));

for slide = 1:length(slideName)
    myName = slideName{slide};
    eval(strcat('temp1=data',myName,';'));
    s1 = size(temp1,2);
    temp2 = temp1;
    %sample1 = datasample(temp1,5000);
    %eval(strcat('sample',myName,'=sample1;'));
 
    for m = 1:length(markers)
       data1 = temp1{:,markers(m)};
       mgate = gates(m);
       [pluscells, gate, peak,lowb,highb] = findgate3(log(data1+5),0,0.1,log(mgate));
       temp2{:,s1+m}=pluscells;
       temp2.Properties.VariableNames(s1+m) = strcat(markers(m),'p');
       gatecounts(slide,m) = mean(pluscells);
    end
    eval(strcat('gate',myName,'=temp2;'));
        
end

figure,imagesctext(gatecounts,16);
set(gca,'xtick',1:length(markers));
set(gca,'xticklabels',markers);
set(gca,'xticklabelrotation',45);
set(gca,'ytick',1:length(slideName));
set(gca,'yticklabels',slideName);
colormap(jet);