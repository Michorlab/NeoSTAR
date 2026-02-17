%% frame-by-frame neighbor analysis between EGFR & SMA cells for Sabin's Xenografts
%% 2017/06/22

myName = 'SCP29T1';
eval(strcat('data1=DATA',myName,';'));
%data1 = DATASCP01BT1;
EGFRcells = data1(data1.EGFR > exp(9),:);
SMAcells = data1(data1.SMA >exp(9),:);

maxframe = max(data1.frame);

cell2cell = cell(maxframe,5);
neighbor = NaN(maxframe,5);
cutoff = 15;


for f=1:maxframe
    EC = EGFRcells(EGFRcells.frame ==f,:);
    SC = SMAcells(SMAcells.frame ==f,:);
    AC = data1(data1.frame ==f,:);
    
    if(length(EC.frame) >100 && length(SC.frame) >100)
       temp1 = pdist2(EC{:,35:36},EC{:,35:36});
       temp1 = temp1(:);
       cell2cell{f,1} = temp1;
       neighbor(f,1) = sum(temp1<cutoff);

       temp1 = pdist2(SC{:,35:36},SC{:,35:36});
       temp1 = temp1(:);
       cell2cell{f,2} = temp1;
       neighbor(f,2) = sum(temp1<cutoff);

       temp1 = pdist2(EC{:,35:36},SC{:,35:36});
       temp1 = temp1(:);
       cell2cell{f,3} = temp1;
       neighbor(f,3) = sum(temp1<cutoff);

       temp1 = pdist2(EC{:,35:36},AC{:,35:36});
       temp1 = temp1(:);
       cell2cell{f,4} = temp1;
       neighbor(f,4) = sum(temp1<cutoff);

       temp1 = pdist2(SC{:,35:36},AC{:,35:36});
       temp1 = temp1(:);
       cell2cell{f,5} = temp1;
       neighbor(f,5) = sum(temp1<cutoff);

       %S1 = datasample(AC,length(EC.frame));
       %S2 = datasample(AC,length(SC.frame));
       %cell2cell{f,2} = pdist2(S1{:,35:36},S2{:,35:36});
       
       %temp1 = cell2mat(cell2cell);
       
    end
    neighbor(f,6) = length(EC.frame);
    neighbor(f,7) = length(SC.frame);
    neighbor(f,8) = length(AC.frame);

end

neighborT = array2table(neighbor);
neighborT.Properties.VariableNames = {'EC2EC','SC2SC','EC2SC','EC2AC','SC2AC','ECn','SCn','ACn'};
neighborT.frame = (1:maxframe)';

neighborT.nEC2E2 = neighborT.EC2EC ./ neighborT.EC2AC;
neighborT.nSC2SC = neighborT.SC2SC ./ neighborT.SC2AC;
neighborT.nEC2SC = neighborT.EC2SC ./ min(neighborT.EC2AC,neighborT.SC2AC);

figure,boxplot(neighborT{:,10:12});
title(myName);
set(gca,'xtick',1:3);
set(gca,'xticklabels',{'E2E','S2S','E2S'});




    
   