%% Analysis NeoSTAR pilot
%  Jerry Lin 2024/02/10


%% Import data

dir1 = dir;
dir1 = struct2table(dir1);
dir1 = dir1.name;

filelist = dir1;
filelist = filelist(3:end);
slideName = cellfun(@(X) X(1:8),filelist,'UniformOutput',false);
 
%% Correct the split-over from TROP2 to CD3e/PD_1

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));

    data1.CD3ep = data1.CD3ep & ~data1.TROP2p;
    data1.PD_1p = data1.PD_1p & ~data1.TROP2p;

    eval(strcat('data',slideName{i},'=data1;'));
end
toc;

%% Double gating

for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));

    for j=1:size(doubleGates,1)
        gate1 = strcat(doubleGates{j,1},'p');
        gate2 = strcat(doubleGates{j,2},'p');
        gatename = strcat(gate1,gate2);
        data1{:,gatename}=data1{:,gate1} & data1{:,gate2};
    end
    eval(strcat('data',slideName{i},'=data1;'));
end
toc;

%% Generate summary tables

% CycIF_alldata;

tic;

for i = 1:maxT
    marker1 = strcat('topic',num2str(i));
    alldata{:,marker1}= false(size(alldata,1),1);
    alldata{alldata.topics==i,marker1} = true;
end
toc;

sumAll = varfun(@mean,alldata,'GroupingVariables','slideName');
sumAll = join(sumAll,slideINFO,'Keys','slideName');
toc;

sumTumor = varfun(@mean,alldata(alldata.Region==1,:),'GroupingVariables','slideName');
sumTumor = join(sumTumor,slideINFO,'Keys','slideName');
toc;

sumBorder = varfun(@mean,alldata(alldata.Region==2,:),'GroupingVariables','slideName');
sumBorder = join(sumBorder,slideINFO,'Keys','slideName');
toc;

%% Boxplot (all intensity)

sum2 = sumAll;
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp{i});
    subplot(3,6,i);
    myboxplot2d(log(sum2{:,marker1}),sum2.Response);
    ylabel('Intensity(log)');
    title(labelp{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);


%% Boxplot (all gated markers)

sum2 = sumTumor;
%sum2 = sum2(~ismember(sum2.Group,'N/A'),:);
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    subplot(3,6,i);
    myboxplot2d(sum2{:,marker1}*100,sum2.Response);
    ytickformat('percentage');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);

%% forestplot (all gated markers)

sum2 = sum;
title1 = 'Tumor only';
odd_ratios = zeros(length(labelp),1);
CI = zeros(length(labelp),2);

for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
end
myforestplot(labelp3,odd_ratios,CI);
 
xlabel('fold difference');
title('Non-Responders/Responders');

%% Boxplot All double gates

sum3 = sumAll;
sum3 = sum3(~ismember(sum3.Group,'N/A'),:);
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(doubleGates,1)
    marker1 = strcat('mean_',doubleGates{i,1},'p',doubleGates{i,2},'p');
    subplot(3,8,i);
    myboxplot2d(sum3{:,marker1}*100,sum3.Group);
    ytickformat('percentage');
    title(strcat(doubleGates{i,1},'+',doubleGates{i,2},'+'),'Interpreter','none');
end
sgtitle(title1);
set(gcf,'color','w');
%% Boxplot (all topics)

sum2 = sumAll;
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:maxT
    marker1 = strcat('mean_topic',num2str(i));
    subplot(3,6,i);
    myboxplot2d(sum2{:,marker1}*100,sum2.Response);
    ytickformat('percentage');
    title(strcat('Topic',num2str(i)));
end

set(gcf,'color','w');
sgtitle(title1);

%% forestplot (all topics)

sum2 = sumAll;
odd_ratios = zeros(length(labelp),1);
CI = zeros(length(labelp),2);
names = cell(length(labelp),1)
for i = 1:maxT
    marker1 = strcat('mean_topic',num2str(i));
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
    names(i)={strcat('Topic',num2str(i))};
end
myforestplot(names,odd_ratios,CI);
 
xlabel('fold difference');
title('Non-Responders/Responders');


%% All markers (IFMs)

tic;
markers = {'CD3e','CD8a','CD45','CD11c','CD68','CD163','CD4','CD20','SMA','FOXP3','CD31','PD_1','TROP2','Pan_CK','Ki_67'};

markers = markers';

temp1 = strcat(markers,'_R1');
temp2 = strcat(markers,'_R2');
allmarkers = vertcat(temp1,temp2);
toc;


%% Generate sumAllset (IFMs)

tic;

sumAll = sortrows(sumAll,"slideName","ascend");
sumBorder = sortrows(sumBorder,"slideName","ascend");
sumTumor = sortrows(sumTumor,"slideName","ascend");

sum1 = sumTumor(:,strcat('mean_',markers','p'));
sum2 = sumBorder(:,strcat('mean_',markers','p'));

sum1.Properties.VariableNames = strcat(markers,'_R1');
sum2.Properties.VariableNames = strcat(markers,'_R2');
sumAllset = horzcat(sumAll,sum1,sum2);
toc;


%% Extensive test for all combinations (IFMs)

k = 4;
list1 = 1:length(allmarkers);
test1 = nchoosek(list1,k);

maxR = 0;
maxI = 0;
f = waitbar(0, 'Starting');
arrayR = zeros(size(test1,1),1);
arrayI = false(size(test1,1),length(allmarkers));

n = size(test1,1);
sum2 = sumAllset;
cutvalue = 50;

for i = 1:size(test1,1)
    markers = allmarkers(test1(i,:));

    tempScore = zeros(size(sum2,1),1);
    
    for j = 1:length(markers)
        list1 = sum2{:,markers{j}};
        cutoff1 = prctile(list1,cutvalue);
        tempScore = tempScore + (list1 > cutoff1);
    end

    r = log(1/myttest2(tempScore,sum2.Response));

    if r > maxR
        maxR = r;
        maxI = i;
    end
    arrayR(i) = r;
    arrayI(i,test1(i,:)) = true;
    
    waitbar(i/n, f, sprintf('Progress: %d %%', floor(i/n*100)));
end
close(f);

%% heatmap of all combinations (IFMs)

figure('units','normalized','outerposition',[0 0 1 1]);

[arrayRs, index]=sortrows(arrayR,'descend');
arrayIs = arrayI(index,:);
i = find(arrayRs>4.6,1,'last');

subplot(4,1,1);
bar(arrayRs);
ylabel('log(1/p)');
ylims = ylim;
line([i i],ylims,'Color','r','LineWidth',2);
ylim(ylims);
set(gca,'xtick',[]);

subplot(4,1,2:4)
imagesc(arrayIs');
colormap(cool);
ylabels = regexprep(allmarkers,'norm_','');
ylabels = regexprep(ylabels,'_R1','_T');
ylabels = regexprep(ylabels,'_R2','_B');
xlims = xlim;
ylims = ylim;
line(xlims,[mean(ylims),mean(ylims)],'Color','k','LineWidth',2);

set(gca,'ytick',1:length(allmarkers));
set(gca,'yticklabels',ylabels);
set(gca,'ticklabelinterpreter','none');
set(gca,'xtick',[]);
xlabel('All combinations');
set(gcf,'color','w');


%% Check result (IFMs)

i =3;
markers = allmarkers(arrayIs(i,:))

tempScore = zeros(size(sumAllset,1),1);

for j = 1:length(markers)
    list1 = sumAllset{:,markers{j}};
    cutoff1 = prctile(list1,cutvalue);
    tempScore = tempScore + (list1 > cutoff1);
end

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);

subplot(1,2,1)
myboxplot2d(tempScore,sumAllset.Response);
ylim([0 5]);
set(gca,'ytick',0:5);
title(strcat('Combination:',num2str(i)));

list1 = histc(tempScore(ismember(sumAllset.Response,'NR')),0:4);
list1 = list1 / sum(list1);
list2 = histc(tempScore(ismember(sumAllset.Response,'R')),0:4);
list2 = list2 / sum(list2);

list1 = horzcat(list1,list2);
subplot(1,2,2);
bar(list1*100);
legend({'NR','R'},'Location','northwest');
set(gca,'xticklabel',0:4);
ylabel('percentage');
ytickformat('percentage');
xlabel('scores');

set(gcf,'color','w')

%% Segmentation and assign TLS for all slides

minD = 100;
minSize = 150;

tic;
for i =1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data0 = eval(strcat('data',slideName{i}));
    
    data0.cellID = (1:size(data0,1))';
    data1 = data0(data0.TLS>0,:);

    pc1 = pointCloud(double([data1.Xt,data1.Yt,ones(size(data1,1),1)]));
    label1=pcsegdist(pc1,minD);
    data1.label1 = label1;

    table1 = tabulate(data1.label1);
    table1 = table1(table1(:,2)> minSize,:);
    data1 = data1(ismember(data1.label1,table1(:,1)),:);
    table1(:,4) = (1:size(table1,1))';
    data1.label2 = zeros(size(data1,1),1);
    for j = 1:size(table1,1)
        data1.label2(data1.label1==table1(j,1)) = table1(j,4);
    end

    
    data0.TLSID = zeros(size(data0,1),1);
    data0.TLSID(data1.cellID) = data1.label2;
    eval(strcat('data',slideName{i},'=data0;'));
    tabulate(data0.TLSID);
    toc;
    %clear data0 data1;
end

%% Calculate aggregation index?

name1 = 'CD163';
marker1 = strcat(name1,'p');

data1 = dataLSP20621;
data2 = data1(data1{:,marker1},:);
data3 = datasample(data1,size(data2,1),'replace',false);
[~,d2] = knnsearch([data2.Xt,data2.Yt],[data2.Xt,data2.Yt],'k',2);
[~,d3] = knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],'k',2);
AI = 1- mean(d2(:,2))/mean(d3(:,2))

%% Calculate aggregation index for all samples/markers

allaggindex = zeros(length(slideName),length(labelp));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s}));
    data1 = eval(strcat('data',slideName{s}));

    parfor i = 1:length(labelp)
        name1 = labelp{i};
        marker1 = strcat(name1,'p');
        data2 = data1(data1{:,marker1},:);
        data3 = datasample(data1,size(data2,1),'replace',false);
        [~,d2] = knnsearch([data2.Xt,data2.Yt],[data2.Xt,data2.Yt],'k',2);
        [~,d3] = knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],'k',2);
        allaggindex(s,i) = 1- mean(d2(:,2))/mean(d3(:,2));
    end
    toc;
end

%% violinplot for topics

% sumX1 = varfun(@mean,alldata,'GroupingVariables',{'slideName','topics'});
figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:16
    marker1 = strcat('mean_',labelp2{i});
    subplot(8,2,i);
    violinplot(sumX1{:,marker1}*100,sumX1.topics);
    ytickformat('percentage');
    ylabel(labelp3{i});
end

%% Ratio test

array1 = zeros(length(labelp),length(labelp));

for i = 1:length(labelp)
    list1 = sumAll{:,strcat('mean_',labelp2{i})};
    for j = 1:length(labelp)
        list2 = sumAll{:,strcat('mean_',labelp2{j})};
        array1(i,j)=myttest2(list1./list2,sumAll.Response);
    end
end


%% Douple positive

array1 = zeros(length(labelp),length(labelp));
tic;
for i = 1:length(labelp)
    parfor j = 1:length(labelp)
        data1 = alldata;
        list1 = alldata{:,labelp2{i}};
        list2 = alldata{:,labelp2{j}};
        data1.temp1 = list1 & list2;
        sumTemp = varfun(@mean,data1,'GroupingVariables','slideName');
        array1(i,j) = myttest2(sumTemp.mean_temp1,sumAll.Response);
    end
    toc;
end

%% All interaction (by gates)

arrayIntGates = zeros(length(labelp2),length(labelp2),length(slideName));

tic;
for s = 1:length(slideName)
    data0 = eval(strcat('data',slideName{s}));
    disp(strcat('Processing:',slideName{s}));
    
    for i = 1:length(labelp2)
        marker1 = labelp2{i};
        parfor j = 1:length(labelp2)
                data1 = data0;
                marker2 = labelp2{j}
                data2 = data1(data1{:,marker1},:);
                data3 = data1(data1{:,marker2},:);
                data2r = datasample(data1,size(data2,1));
                data3r = datasample(data1,size(data3,1));
                if (size(data2,1)>0 & size(data3,1)>0)
                    [~,d23]= knnsearch([data3.Xt,data3.Yt],[data2.Xt,data2.Yt]);
                    [~,d32]= knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt]);
                    [~,d23r]= knnsearch([data3r.Xt,data3r.Yt],[data2.Xt,data2.Yt]);
                    [~,d32r]= knnsearch([data2r.Xt,data2r.Yt],[data3.Xt,data3.Yt]);
                    arrayIntGates(i,j,s)= mean(d23r)*mean(d32r)/mean(d23)/mean(d32);
                end
         end
         toc;
    end
end

%% Gate2gate inteaciton Boxplot (all)

i = 17;

figure('units','normalized','outerposition',[0 0 1 1]);
for j = 1:length(labelp2)
    subplot(3,6,j);
    myboxplot2d(squeeze(arrayIntGates(i,j,:)),sumAll.Response);
    %ytickformat('percentage');
    title(labelp3{j});   
end
set(gcf,'color','w');
sgtitle(labelp3{i})

%% P values for all interaction

arrayIntPvalue = zeros(length(labelp2),length(labelp2));
tic;
for i = 1:length(labelp2)
    for j = 1:length(labelp2)
            list1 = arrayIntGates(i,j,ismember(sumAll.Response,'R'));
            list2 = arrayIntGates(i,j,ismember(sumAll.Response,'NR'));
            arrayIntPvalue(i,j)=myttest2(squeeze(arrayIntGates(i,j,:))*100,sumAll.Response);
    end
    toc;
end

%% All pairwise interaction (second method)

cutoff = 20;

allint1 = zeros(length(labelp2),length(labelp2),length(slideName));
allint1r = zeros(length(labelp2),length(labelp2),length(slideName));

tic;
for s = 1:length(slideName)
    data0 = eval(strcat('data',slideName{s}));
    disp(strcat('Processing:',slideName{s}));

    for i = 1:length(labelp2)
        marker1 = labelp2{i};
        parfor j = 1:length(labelp2)
            data1 = data0;
            marker2=labelp2{j};
            dataP1 = data1(data1{:,marker1},:);
            dataP2 = data1(data1{:,marker2},:);
            dataP2r = datasample(data1,size(dataP2,1));
            [~,d1]=knnsearch([dataP2.Xt,dataP2.Yt],[dataP1.Xt,dataP1.Yt],'k',2);
            [~,d2]=knnsearch([dataP2r.Xt,dataP2r.Yt],[dataP1.Xt,dataP1.Yt],'k',2);
            if ~isempty(d1)
                allint1(i,j,s) = mean(d1(:,2)<cutoff);
                allint1r(i,j,s) = mean(d2(:,2)<cutoff);
            end
        end
    end
    toc;
end

%% P values for all interaction (second method)

temp1 = allint1r ./ allint1;
allIntPvalue = zeros(length(labelp2),length(labelp2));
allIntRatio = zeros(length(labelp2),length(labelp2));
tic;
for i = 1:length(labelp2)
    for j = 1:length(labelp2)
        [allIntPvalue(i,j),allIntRatio(i,j)] = myttest2(squeeze(temp1(i,j,:)),sumAll.Response);
    end
    toc;
end

%% Selected inteaciton Boxplot (Second method)

i = find(ismember(labelp,'CD68'));
temp1 = allint1r ./ allint1;

figure('units','normalized','outerposition',[0 0 1 1]);
for j = 1:length(labelp2)
    subplot(3,6,j);
    myboxplot2d(squeeze(temp1(i,j,:)),sumAll.Response);
    %ytickformat('percentage');
    title(labelp3{j});   
end
set(gcf,'color','w');
sgtitle(labelp3{i})

%% Selected inteaciton  single Boxplot (Second method) (Fig4B)

i = find(ismember(labelp,'CD68'));
j = find(ismember(labelp,'CD20'));

temp1 = allint1r ./ allint1;

figure,myboxplot2dj(squeeze(temp1(i,j,:)),sumAll.Response);
title(strcat(labelp3{i},'::',labelp3{j}));  
ylabel('Normalized Interaction');
set(gcf,'color','w');
set(gca,'xticklabels',{'pCR','RD'});

%% circular graph of responder 
temp1 = allint1r ./ allint1;
temp2 = temp1(:,:,ismember(sumAll.Response,'R'));
temp2 = mean(temp2,3);
temp2(isinf(temp2))=0;
temp2 = temp2-2;
temp2(temp2<1.75) = 0;
figure,circularGraph(temp2,'Colormap',repmat([1,0,0],17,1),'Label',labelp4);
text(-0.35,-1.7,'Responder','FontSize',16);

%% circular graph of non-responder 
temp1 = allint1r ./ allint1;
temp2 = temp1(:,:,~ismember(sumAll.Response,'R'));
temp2 = mean(temp2,3);
temp2(isinf(temp2))=0;
temp2(isnan(temp2))=0;
temp2 = temp2-2;
temp2(temp2<1.5) = 0;
figure,circularGraph(temp2,'Colormap',repmat([0,0,1],17,1),'Label',labelp4);
text(-0.55,-1.7,'Non-Responder','FontSize',16);


%% Import TROP2-Ecad Overlap data;

arrayTrop2Ecad = zeros(length(slideName),3);

for i =1:length(slideName)
    filename = strcat('Results-',slideName{i},'.csv');
    temp1 = readtable(filename);
    arrayTrop2Ecad(i,1) = temp1.Mean(1);
    arrayTrop2Ecad(i,2) = temp1.Mean(2);
    arrayTrop2Ecad(i,3) = temp1.Mean(2)/(temp1.Mean(1)+temp1.Mean(2));
end

%% Calculate area 

arrayArea = zeros(length(slideName),1);

tic;
for s = 1:length(slideName)
    data0 = eval(strcat('data',slideName{s}));
    disp(strcat('Processing:',slideName{s}));
    sum0 = varfun(@mean,data0,'GroupingVariables','frame');
    arrayArea(s) = sum(sum0.GroupCount>20)*0.15*0.15;
    toc;
end

%% Boxplot (all gated markers, adjusted by area)

sum2 = sumAll;
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    subplot(3,6,i);
    myboxplot2d(sum2{:,marker1} .* sum2.GroupCount ./arrayArea,sum2.Response);
    %ytickformat('percentage');
    ylabel('cells/mm^2');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);

%% Boxplot All double gates (adjusted by area)

sum3 = sumAll;
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(doubleGates,1)
    marker1 = strcat('mean_',doubleGates{i,1},'p',doubleGates{i,2},'p');
    subplot(3,7,i);
    myboxplot2d(sum3{:,marker1}.*sum3.GroupCount./arrayArea,sum3.Response);
    %ytickformat('percentage');
    ylabel('cells/mm^2');
    title(strcat(doubleGates{i,1},'+',doubleGates{i,2},'+'),'Interpreter','none');
end
sgtitle(title1);
set(gcf,'color','w');

%% assign cell idx

data1 = alldata;

%label_cells = tableCelltype.Properties.VariableNames;
%label_cells = label_cells(2:17);

array1 = data1{:,labelp2};
idx1 = zeros(size(array1,1),1);

tic;
parfor i =1:size(array1,1)
    temp1 = array1(i,:);
    temp1 = num2str(temp1);
    temp1(isspace(temp1))='';
    idx1(i) = bin2dec(temp1);
end
toc;

data1.cell_idx = idx1;

clear array1 temp1 idx1;

%% check cell type

temp1 = dec2bin(16520,17);
temp1 = temp1';
temp1 = str2num(temp1);
temp1 = temp1 > 0;
labelp3(temp1)

%% check tumor heterogenity

labelX = {'E_cadherin','Pan_CK','TROP2'};
labelX2 = strcat(labelX,'p');
labelX3 = strcat(labelX,'+');

data1 = dataLSP20621;

array1 = data1{:,labelX2};
idx1 = zeros(size(array1,1),1);

tic;
parfor i =1:size(array1,1)
    temp1 = array1(i,:);
    temp1 = num2str(temp1);
    temp1(isspace(temp1))='';
    idx1(i) = bin2dec(temp1);
end
toc;

data1.tumor_idx = idx1;

%% check tumor type

temp1 = dec2bin(7,3);
temp1 = temp1';
temp1 = str2num(temp1);
temp1 = temp1 > 0;
strcat(labelX3(temp1)')

%% Apply tumor heterogenity to all samples

labelX = {'E_cadherin','Pan_CK','TROP2'};
labelX2 = strcat(labelX,'p');
labelX3 = strcat(labelX,'+');

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s}));
    data1 = eval(strcat('data',slideName{s}));
    
    array1 = data1{:,labelX2};
    idx1 = zeros(size(array1,1),1);
    
    parfor i =1:size(array1,1)
        temp1 = array1(i,:);
        temp1 = num2str(temp1);
        temp1(isspace(temp1))='';
        idx1(i) = bin2dec(temp1);
    end
    toc;
    data1.tumor_idx = idx1;
    eval(strcat('data',slideName{s},'=data1;'));
end

%% Generate sumT 
tic;
for i = 0:7
    marker1 = strcat('tumor_idx',num2str(i));
    alldata{:,marker1} = false(size(alldata,1),1);
    alldata{alldata.tumor_idx==i,marker1} = true;
end

sumT = varfun(@mean,alldata,'GroupingVariables','slideName');
sumT = join(sumT,slideINFO,'keys','slideName');
toc;

%% check all tumor_idx

figure('units','normalized','outerposition',[0 0 1 1]);

for i  = 0:7
    subplot(2,4,i+1);
    name1 = strcat('tumor_idx',num2str(i));
    marker1 = strcat('mean_',name1);
    myboxplot2d(sumT{:,marker1} .* 100,sumT.Response);
    ytickformat('percentage');
    if i == 0
        title1 = 'All neg.';
    else
        title1 = listTumor{i};
    end
    title(title1);
end

%% Apply ROI to all

dir1 = dir;
dir1 = struct2table(dir1);
listROI = dir1.name;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));
    data1.ROI = zeros(size(data1,1),1);
    %data1.ROIname = repmat({'none'},size(data1,1),1);
    filename = listROI{i};
    if exist(filename,'file')
        data1 = CycIF_assignOmeroROIsingle2(data1,filename,0.325,false,'ROI',true);
    end
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Output Topic maps

tic;
for i =1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));

    figure('units','normalized','outerposition',[0 0 1 1]);

    CycIF_tumorview(data1,'topics',5,0);
    daspect([1 1 1]);
    set(gcf,'color','w');
    title(slideName{i});
    filename = strcat(slideName{i},'.png');
    saveas(gcf,filename);
    close;
    toc;
end

%% Calculate TROP2 autocorrelation in E_cadherin+ cells

max_n = 2000;
r_max = zeros(length(slideName),1);
d_quad = zeros(length(slideName),1);
r_all = zeros(length(slideName),max_n);

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1 = data1(data1.E_cadherinp,:);
    
    [idx,d] = knnsearch([data1.Xt,data1.Yt],[data1.Xt,data1.Yt],'k',max_n);
    test1 = data1.TROP2(idx);
    test2 = corr(test1);
    r_max(i) = test2(1,max_n);
    r_all(i,:)= test2(1,:);
    d_quad(i) = find(test2(1,:)>(1-r_max(i))/4+r_max(i),1,'last');
    toc;
end

%% Test all distance;


all_r_pvalue = zeros(max_n,1);

tic;
for i = 1:max_n
    all_r_pvalue(i) = myttest2(r_all(:,i),sumAll.Response);
end
toc;

%% Plot correlation landscale

figure;

plot(smoothdata(r_all(ismember(sumAll.Response,'NR'),1:1445)','SmoothingFactor',0.2),'Color','b','LineWidth',1);
hold on;
plot(smoothdata(r_all(ismember(sumAll.Response,'R'),1:1445)','SmoothingFactor',0.2),'Color','r','LineWidth',1);

%% forest (all gated markers)

sum2 = sumAll;
%title1 = 'Whole sample';
odd_ratios = zeros(length(labelp),1);
CI = zeros(length(labelp),2);

for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
end
myforestplot(labelp3,odd_ratios,CI);
 
xlabel('fold difference');
title('RD/pCR (whole sample)','FontSize',18);

%% calculate all cluster score

arrayClusterScore = zeros(length(labelp),length(labelp),length(slideName));

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    for i = 1:length(labelp)
        marker1 = labelp{i};
        for j = 1:length(labelp)
            marker2 = labelp{j};
           
            gate1 = strcat(marker1,'p');
            gate2 = strcat(marker2,'p');
            
            data2 = data1(data1{:,gate1},:);
            data3 = data1(data1{:,gate2},:);
            data4 = datasample(data1,size(data3,1));
            [~,d1] = knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],'K',2);
            [~,d2] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'K',2);
            arrayClusterScore(i,j,s) = mean(d2(:,2))/mean(d1(:,2));
        end
    end
    toc;
end

%% check cluster score per sample

i = 1;

figure,imagesc(arrayClusterScore(:,:,i));
colormap(jet);
colorbar;

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
caxis([0 3]);
title(slideName{i},'FontSize',18);
set(gcf,'color','w');

%% Test all pvalues

arrayClusterPvalues = zeros(length(labelp),length(labelp));

tic;
for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayClusterPvalues(i,j)=myttest2(squeeze(arrayClusterScore(i,j,:)),sumAll.Response);
    end
end
toc;

figure,imagesc(arrayClusterPvalues<0.01);

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);

%% check 

name1 = 'CD8a';
name2 = 'SMA';

i = find(ismember(labelp,name1));
j = find(ismember(labelp,name2));

figure,myboxplot2d(squeeze(arrayClusterScore(i,j,:)),sumAll.Response);
title(strcat(name1,'-',name2),'FontSize',18);

%% Define TROP2 high cells

for i = 1:size(gateTable,1)
    disp(strcat('Processing:',gateTable.slideName{i},'...',num2str(i)));
    data1 = eval(strcat('data',gateTable.slideName{i}));
    data1.TROP2h = data1.TROP2 > exp(gateTable.TROP2(i)+1);
    eval(strcat('data',gateTable.slideName{i},'=data1;'));
end

%% Define TROP2c (TROP2 category, 0=other 1=negative 2=low 3=high)

for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data1.TROP2c = double(data1.TROP2p);
    data1.TROP2c(data1.TROP2c==1) = 2;
    data1.TROP2c(data1.TROP2h) = 3;
    data1.TROP2l = data1.TROP2c ==2;
    data1.TROP2n = data1.E_cadherinp & ~data1.TROP2p;
    data1.TROP2c(data1.TROP2n) = 1;
    tabulate(data1.TROP2c);
    eval(strcat('data',gateTable.slideName{s},'=data1;'));
end


%% Calculate Entropy (E_cadherin+ & TROP2+)

sfactor = 25;

arrayEntropyEcad = zeros(length(slideName),1);
arrayEntropyTROP2P = zeros(length(slideName),1);
arrayEntropyTROP2L = zeros(length(slideName),1);
arrayEntropyTROP2H = zeros(length(slideName),1);
arrayEntropyTROP2N = zeros(length(slideName),1);
arrayEntropyTROP2C = zeros(length(slideName),1);

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    marker1 = 'E_cadherinp';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyEcad(s) = entropy(temp1);
    marker1 = 'TROP2p';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyTROP2P(s) = entropy(temp1);
    marker1 = 'TROP2h';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyTROP2H(s) = entropy(temp1);
    marker1 = 'TROP2l';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyTROP2L(s) = entropy(temp1);
    marker1 = 'TROP2n';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyTROP2N(s) = entropy(temp1);
    marker1 = 'TROP2c';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyTROP2C(s) = entropy(temp1);

    toc;
end


%% calculate all cluster score

arrayIntTROP2N = zeros(length(labelp),length(slideName));

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data2 = data1(data1.TROP2n,:);

    for j = 1:length(labelp)
        marker2 = labelp{j};
        gate2 = strcat(marker2,'p');
        data3 = data1(data1{:,gate2},:);
        data4 = datasample(data1,size(data3,1));
        [~,d1] = knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],'K',1);
        [~,d2] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'K',1);
        arrayIntTROP2N(j,s) = mean(d2)/mean(d1);
    end
    toc;
end

%------ Check results------

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    subplot(3,6,i);
    myboxplot2d(arrayIntTROP2N(i,:),sumAllsample.Response);
    ylabel('Interaction score');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('TROP2 negative');

%%  Entropy (TROP2 category)

arrayEntropyTROP2C = zeros(length(slideName),1);

sfactor = 25;
tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    marker1 = 'TROP2c';
    temp1 = accumarray([floor(data1.Xt/sfactor)+1,floor(data1.Yt/sfactor)+1],data1{:,marker1},[floor(max(data1.Xt)/sfactor)+1,floor(max(data1.Yt)/sfactor)+1]);
    arrayEntropyTROP2C(s) = entropy(temp1);
    toc;
end

figure,myboxplot2d(arrayEntropyTROP2C,sumAllsample.Response);

%% Assign VESp and CD31pVESn

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    data1.VESp = data1.VesID>0;
    data1.CD31pVESn = data1.CD31p & ~data1.VESp;
    data1.CD31pVESp = data1.CD31p & data1.VESp;
    eval(strcat('data',name1,'=data1;'));
end
toc;

%% Distance between Immune cells and VESp

%arrayDisVESp = zeros(length(labelp),length(slideName));
%arrayDisCD31pVESn = zeros(length(labelp),length(slideName));
%arrayDisCD31pVESp = zeros(length(labelp),length(slideName));
arrayDisSMApVESp = zeros(length(labelp),length(slideName));
arrayDisCD31p = zeros(length(labelp),length(slideName));

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data2 = data1(data1.SMAp & data1.VESp,:);
    data3 = data1(data1.CD31p,:);
    for i = 1:length(labelp2)
        data4 = data1(data1{:,labelp2{i}},:);
        data5 = datasample(data1,size(data4,1));
        
        [~,d1] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'k',2);
        [~,d2] = knnsearch([data2.Xt,data2.Yt],[data5.Xt,data5.Yt],'k',2);
        arrayDisSMApVESp(i,s)= mean(d2(:,2))/mean(d1(:,2));

        [~,d1] = knnsearch([data3.Xt,data3.Yt],[data4.Xt,data4.Yt],'k',2);
        [~,d2] = knnsearch([data3.Xt,data3.Yt],[data5.Xt,data5.Yt],'k',2);
        arrayDisCD31p(i,s)= mean(d2(:,2))/mean(d1(:,2));
    end
    toc;
end

%% Boxplot (arrayDisCD31pVESp)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisCD31pVESp(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To Vessel CD31+');

%% Boxplot (arrayDisCD31pVESn)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisCD31pVESn(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To non-Vessel CD31+');

%% Boxplot (arrayDisCD31p)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisCD31p(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To all CD31+');

%% Boxplot (arrayDisSMApVESp)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisSMApVESp(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To Vessel SMA+');


%% Boxplot (arrayDisVESp versus arrayDisCD31pVESn)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    list1 = arrayDisVESp(i,:);
    list2 = arrayDisCD31pVESn(i,:);
    list3 = vertcat(list1,list2);list3 = list3';

    myboxplot2d(list3);
    set(gca,'xticklabels',{'VES.','Imm. CD31'});
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');

%% Distance between Immune cells and TROP2 high/low

arrayDisTROP2h = zeros(length(labelp),length(slideName));
arrayDisTROP2l = zeros(length(labelp),length(slideName));

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data2 = data1(data1.TROP2h,:);
    data3 = data1(data1.TROP2l,:);
    for i = 1:length(labelp2)
        data4 = data1(data1{:,labelp2{i}},:);
        data5 = datasample(data1,size(data4,1));
        
        [~,d1] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'k',1);
        [~,d2] = knnsearch([data2.Xt,data2.Yt],[data5.Xt,data5.Yt],'k',1);
        arrayDisTROP2h(i,s)= mean(d2)/mean(d1);

        [~,d1] = knnsearch([data3.Xt,data3.Yt],[data4.Xt,data4.Yt],'k',1);
        [~,d2] = knnsearch([data3.Xt,data3.Yt],[data5.Xt,data5.Yt],'k',1);
        arrayDisTROP2l(i,s)= mean(d2)/mean(d1);
    end
    toc;
end

%% Boxplot (arrayDisTROP2h)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisTROP2h(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To TROP2 high');

%% Boxplot (arrayDisTROP2l)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisTROP2l(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To TROP2 low');

%% Boxplot (arrayDisTROP2h versus arrayDisTROP2low)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    list1 = arrayDisTROP2h(i,:);
    list2 = arrayDisTROP2l(i,:);
    list3 = vertcat(list1,list2);list3 = list3';

    myboxplot2d(list3);
    set(gca,'xticklabels',{'High','Low'});
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To TROP2+');

%%  Ki67+ Tumor (pCR/RD) (Fig2J)

sum1 = sumTumor;

figure,myboxplot2dj(sum1.mean_E_cadherinpKi_67p*100,sum1.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});
title('TROP2+ tumor cells','FontSize',16);

%%  Ki67+ Tumor (Waterfall plot)

sum1 = sumTumor;
name1 = 'TROP2+ Tumor';
marker1 = 'mean_E_cadherinpTROP2p';

sum1 = sortrows(sum1,marker1,"ascend");

figure('units','normalized','outerposition',[0.5 0 0.5 1]);

b= barh(sum1{:,marker1}*100,'FaceColor',[0,0,1],'EdgeColor',[1,1,1],'LineWidth',0.1);
hold on;
b.FaceColor = 'flat';
xtickformat('percentage');

for i = 1:size(sum1,1)
    b.CData(i,:) = [0 0 1];
end

for i = 1:size(sum1,1)
    if ismember(sum1.Response(i),'R')
        b.CData(i,:) = [1 0 0];
    end
end


h = zeros(2, 1);
h(1) = barh(NaN,'Red');
h(2) = barh(NaN,'Blue');
legend(h, {'pCR','RD'},'Location','southeast','FontSize',16);

title(name1,'FontSize',18,'Interpreter','none');
set(gcf,'color','w');
set(gca,'ytick',[]);
ylabel('samples');

%%  Ki67+ CD8 cells (Fig3G)

sum1 = sumTumor;

figure,myboxplot2dj(sum1.mean_CD8apKi_67p./sum1.mean_CD8ap*100,sum1.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});
title('%Ki67+ in CD8+ cells','FontSize',16);


%%  PD1+ CD4/CD8 cells

sum1 = sumAll;

figure,myboxplot2d(sum1.mean_CD4pPD_1p./sum1.mean_CD4p*100,sum1.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});
title('PD1+ CD8 cells','FontSize',16);


%%  Treg cells

sum1 = sumAll;

figure,myboxplot2d(sum1.mean_CD4pFOXP3p./sum1.mean_CD4p,sum1.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});
title('Treg (CD4+/FOXP3+)','FontSize',16);

%%  CD8a/CD4 ratio (Fig3X??)

sum1 = sumAll_new;

figure,myboxplot2dj(sum1.mean_CD8ap./sum1.mean_CD4p,sum1.Response);
%ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',12);
title('CD8/CD4 ratio','FontSize',16);

%% forest (all gated markers, only immune)

sum2 = sumImmune;
title1 = 'Not-Epi';
odd_ratios = zeros(length(labelp),1);
CI = zeros(length(labelp),2);

for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
end
myforestplot(labelp3,odd_ratios,CI);
 
xlabel('fold difference');
title(strcat('RD/pCR (',title1,')'),'FontSize',18);

%% generate double gates list
labeldp = cell(length(doubleGates),1);
labeldp2 = cell(length(doubleGates),1);

for i = 1:size(doubleGates,1)
    labeldp{i} = strcat(doubleGates{i,1},'p',doubleGates{i,2},'p');
    labeldp2{i} = strcat(doubleGates{i,1},'+',doubleGates{i,2},'+');
end

%% forest (all double gates)

%sum2 = sumTumor_new;
sum2 = sumNotEpi;
title1 = 'Not-Epi';

odd_ratios = zeros(length(labeldp),1);
CI = zeros(length(labeldp),2);

for i = 1:length(labeldp)
    marker1 = strcat('mean_',labeldp{i});
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
end
myforestplot(labeldp2,odd_ratios,CI);
 
xlabel('fold difference');
title(strcat('RD/pCR (',title1,')'),'FontSize',18);

%% All ratios
sum1 = sumAll;
allRatio = zeros(length(labelp),length(labelp),size(sum1,1));
allRatioP = zeros(length(labelp),length(labelp));

for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    list1 = sum1{:,marker1};
    for j = 1:length(labelp)
        marker2 = strcat('mean_',labelp2{j});
        list2 = sum1{:,marker2};
        allRatio(i,j,:)=list1./list2;
        allRatioP(i,j) = myttest2(list1./list2,sum1.Response);
    end
end

figure,imagesc(allRatioP<0.05);
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);


%% Neighborhood analysis

arrayNeighbor = zeros(length(labelp),length(labelp),length(slideName));
cutoff = 20;

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data0 = eval(strcat('data',name1));
    for i = 1:length(labelp2)
        data1 = data0(data0{:,labelp2{i}},:);
        parfor j= 1:length(labelp2)
            data2 = data0(data0{:,labelp2{j}},:);
            [~,d1]= knnsearch([data1.Xt,data1.Yt],[data2.Xt,data2.Yt],'k',2);
            [~,d2]= knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',2);
            arrayNeighbor(i,j,s)=sum(d1(:,2)<cutoff)*sum(d2(:,2)<cutoff)/size(data1,1)/size(data2,1);
        end
    end
    toc;
end

%% Neighborhood analysis (P value)

arrayNeighborP = zeros(length(labelp),length(labelp));

for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayNeighborP(i,j) = myttest2(squeeze(arrayNeighbor(i,j,:)),sumAll.Response);
    end
end

figure,imagesc(arrayNeighborP<0.05);

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);


%% Neighborhood analysis (Trop2 high/Low)

arrayNeighborTROP2h = zeros(length(labelp),length(slideName));
arrayNeighborTROP2l = zeros(length(labelp),length(slideName));

cutoff = 20;

for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data0 = eval(strcat('data',name1));
    data2 = data0(data0.TROP2h,:);
    data3 = data0(data0.TROP2l,:);

    parfor i = 1:length(labelp2)
        data1 = data0(data0{:,labelp2{i}},:);
        [~,d12]= knnsearch([data1.Xt,data1.Yt],[data2.Xt,data2.Yt],'k',2);
        [~,d21]= knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',2);
        arrayNeighborTROP2h(i,s)= sum(d12(:,2)<cutoff)*sum(d21(:,2)<cutoff)/size(data1,1)/size(data2,1);

        [~,d13]= knnsearch([data1.Xt,data1.Yt],[data3.Xt,data3.Yt],'k',2);
        [~,d31]= knnsearch([data3.Xt,data3.Yt],[data1.Xt,data1.Yt],'k',2);
        arrayNeighborTROP2l(i,s)= sum(d13(:,2)<cutoff)*sum(d31(:,2)<cutoff)/size(data1,1)/size(data3,1);
    end
    toc;
end


%% Boxplot (TROP2 high/low)

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayNeighborTROP2l(i,:),sumAll.Response);
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('TROP2 low','FontSize',18);


%% Neighborhood analysis (TROP2n)

arrayNeighborTROP2n = zeros(length(labelp),length(slideName));

cutoff = 20;

for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data0 = eval(strcat('data',name1));
    data2 = data0(data0.TROP2n,:);
    
    parfor i = 1:length(labelp2)
        data1 = data0(data0{:,labelp2{i}},:);
        [~,d12]= knnsearch([data1.Xt,data1.Yt],[data2.Xt,data2.Yt],'k',2);
        [~,d21]= knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',2);
        arrayNeighborTROP2n(i,s)= sum(d12(:,2)<cutoff)*sum(d21(:,2)<cutoff)/size(data1,1)/size(data2,1);
    end
    toc;
end

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayNeighborTROP2n(i,:),sumAll.Response);
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('TROP2 neg.','FontSize',18);

%% Neighborhood analysis (Ves & CD31)

arrayNeighborVes = zeros(length(labelp),length(slideName));
arrayNeighborCD31nv = zeros(length(labelp),length(slideName));

cutoff = 20;

for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data0 = eval(strcat('data',name1));
    data2 = data0(data0.VesID>0,:);
    data3 = data0(data0.CD31p & data0.VesID==0,:);

    parfor i = 1:length(labelp2)
        data1 = data0(data0{:,labelp2{i}},:);
        [~,d12]= knnsearch([data1.Xt,data1.Yt],[data2.Xt,data2.Yt],'k',2);
        [~,d21]= knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',2);
        arrayNeighborVes(i,s)= sum(d12(:,2)<cutoff)*sum(d21(:,2)<cutoff)/size(data1,1)/size(data2,1);

        [~,d13]= knnsearch([data1.Xt,data1.Yt],[data3.Xt,data3.Yt],'k',2);
        [~,d31]= knnsearch([data3.Xt,data3.Yt],[data1.Xt,data1.Yt],'k',2);
        arrayNeighborCD31nv(i,s)= sum(d13(:,2)<cutoff)*sum(d31(:,2)<cutoff)/size(data1,1)/size(data3,1);
    end
    toc;
end

%% Boxplot (Ves/CD31)

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayNeighborCD31nv(i,:),sumAll.Response);
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('non-Vessel CD31','FontSize',18);

%% violinplot for kmeans

% sumX1 = varfun(@mean,alldata,'GroupingVariables',{'slideName','topics'});
figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:18
    marker1 = strcat('mean_',labelp2{i});
    subplot(9,2,i);
    violinplot(sumT{:,marker1}*100,sumT.kmean25);
    ytickformat('percentage');
    ylabel(labelp3{i});
end

%% TLS counts (Fig4E)

% CycIF_Generalsegmentation_20240804
% minDis =50, minCount = 250;

% sumCD20p = varfun(@mean,alldata,'GroupingVariables',{'slideName','CD20pID'});
% sumCD20p = varfun(@mean,sumCD20p,'GroupingVariables','slideName');
% sumCD20p = join(sumCD20p,slideINFO,'keys','slideName');

figure,myboxplot2dj(sumCD20p.GroupCount-1,sumCD20p.Response);
title('TLS counts','FontSize',18);
set(gca,'xticklabels',{'pCR','RD'});
ylabel('number of TLS');

%% TLS counts in four group

sum1 = sumCD20p;
sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;

figure,myboxplot3d(sum1.GroupCount,sum1.Group3);
title('TLS counts','FontSize',14);
set(gca,'xticklabels',listGroup3);
ylabel('number of TLS');


%% TLS counts in three group

sum1 = sumCD20p;


figure,myboxplot4d(sum1.GroupCount-1,sum1.Group3);
title('TLS counts','FontSize',14);
set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('number of TLS');
%% TROP2 ITH (Fig2X??)

alldata.TROP2z = zscore(alldata.TROP2);
alltumor = alldata(alldata.E_cadherinp,:);

sumAll = sortrows(sumAll,"slideName","ascend");

figure,boxplot(zscore(alltumor.TROP2z),alltumor.slideName,'symbol','');
h = findobj(gca, 'Tag', 'Box');  % Find all box objects
set(gca,'xticklabels',sumAll.PatientID)

for j = 1:length(h)
    if ismember(sumAll.Response(j),'R')
        patch(get(h(j), 'XData'), get(h(j), 'YData'), [1,0,0], 'FaceAlpha', 0.5);
    else
        patch(get(h(j), 'XData'), get(h(j), 'YData'), [0,0,1], 'FaceAlpha', 0.5);
    end
end

title('TROP2 ITH');
ylabel('zscore(TROP2)');

%% TROP2 ITH (log intensity)

sumAll = sortrows(sumAll,"slideName","ascend");

figure,boxplot(log(alltumor.TROP2+5),alltumor.slideName,'symbol','');
ylim([5 11]);
h = findobj(gca, 'Tag', 'Box');  % Find all box objects
set(gca,'xticklabels',sumAll.PatientID)
ylabel('Intensity(log)')
title('TROP2 expression (ITH)','FontSize',18);

for j = 1:length(h)
    if ismember(sumAll.Response(j),'R')
        patch(get(h(j), 'XData'), get(h(j), 'YData'), [1,0,0], 'FaceAlpha', 0.5);
    else
        patch(get(h(j), 'XData'), get(h(j), 'YData'), [0,0,1], 'FaceAlpha', 0.5);
    end
end

%% TROP2 IQR

sumTemp = varfun(@std,alltumor,'GroupingVariables','slideName');
sumTemp = join(sumTemp,slideINFO,'keys','slideName');
figure,myboxplot2dj(log(sumTemp.std_TROP2),sumTemp.Response);
set(gca,'xticklabels',{'pCR','RD'});
title('TROP2 std','FontSize',16);

%% TROP2 std (Fig 2X)
sumTemp = varfun(@std,alltumor,'GroupingVariables','slideName');
sumTemp = join(sumTemp,slideINFO,'keys','slideName');
figure,myboxplot2d(sumTemp.std_TROP2,sumTemp.Response);
set(gca,'xticklabels',{'pCR','RD'});
title('TROP2 Std.','FontSize',16);


%% forest (topics)

sum2 = sumAll;
%title1 = 'Whole sample';
odd_ratios = zeros(maxT,1);
CI = zeros(maxT,2);
labels = cell(maxT,1)
for i = 1:maxT
    marker1 = strcat('mean_topic',num2str(i));
    labels{i} = strcat('topic',num2str(i));
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
end
myforestplot(labels,odd_ratios,CI);
 
xlabel('fold difference');
title('RD/pCR (whole sample)','FontSize',18);

xscale('log');
ylim([0 17]);

%% Gating TROP2 automatically

gateTROP2k = zeros(length(slideName),1);
k = 2;

tic; 
for s = 1:length(slideName)
%data1 = dataLSP20627;
    disp(strcat('Processsing:',slideName{s},'.....',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    [idx,c]=kmeans(log(data1.TROP2+100),k,'MaxIter',100,'Replicate',5);
    idxnew = zeros(length(idx),1);
    [~, order] = sort(c,'ascend');
    for i = 1:k
        idxnew(idx==order(i))=i;
    end
    data1.TROP2pk = idxnew ==k;
    gateTROP2k(s)= log(min(data1.TROP2(data1.TROP2pk)));
    
    data1.TROP2hk = data1.TROP2 > exp(gateTROP2k(s)+1);
    data1.TROP2lk = data1.TROP2pk & ~data1.TROP2hk;
    data1.TROP2ck = zeros(size(data1,1),1);
    data1.TROP2ck(data1.E_cadherinp)=1;
    data1.TROP2ck(data1.TROP2lk)=2;
    data1.TROP2ck(data1.TROP2hk)=3;
    tabulate(data1.TROP2ck);
    eval(strcat('data',gateTable.slideName{s},'=data1;'));
    
    toc;
end


%% test

k = 2;
data1 = dataLSP22336;

[idx,c]=kmeans(log(data1.TROP2+5),k,'MaxIter',200,'Replicate',5);
idxnew = zeros(length(idx),1);
[~, order] = sort(c,'ascend');
for i = 1:k
    idxnew(idx==order(i))=i;
end
data1.TROP2pk = idxnew ==k;


%% Neighborhood analysis (Trop2 Ki67+/-)

arrayNeighborTROP2pp = zeros(length(labelp),length(slideName));
arrayNeighborTROP2np = zeros(length(labelp),length(slideName));

%cutoff = 20;

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data0 = eval(strcat('data',name1));
    data2 = data0(data0.TROP2p & data0.Ki_67p,:);
    data3 = data0(data0.TROP2p & ~data0.Ki_67p,:);

    for i = 1:length(labelp2)
        data1 = data0(data0{:,labelp2{i}},:);
        data4 = datasample(data0,size(data1,1),'replace',false);
        [~,d12]= knnsearch([data1.Xt,data1.Yt],[data2.Xt,data2.Yt],'k',2);
        [~,d42]= knnsearch([data4.Xt,data4.Yt],[data2.Xt,data2.Yt],'k',2);
        arrayNeighborTROP2pp(i,s)= mean(d42(:,2)) ./ mean(d12(:,2));

        [~,d13]= knnsearch([data1.Xt,data1.Yt],[data3.Xt,data3.Yt],'k',2);
        [~,d43]= knnsearch([data4.Xt,data4.Yt],[data3.Xt,data3.Yt],'k',2);
        arrayNeighborTROP2np(i,s)= mean(d43(:,2)) ./ mean(d13(:,2));
    end
    toc;
end

%% Boxplot (TROP2 Ki67+/- Neighbor)

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayNeighborTROP2pp(i,:),sumAll.Response);
    title(labelp3{i},'Interpreter','none');
    set(gca,'xticklabels',{'pCR','RD'})
end

set(gcf,'color','w');
sgtitle('To TROP2+Ki67+','FontSize',14);

%% Distance between Immune cells and TROP2 high/low

arrayDisTROP2pp = zeros(length(labelp),length(slideName));
arrayDisTROP2np = zeros(length(labelp),length(slideName));

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data1 = data1(data1.ROI==0,:);
    data2 = data1(data1.TROP2p & data1.Ki_67p,:);
    data3 = data1(data1.TROP2p & ~data1.Ki_67p,:);
    for i = 1:length(labelp2)
        data4 = data1(data1{:,labelp2{i}},:);
        data5 = datasample(data1,size(data4,1));
        
        [~,d1] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'k',1);
        [~,d2] = knnsearch([data2.Xt,data2.Yt],[data5.Xt,data5.Yt],'k',1);
        arrayDisTROP2pp(i,s)= mean(d2)/mean(d1);

        [~,d1] = knnsearch([data3.Xt,data3.Yt],[data4.Xt,data4.Yt],'k',1);
        [~,d2] = knnsearch([data3.Xt,data3.Yt],[data5.Xt,data5.Yt],'k',1);
        arrayDisTROP2np(i,s)= mean(d2)/mean(d1);
    end
    toc;
end

%% Boxplot (arrayDisTROP2pp)

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(arrayDisTROP2np(i,:),sumAll.Response);
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle('To TROP2+Ki67-');

%% Paired boxplot

sum2 = sumAll;
groupc = repmat({'b'},size(sum2,1),1);

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2p([arrayNeighborTROP2np(i,:)',arrayNeighborTROP2pp(i,:)'],groupc);
    set(gca,'xticklabels',{'Ki67-','Ki67+'});
    ylabel('Interacting Strength');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');

%% Distance 

arrayDisTROP2np2VES = zeros(length(slideName),1);

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data2 = data1(data1.TROP2p & ~data1.Ki_67p,:);
  
    data4 = data1(data1.VesID>0,:);
    data5 = datasample(data1,size(data4,1));
        
    [~,d1] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'k',1);
    [~,d2] = knnsearch([data2.Xt,data2.Yt],[data5.Xt,data5.Yt],'k',1);
    arrayDisTROP2np2VES(s)= mean(d2)/mean(d1);
end

figure,myboxplot2d(arrayDisTROP2pp2VES,sumAll.Response);

%% Distance 

arrayDisTROP2np2CD31nv = zeros(length(slideName),1);

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));
    data2 = data1(data1.TROP2p & data1.Ki_67p,:);
  
    data4 = data1(data1.CD31p & data1.VesID==0,:);
    data5 = datasample(data1,size(data4,1));
        
    [~,d1] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'k',1);
    [~,d2] = knnsearch([data2.Xt,data2.Yt],[data5.Xt,data5.Yt],'k',1);
    arrayDisTROP2np2CD31nv(s)= mean(d2)/mean(d1);
end

figure,myboxplot2d(arrayDisTROP2pp2CD31nv,sumAll.Response);

%% Visualization

temp1 = allint1r ./ allint1;

allintMeanR = mean(temp1(:,:,ismember(sumAll.Response,'R')),3);
allintMeanNR = mean(temp1(:,:,ismember(sumAll.Response,'NR')),3);

figure,mydotplot(allIntPvalue,allIntRatio,labelp3,labelp3);

%% Vessel counts 

% CycIF_Generalsegmentation_20240804

%alldata.CD31nv = alldata.CD31p & alldata.VesID==0;

sumVES = varfun(@mean,alldata,'GroupingVariables',{'slideName','VesID'});
sumVES = varfun(@mean,sumVES,'GroupingVariables','slideName');
sumVES = join(sumVES,slideINFO,'keys','slideName');
figure,myboxplot2d(sumVES.GroupCount-1,sumVES.Response);
title('Vessel counts','FontSize',16);
set(gca,'xticklabels',{'pCR','RD'});
ylabel('number of Vessels');

%% Calculate Distance to nearest vessels & TROP2+ cells

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'.....',num2str(s)))
    data1 = eval(strcat('data',slideName{s}));
    data1.CD31nv = data1.CD31p & data1.VesID==0;
    data1.CD31v = data1.CD31p & data1.VesID>0;

    data2 = data1(data1.VesID>0,:);
    [~,d]=knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',1);
    data1.DisVes = d;

    data2 = data1(data1.TROP2p,:);
    [~,d]=knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',1);
    data1.DisTROP2 = d;
    eval(strcat('data',slideName{s},'=data1;'));
    toc;
end

%% Calculate allDisVes

allDisVes = zeros(length(slideName),length(labelp));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'.....',num2str(s)))
    data1 = eval(strcat('data',slideName{s}));
    for i = 1:length(labelp)
        allDisVes(s,i)=mean(data1.DisVes(data1{:,labelp2{i}}));
    end
    toc;
end

%% Calculate ROC curve (CD31+Ki67)

list1 = sumTumor.mean_Ki_67p;
list2 = sumAll.mean_CD31p;

xd = (max(list1)-min(list1))/50;
x = min(list1)-xd:xd:max(list1)+xd;

yd = (max(list2)-min(list2))/50;
y = min(list2)-yd:yd:max(list2)+yd;

tpr_x = zeros(length(x),1);
fpr_x = zeros(length(x),1);
tpr_y = zeros(length(y),1);
fpr_y = zeros(length(y),1);
tpr_xy = zeros(length(x),length(y));
fpr_xy = zeros(length(x),length(y));

for i = 1:length(x)
    [tpr_x(i),fpr_x(i)]=myTPR_FPR(ismember(sumTumor.Response,'R'),list1>=x(i));
end
auc_x = abs(trapz(fpr_x,tpr_x));

for j = 1:length(y)
    [tpr_y(j),fpr_y(j)]=myTPR_FPR(ismember(sumAll.Response,'R'),list2<=y(j));
end
auc_y = abs(trapz(fpr_y,tpr_y));

for i = 1:length(x)
    for j=1:length(y)
        [tpr_xy(i,j),fpr_xy(i,j)]=myTPR_FPR(ismember(sumAll.Response,'R'),list1>=x(i) & list2<=y(j));
    end
end
fpr_xy = fpr_xy(:);
tpr_xy = tpr_xy(:);
tpr_xy = groupsummary(tpr_xy,fpr_xy,'max');
fpr_xy = unique(fpr_xy);

auc_xy = abs(trapz(fpr_xy,tpr_xy));

figure,plot(fpr_x(:),tpr_x(:),'r-','LineWidth',1);
hold on;
plot(fpr_y(:),tpr_y(:),'g-','LineWidth',1);
plot(fpr_xy(:),tpr_xy(:),'b-','LineWidth',1);

% p = polyfit(fpr,tpr,3);
% y_fit = polyval(p,fpr);
% plot(fpr, y_fit, 'r-', 'LineWidth', 2);
hold on;
plot([0,1],[0,1],'k--','LineWidth',0.5);
xlabel('FPR');
ylabel('TPR');
legend(strcat('Ki67(auc=',num2str(auc_x),')'),strcat('CD31(auc=',num2str(auc_y),')'),strcat('CD31+Ki67(auc=',num2str(auc_xy),')'),'Location','best');
grid on;


%% Calculate ROC curve (CD31 only

list2 = sumAll.mean_topic3;

yd = (max(list2)-min(list2))/50;
y = min(list2)-yd:yd:max(list2)+yd;

tpr_y = zeros(length(y),1);
fpr_y = zeros(length(y),1);

for j = 1:length(y)
    [tpr_y(j),fpr_y(j)]=myTPR_FPR(ismember(sumAll.Response,'R'),list2>=y(j));
end


figure,
plot(fpr_y(:),tpr_y(:),'g-','LineWidth',1);
auc_y = abs(trapz(fpr_y,tpr_y));

% p = polyfit(fpr,tpr,3);
% y_fit = polyval(p,fpr);
% plot(fpr, y_fit, 'r-', 'LineWidth', 2);
hold on;
plot([0,1],[0,1],'k--','LineWidth',0.5);
xlabel('FPR');
ylabel('TPR');
legend('Topic 9','Location','best');
grid on;
title(strcat('auc=',num2str(auc_y)),'FontSize',16);

%% 

arrayDP = zeros(length(labelp),length(labelp),length(slideName));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s}));
    data1 = eval(strcat('data',slideName{s}));
    for i = 1:length(labelp)
        for j = i+1:length(labelp)
            arrayDP(i,j,s)= mean(data1{:,labelp2{i}} & data1{:,labelp2{j}})/(mean(data1{:,labelp2{i}})*mean(data1{:,labelp2{j}}));
        end
    end
    toc;
end

%% Cell type clustering (Kmeans)

% CycIF_alldata;

%idx = kmeans(alldata{:,labelp2},7,'MaxIter',200,'Display','final','Replicates',5);
%alldata.kmean7 = idx;

temp1 = varfun(@mean,alldata(:,[1:98,100:end]),'GroupingVariables','kmean11');
figure,imagesc(temp1{:,strcat('mean_',labelp2)});
set(gca,'xtick',1:length(labelp2));
set(gca,'xticklabels',labelp3);

colormap(jet);

%% Reassign cell tyeps (kmean11)

% 1 = 1 tumor 1 (Ki67+)
% 4 = 2 tumor 2 (Ki67-)
% 8 = 3 tumor 3 (TROP2-)
% 11 = 4 tumor 4 (Mix)
% 2 = 5 Immune 1 (CD4 T)
% 5 = 6 Immune 2 (CD8 T)
% 6 = 7 Immune 3 (Myeloid)
% 9 = 8 Immune 4 (CD163+)
% 10 = 9 Immune 5 (B cells)
% 7 = 10 Stroma 1 (CD31+)
% 3 = 11 Stroma 2 (Other)

alldata.CellType = zeros(size(alldata,1),1);

Typelist = [1 5 11 2 6 7 10 3 8 9 4];

for i = 1:11
    alldata.CellType(alldata.kmean11==i)= Typelist(i);
end

temp1 = varfun(@mean,alldata(:,[1:98,100:end]),'GroupingVariables','CellType');
figure,imagesc(temp1{:,strcat('mean_',labelp2)});
set(gca,'xtick',1:length(labelp2));
set(gca,'xticklabels',labelp3);
colormap(jet);

% ----- Assign Cell Category (1:tumor 2:Immune 3: Stroma)

alldata.CellCat = zeros(size(alldata,1),1);
Catlist = [1 1 1 1 2 2 2 2 2 3 3];

for i = 1:11
    alldata.CellCat(alldata.CellType==i)= Catlist(i);
end

temp1 = varfun(@mean,alldata(:,[1:98,100:end]),'GroupingVariables','CellType');
figure,imagesc(temp1{:,strcat('mean_',labelp2)});
set(gca,'xtick',1:length(labelp2));
set(gca,'xticklabels',labelp3);
colormap(jet);

%% Assign to Celltype all slides

labelCT = {'Tumor1','Tumor2','Tumor3','Tumor4','Immune1','Immune2','Immune3','Immune4','Immune5','Stroma1','Stroma2'};
labelCTname = {'Tumor1(Ki67+)','Tumor2(Ki67-)','Tumor3(TROP2-)','Tumor4(MIX)','Immune1(CD4 T)','Immune2(CD8 T)','Immune3(Myeloid)','Immune4(CD163+)','Immune5(B cells)','Stroma1(CD31+)','Stroma2(Other)'};

labelCC = {'Tumor','Immune','Stroma'};

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'...',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    data2 = alldata(ismember(alldata.slideName,slideName{s}),:);
    data1.CellType = data2.CellType;
    data1.CellCat = data2.CellCat;
    for i = 1:length(labelCT)
        data1{:,labelCT{i}} = data1.CellType ==i;
    end
    for j = 1:length(labelCC)
        data1{:,labelCC{j}} = data1.CellCat == j;
    end
    eval(strcat('data',slideName{s},'=data1;'));
    toc;
end

%% Plot cell type heatmap

% CycIF_alldata;

temp1 = varfun(@mean,alldata(:,1:end-1),'GroupingVariables','CellType');  % remove slideName
figure,imagesc(zscore(temp1{:,strcat('mean_',labelp2)}));
colormap(jet);

set(gca,'xtick',1:length(labelp2));
set(gca,'xticklabels',labelp3);
set(gca,'yticklabels',labelCTname);
set(gca,'ytick',1:length(labelCT));


%% Stacked bar plot (cell types)

%CycIF_alldata;
temp1 = varfun(@mean,alldata,'GroupingVariables','slideName');

figure('units','normalized','outerposition',[0 0 1 1]);

colors = jet(11);
% ---- legend----
subplot(1,20,1);
h1=bar(1,ones(11,1), 'stacked');
for k = 1:numel(h)
    h1(k).FaceColor = colors(k,:);
    %h1(k).EdgeColor = [1 1 1];
end
xlim([0.2 1.8]);
ylim([0 11]);
set(gca,'ytick',0.5:1:10.5);
set(gca,'yticklabels',labelCT);
%set(gca,'ytick',[]);
set(gca,'XColor','none');
%set(gca,'YColor','none');

% ---- Actual plot---
subplot(1,20,2:20);

h2= bar(temp1{:,strcat('mean_',labelCT)},'stacked');
ylim([0 1]);
colors = jet(11);

for k = 1:numel(h)
    h2(k).FaceColor = colors(k,:);
    %h2(k).EdgeColor = [1 1 1];
end
ylim([0 1]);
set(gca,'ytick',[]);
xlim([0.5 48.5]);
set(gca,'xtick',[]);
xlabel('Samples');

%% forest (cell types)

%temp1 = varfun(@mean,alldata,'GroupingVariables','slideName');
temp2 = varfun(@mean,alldata(alldata.Region>0,:),'GroupingVariables','slideName');

sum2 = temp2;
sum2 = join(sum2,slideINFO,'keys','slideName');

odd_ratios = zeros(length(labelCT),1);
CI = zeros(length(labelCT),2);
%labels = cell(length(labelCT),1);

for i = 1:length(labelCT)
    marker1 = strcat('mean_',labelCT{i});
    %labels{i} = strcat('topic',num2str(i));
    [odd_ratios(i),CI(i,1),CI(i,2)]=mygroupfold(flip(sum2{:,marker1}*100),flip(sum2.Response),false);
end
myforestplot(labelCTname,odd_ratios,CI);
 
xlabel('fold difference');
title('RD/pCR (Tumor region)','FontSize',18);

%xscale('log');
ylim([0 length(labelCT)+1]);

%% Triparticle analysis (all combs)

n = 11;  % range of numbers
r = 3;   % number to choose

% Generate combinations with repetition
combs = nchoosek(1:n + r - 1, r);   % choose positions
combs = combs - repmat(0:r-1, size(combs,1), 1);  % adjust values
combsname = cell(size(combs,1),1);

parfor i = 1:size(combs,1)
    combsname{i} = strjoin(string(combs(i,:)),':');
end
combsname = cellstr(combsname);

%% Triparticle analysis (sample list)

allTripart = zeros(length(combsname),length(slideName));

tic;
for s = 1:length(slideName)
    %data1 = dataLSP20621;
    data1 = eval(strcat('data',slideName{s}));
    disp(strcat('Processing:',slideName{s},'...',num2str(s)));

    idx = knnsearch([data1.Xt,data1.Yt],[data1.Xt,data1.Yt],'k',3);
    idx = data1.CellType(idx);
    idx = sort(idx,2);
    idx = sort(idx,1);
    
    temp1 = cell(size(idx,1),1);
    parfor i = 1:size(idx,1)
        temp1{i} = strjoin(string(idx(i,:)),':');
    end
    temp1 = cellstr(temp1);
    temp2 = tabulate(temp1);
    
    for i = 1:size(temp2,1)
        combname = find(ismember(combsname,temp2{i,1}));
        allTripart(combname,s)=temp2{i,3};
    end
    toc;
end

%% Triparticle analysis (P values)

allTripPvalue = ones(length(combsname),1);

for i = 1:length(combsname)
    if sum(allTripart(i,:))>5
        allTripPvalue(i)= myttest2(allTripart(i,:),sumAll.Response);
    end
end

%% Plot triparticle interaction
i = 67;

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);

% ---- Plot triangle ----
subplot(1,2,1)
x = [0, 1, 0.5];
y = [0, 0, sqrt(3)/2];
plot([x x(1)], [y y(1)], 'k-', 'LineWidth', 1.5); % Close the triangle by repeating the first point
axis equal
hold on

labels = labelCTname(combs(i,:));
text(0,0,labels{1},'FontSize', 12, 'FontWeight', 'bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(1,0,labels{2},'FontSize', 12, 'FontWeight', 'bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(0.5,0.866,labels{3},'FontSize', 12, 'FontWeight', 'bold','HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

hold off
box off;
axis off;
xlim([-0.25 1.25]);
ylim([-0.25 1.25]);

% ---- Boxplot ----
subplot(1,2,2);
myboxplot2d(allTripart(i,:),sumAll.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});


%% calculate all cluster score  (Cell Type)

arrayClusterScoreCT = zeros(length(labelCT),length(labelCT),length(slideName));

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    for i = 1:length(labelCT)
        gate1 = labelCT{i};
        for j = 1:length(labelCT)
            gate2 = labelCT{j};
           
            %gate1 = strcat(marker1,'p');
            %gate2 = strcat(marker2,'p');
            
            data2 = data1(data1{:,gate1},:);
            data3 = data1(data1{:,gate2},:);
            data4 = datasample(data1,size(data3,1));
            [~,d1] = knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],'K',2);
            [~,d2] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'K',2);
            arrayClusterScoreCT(i,j,s) = mean(d2(:,2))/mean(d1(:,2));
        end
    end
    toc;
end

%% Test all pvalue (cluster Cell Type)

arrayClusterCTPvalues = ones(length(labelCT),length(labelCT));
arrayClusterCTfolds = zeros(length(labelCT), length(labelCT));
tic;
for i = 1:length(labelCT)
    for j = 1:length(labelCT)
        [arrayClusterCTPvalues(i,j),arrayClusterCTfolds(i,j)]=myttest2(squeeze(arrayClusterScoreCT(i,j,:)),sumAll_new.Response);
    end
end
toc;



%% Plot interaction (Cell type)

Xpos = repmat((1:length(labelCT))',1,length(labelCT));
Ypos = repmat(1:length(labelCT),length(labelCT),1);

figure,scatter(Xpos(:),Ypos(:),-log(arrayClusterCTPvalues(:))*15,arrayClusterCTfolds(:),'fill');
xlim([0.5 length(labelCT)+0.5]);
ylim([0.5 length(labelCT)+0.5]);

set(gca,'xtick',1:length(labelCT));
set(gca,'xticklabels',labelCTname);
set(gca,'ytick',1:length(labelCT));
set(gca,'yticklabels',labelCTname);
colormap(redbluecmap);
colorbar;
caxis([0 2]);

%% plot single interaction (Cell type)

X = 6;
Y = 9;

figure,myboxplot2d(squeeze(arrayClusterScoreCT(X,Y,:)),sumAll.Response);
title(strcat(labelCTname{X},'::',labelCTname{Y}),'FontSize',16);
xlims = xlim;
line(xlims,[1 1],'LineWidth',0.5,'LineStyle','--');
set(gca,'xticklabels',{'pCR','RD'});

%% Ripley's K/G function (TROP2+Ki67+ or TROP2+Ki67-)

listR = 0:1:25;

arrayG_TROP2Ki67 = zeros(length(slideName),length(listR));
%arrayK_TROP2Ki67 = zeros(length(slideName),length(listR));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'....',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    %data1.TROP2pKi_67n = data1.TROP2p & ~data1.Ki_67p;
    %eval(strcat('data',slideName{s},'=data1;'));
    data2 = data1(data1.TROP2pKi_67p,:);
    data3 = data1(data1.TROP2pKi_67n,:);
    parfor i = 1:length(listR)
        arrayG_TROP2Ki67(s,i)=bivariateG([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],listR(i));
        %arrayK_TROP2Ki67(s,i)=bivariateRipleyK([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],listR(i),arrayArea(s)*10000);
    end
    toc;
end

%% check 

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(listR_old)
    subplot(3,9,i);
    myboxplot2d(arrayG_TROP2LH(:,i),sumAll.Response);
    title(strcat('r=',num2str(listR(i))));
end
set(gcf,'color','w');

%% Ripley's G function (TROP2lk or TROP2hk)

listR = 0:1:100;

arrayG_TROP2LH = zeros(length(slideName),length(listR));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'....',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    data2 = data1(data1.TROP2lk,:);
    data3 = data1(data1.TROP2hk,:);
    parfor i = 1:length(listR)
        arrayG_TROP2LH(s,i)=bivariateG([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],listR(i));
    end
    toc;
end

%% check 

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(listR)
    subplot(3,9,i);
    myboxplot2d(arrayG_TROP2Ki67(:,i),sumAll.Response);
    title(strcat('r=',num2str(listR(i))));
end
set(gcf,'color','w');

%% Ripley's K function (TROP2lk or TROP2hk)

listR = 0:0.4:10;
arrayK_TROP2LH = zeros(length(slideName),length(listR));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'....',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    data2 = data1(data1.TROP2lk,:);
    data3 = data1(data1.TROP2hk,:);
    parfor i = 1:length(listR)
        arrayK_TROP2LH(s,i)=bivariateRipleyK([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],listR(i),arrayArea(s)*10000);
    end
    toc;
end

%% Ripley's K function (TROP2Ki67+ or TROP2Ki67-)

listR = 0:5:100;
arrayK_TROP2Ki67 = zeros(length(slideName),length(listR));

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'....',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    data2 = data1(data1.TROP2pKi_67p,:);
    data3 = data1(data1.TROP2p & ~data1.Ki_67p,:);
    parfor i = 1:length(listR)
        arrayK_TROP2Ki67(s,i)=bivariateRipleyK([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],listR(i),arrayArea(s)*10000);
    end
    toc;
end

%%   
figure;
for i = 1:length(listR)
    subplot(4,13,i);
    myboxplot2d(arrayK_TROP2LH(:,i),sumAll.Response);
end

%% Sandro's plot

figure('units','normalized','outerposition',[0.5 0.5 0.25 0.5]);
myboxplot2dj(sumTumor.mean_CD31p*100,sumTumor.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);
title('CD31+ cells','FontSize',16);
set(gcf,'color','w');


%%  TRO2/Ki67

sum1 = sumAll;

figure('units','normalized','outerposition',[0.5 0.5 0.25 0.5]);

myboxplot2d(sum1.mean_TROP2p./sum1.mean_Ki_67p,sum1.Response);
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);
title('TROP2+/Ki67+ ratio','FontSize',16);
set(gcf,'color','w');

%% Test all pvalue (cluster score)

temp1 = allint1r ./ allint1;
temp1 = temp1';

arrayClusterPvalues = ones(length(labelp),length(labelp));
arrayClusterFolds = zeros(length(labelp), length(labelp));
tic;
for i = 1:length(labelp)
    for j = 1:length(labelp)
        [arrayClusterPvalues(i,j),arrayClusterFolds(i,j)]=myttest2(squeeze(temp1(i,j,:)),sumAll.Response);
    end
end
toc;

%% Plot interaction (Marker)  (Fig4A)

Xpos = repmat((1:length(labelp))',1,length(labelp));
Ypos = repmat(1:length(labelp),length(labelp),1);

figure,scatter(Ypos(:),Xpos(:),-log(arrayClusterPvalues(:))*5,arrayClusterFolds(:),'fill');
xlim([0.5 length(labelp)+2.5]);
ylim([0.5 length(labelp)+0.5]);
hold on;
scatter(length(labelp)+1,2,-log(0.01)*5,'k');
text(length(labelp)+0.5,1,'p=0.05','FontSize',6);
scatter(length(labelp)+1,4,-log(0.001)*5,'k');
text(length(labelp)+0.5,3,'p=0.01','FontSize',6);
scatter(length(labelp)+1,6,-log(0.0001)*5,'k');
text(length(labelp)+0.5,5,'p=0.001','FontSize',6);

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
colormap(redbluecmap);
colorbar;
caxis([0 2]);

%% Selected inteaciton  single Boxplot (Second method)

i = find(ismember(labelp,'CD8a'));
j = find(ismember(labelp,'CD68'));

temp1 = allint1r ./ allint1;

figure('units','normalized','outerposition',[0.5 0.5 0.25 0.5]);

myboxplot2dj(squeeze(temp1(i,j,:)),sumAll.Response);
title(strcat(labelp{i},'-',labelp{j},' interaction'),'FontSize',16);  
ylabel('Normalized Interaction','FontSize',12);
set(gcf,'color','w');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

%% Plot curve for arrryG

arrayNR = arrayG_TROP2Ki67(ismember(sumAll.Response,'NR'),:);
arrayR = arrayG_TROP2Ki67(ismember(sumAll.Response,'R'),:);

figure,plot(arrayNR','Color','b','LineWidth',0.75);
hold on;
plot(arrayR','Color','r','LineWidth',1);
%legend('pCR','RD');
xlim([0 26]);
ylims = ylim;
line([14 14],ylims,'LineStyle','--','Color','k');

line([1 3],[0.7 0.7],'Color','r');
text(3,0.7,'pCR');
line([1 3],[0.6 0.6],'Color','b');
text(3,0.6,'RD');

ylabel('Ripley G score');
xlabel('Distance(um)');
title('TROP2+Ki67+ versus TROP2+Ki67-','FontSize',16);

%% 
i = 14;
figure,myboxplot2d(arrayG_TROP2Ki67(:,i),sumAll.Response);
ylabel('Ripley G score');
set(gca,'xticklabels',{'pCR','RD'});

%% Calculate TROP2 curve for all tumors

xi = zeros(100,length(slideName));
f = zeros(100,length(slideName));

for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'.....',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    data2 = data1(data1.E_cadherinp,:);
    [xi(:,s),f(:,s)]=ksdensity(log(data2.TROP2));
end

xi_r = xi(:,ismember(sumAll.Response,'R'));
f_r = f(:,ismember(sumAll.Response,'R'));
figure,plot(f_r,xi_r,'LineWidth',1.5,'Color','r');
hold on;

xi_nr = xi(:,ismember(sumAll.Response,'NR'));
f_nr = f(:,ismember(sumAll.Response,'NR'));
plot(f_nr,xi_nr,'LineWidth',1,'Color','b');
hold off;

xlim([5 12]);

%% Display LDA topics (Fig4F)

figure
for topicIdx = 1:maxT
    subplot(4,ceil(maxT/4),topicIdx)
    temp1 = table;
    temp1.Word = labelp3;
    
    temp1.Count = ldaAll.TopicWordProbabilities(:,topicIdx);
    wordcloud(temp1,'Word','Count');
    title("Topic: " + topicIdx)
end
toc;

%% Disorganized CD31+ cells (Fig5C_right)

figure,myboxplot2dj(sumAll_2025.mean_CD31nv*100,sumAll_2025.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});

%% Total CD31+ cells (Fig 5A)

figure,myboxplot2dj(sumAll_2025.mean_CD31p*100,sumAll_2025.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});

%% Vessel counts (Fig5C_left)

figure,myboxplot2dj(sumVES.GroupCount,sumVES.Response);
ylabel('Counts');
set(gca,'xticklabels',{'pCR','RD'});

%% TROP2 size & count (Fig2E)

% sumA = varfun(@mean,alldata,'GroupingVariables',{'slideName','TROP2pID'});
% sumTROP2pID = varfun(@mean,sumA,'GroupingVariables','slideName');
% sumTROP2pID = join(sumTROP2pID,slideINFO,'keys','slideName');

figure,myboxplot2dj(sumTROP2pID.GroupCount,sumTROP2pID.Response);
title('TROP2 cluster counts');
set(gca,'xticklabels',{'pCR','RD'});

figure,myboxplot2dj(sumTROP2pID.mean_GroupCount,sumTROP2pID.Response);
title('TROP2 cluster size (by cells)');
ylabel('Cell counts');
set(gca,'xticklabels',{'pCR','RD'});

%% Replot cluster size with color scale (Fig2K)

% P09 = LSP20623
% P37 = LSP22324

data1 = dataLSP20623;
sum1 = varfun(@mean,data1,'GroupingVariables','TROP2pID');
sum1 = sum1(:,1:2);
data1 = join(data1,sum1,'keys','TROP2pID');
data1.GroupCount(data1.TROP2pID==0) = 0;
%  --- Plot ---
figure,scatter(data1.Xt,data1.Yt,3,repmat([0.7,0.7,0.7],size(data1,1),1),'fill');
hold on;
data2 = data1(data1.GroupCount>0,:);
CycIF_tumorview(data2,'GroupCount',7,0);
daspect([1 1 1]);
colormap(cool);
caxis([1 2000]);
title('P37 (RD)','FontSize',14);

%% Plot relationship between TROP2 cluster size and Ki67+ (Fig 2XX)
colors = repmat([1,0,0],size(sumAll,1),1);

for i = 1:size(sumAll,1)
    if ismember(sumAll.Response(i),'NR')
        colors(i,:) = [0,0,1];
    end
end

figure,scatter(sumAll.mean_Ki_67p*100,sumTROP2pID.mean_GroupCount,40,colors,'fill');
xtickformat('percentage');
xlabel('Ki67+');
ylabel('TROP2+ cluster size');

%% Ki67 in pCR/RD (Fig 2F)

figure,myboxplot2dj(sumTumor.mean_Ki_67p*100,sumTumor.Response);
ytickformat('percentage');
set(gca,'xticklabels',{'pCR','RD'});
ylabel('%Ki67+ cells');
title('Ki67+ Tumor cells','FontSize',16);

%% Ki67 in four groups (label)
sum1 = sumTumor;
sum2 = sum1(sum1.Recurr,:);
figure,myboxplot4d(sum1.mean_Ki_67p*100,sum1.GroupOrder);
hold on;
scatter(sum2.GroupOrder,sum2.mean_Ki_67p*100,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.mean_Ki_67p*100;
txt = sum2.PatientID;

for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

ytickformat('percentage');
set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('%Ki67+ cells');
title('Ki67+ Tumor cells','FontSize',16);

%% Ki67 in Three groups
sum1 = sumTumor;

sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;
figure,myboxplot3d(sum1.mean_Ki_67p*100,sum1.Group3,10);
ytickformat('percentage');
set(gca,'xticklabels',listGroup3);
ylabel('%Ki67+ cells');
title('Ki67+ Tumor cells','FontSize',16);

%% CD31 in Three groups
sum1 = sumAll_2025;

sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;
figure,myboxplot3d(sum1.mean_CD31p*100,sum1.Group3,12);
ytickformat('percentage');
set(gca,'xticklabels',listGroup3);
ylabel('% positive cells');
title('CD31+ cells','FontSize',16);

%% CD31 in four groups
sum1 = sumAll_2025;
sum2 = sum1(sum1.Recurr,:);

figure,myboxplot4d(sum1.mean_CD31p*100,sum1.GroupOrder);
hold on;
scatter(sum2.GroupOrder,sum2.mean_CD31p*100,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.mean_CD31p*100;
txt = sum2.PatientID;

for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

ytickformat('percentage');
set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('% postive cells');
title('CD31+ cells','FontSize',14);

%% Ki67/CD8+ in four groups

sum1 = sumTumor;

figure,myboxplot4d(sum1.mean_CD8apKi_67p./sum1.mean_CD8ap*100,sum1.GroupOrder);

% --- label recurr ---
sum2 = sum1(sum1.Recurr,:);
hold on;
scatter(sum2.GroupOrder,sum2.mean_CD8apKi_67p./sum2.mean_CD8ap*100,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.mean_CD8apKi_67p./sum2.mean_CD8ap*100;
txt = sum2.PatientID;
for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

ytickformat('percentage');
set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('%Ki67+ cells');
title('Proliferative CD8+ cells','FontSize',14);

%% Ki67/CD8+ in three groups

sum1 = sumTumor;
sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;

figure,myboxplot3d(sum1.mean_CD8apKi_67p./sum1.mean_CD8ap*100,sum1.Group3);
ytickformat('percentage');
set(gca,'xticklabels',listGroup3);
ylabel('%Ki67+ cells');
title('Proliferative CD8+ cells','FontSize',14);



%% Vessel count in four group

sum1 = sumVES;
figure,myboxplot4d(sum1.GroupCount,sum1.GroupOrder);

% --- label recurr ---
sum2 = sum1(sum1.Recurr,:);
hold on;
scatter(sum2.GroupOrder,sum2.GroupCount,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.GroupCount;
txt = sum2.PatientID;
for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

ylabel('Counts');
set(gca,'xticklabels',tableGroupOrder.Group);
title('Vessel count','FontSize',14);


%% CD8a/CD4 in four groups

sum1 = sumAll_2025;
figure,myboxplot4d(sum1.mean_CD8ap./sum1.mean_CD4p,sum1.GroupOrder);

% --- label recurr ---
sum2 = sum1(sum1.Recurr,:);
hold on;
scatter(sum2.GroupOrder,sum2.mean_CD8ap./sum2.mean_CD4p,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.mean_CD8ap./sum2.mean_CD4p;
txt = sum2.PatientID;
for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('Ratio');
title('CD8/CD4 ratio','FontSize',14);

%% CD8a/CD4 in three groups

sum1 = sumAll_2025;
sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;

figure,myboxplot3d(sum1.mean_CD8ap./sum1.mean_CD4p,sum1.Group3);

set(gca,'xticklabels',listGroup3);
ylabel('Ratio');
title('CD8/CD4 ratio','FontSize',14);

%% Topics in four group

sum1 = sumAll;
figure,myboxplot4d(sum1.mean_topic10*100,sum1.GroupOrder);

% --- label recurr ---
sum2 = sum1(sum1.Recurr,:);
hold on;
scatter(sum2.GroupOrder,sum2.mean_topic10*100,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.mean_topic10*100;
txt = sum2.PatientID;
for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

ytickformat('percentage');
set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('% positive cells');
title('Topic 10','FontSize',14);


%% Topics in three group

sum1 = sumAll;
sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;

figure,myboxplot3d(sum1.mean_topic10*100,sum1.Group3);
ytickformat('percentage');
set(gca,'xticklabels',listGroup3);
ylabel('% positive cells');
title('Topic 10','FontSize',14);

%% TROP2 cluster size in four group
sum1 = sumTROP2pID;

figure,myboxplot4d(sum1.mean_GroupCount,sum1.GroupOrder);

% --- label recurr ---
sum2 = sum1(sum1.Recurr,:);
hold on;
scatter(sum2.GroupOrder,sum2.mean_GroupCount,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.mean_GroupCount;
txt = sum2.PatientID;
for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

title('TROP2 cluster size','FontSize',14);
ylabel('Cell counts');
set(gca,'xticklabels',tableGroupOrder.Group);

%% TROP2 cluster size in three group
sum1 = sumTROP2pID;
sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;

figure,myboxplot3d(sum1.mean_GroupCount,sum1.Group3);
title('TROP2 cluster counts','FontSize',14);
ylabel('Cell counts');
set(gca,'xticklabels',listGroup3);

%% TROP2 cluster count in four group
sum1 = sumTROP2pID;

figure,myboxplot4d(sum1.GroupCount,sum1.GroupOrder);
title('TROP2 cluster size (by cells)','FontSize',14);
ylabel('Cluster counts');
set(gca,'xticklabels',tableGroupOrder.Group);

%% TROP2 cluster size in three group
sum1 = sumTROP2pID;

sum1.Group3 = sum1.GroupOrder;
sum1.Group3(sum1.Group3==4)=3;

figure,myboxplot3d(sum1.mean_GroupCount,sum1.Group3);
title('TROP2 cluster size (by cells)','FontSize',14);
ylabel('Cell counts');
set(gca,'xticklabels',listGroup3);

%% CD31 in four groups

sum1 = sumAll_2025;

figure,myboxplot4d(sum1.mean_CD31p*100,sum1.GroupOrder);
ytickformat('percentage');
set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('% positive cells');
title('CD31+ cells','FontSize',16);


%% Ecad cluster size & count  (FigS2E)

figure,myboxplot2dj(sumECADpID.mean_GroupCount,sumECADpID.Response);
set(gca,'xticklabels',{'pCR','RD'});
ylabel('Cell count');
title('Ecad+ cluster size');
title('Ecad+ cluster size','FontSize',16);

figure,myboxplot2dj(sumECADpID.GroupCount,sumECADpID.Response);
ylabel('number of clusters');
set(gca,'xticklabels',{'pCR','RD'});
title('Ecad+ cluster count','FontSize',16);

%% Calculate cluster size and Ki67+ ratio

sumLog2TROP2pID = table;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'.....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    sum1 = varfun(@mean,data1,'GroupingVariables','TROP2pID');
    sum1 = sum1(:,1:2);
    data1 = join(data1,sum1,'keys','TROP2pID');
    data1.GroupCount(data1.TROP2pID==0) = 0;
    data2 = data1(data1.GroupCount>0,:);
    data2.log2_TROP2pID = round(log2(data2.GroupCount));
    sum2 = varfun(@mean,data2,'GroupingVariables','log2_TROP2pID');
    sum2.slideName = repmat(slideName(i),size(sum2,1),1);
    if isempty(sumLog2TROP2pID)
        sumLog2TROP2pID = sum2;
    else
        sumLog2TROP2pID = vertcat(sumLog2TROP2pID,sum2);
    end
    toc;
end

sumLog2TROP2pID = join(sumLog2TROP2pID,slideINFO,'keys','slideName');

%% Plot size versus Ki67+ ratio boxplot (Fig2X??)
sum2 = sumLog2TROP2pID;
sum2 = sum2(ismember(sum2.Response,'NR'),:);
figure,boxplot(sum2.mean_Ki_67p*100,sum2.log2_TROP2pID);
ytickformat('percentage');
xlabel('log2(Cluster Size)');
title('Cluster size versus Ki67+ ratio(RD)','FontSize',14);

%% Plot size versus Ki67+ ratio scatter (Fig2X??)
sum2 = sumLog2TROP2pID;

colors = repmat([1,0,0],size(sum2,1),1);
for i = 1:size(sum2,1)
    if ismember(sum2.Response(i),'NR')
        colors(i,:) = [0,0,1];
    end
end

figure,scatter(sum2.log2_TROP2pID,sum2.mean_Ki_67p*100,30,colors,'fill');
ytickformat('percentage');
xlabel('log2(Cluster Size)');
title('Cluster size versus Ki67+ ratio','FontSize',14);

%% Delaunay plots (Fig 2P);

figure,CycIF_tumorview(data1,'TROP2pID',11,0);
hold on;

for i = 1:max(data1.TROP2pID)
    data2 = data1(data1.TROP2pID==i,:);
    if size(data2,1)>3
        DT = delaunayTriangulation([data2.Xt,data2.Yt]);
        triplot(DT,'m');
    end
end

%% Select Topics and Ki67+  (Fig 4L)

sumT2 = sumT1(ismember(sumT1.topics,[11,13,9,10,8]),:);

sumT2.topicgroup = zeros(size(sumT2,1),1);
sumT2.topicgroup(ismember(sumT2.topics,[11,13]))=1;
sumT2.topicgroup(ismember(sumT2.topics,[9,10]))=2;
sumT2.topicgroup(ismember(sumT2.topics,[8]))=3;

figure,myboxplot3(sumT2.mean_Ki_67p*100,sumT2.topicgroup);
ytickformat('percentage');
set(gca,'xticklabels',{'Topic 11/13','Topic 9/10','Topic 8'});
title('Ki67+ cells','FontSize',14);


%% Selected inteaciton  single Boxplot (Fig 4B)

i = find(ismember(labelp,'CD8a'));
j = find(ismember(labelp,'CD68'));

temp1 = allint1r ./ allint1;

figure('units','normalized','outerposition',[0.5 0.5 0.25 0.5]);

myboxplot2dj(squeeze(temp1(i,j,:)),sumAll.Response);
title(strcat(labelp{i},'-',labelp{j},' interaction'),'FontSize',16);  
ylabel('Normalized Interaction','FontSize',12);
set(gcf,'color','w');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

%% TROP2 percentage (Orion versus IHC, Fig S5A)

sum2 = sumAll;
figure, boxplot(sum2.mean_TROP2p*100,sum2.TROP2_percent);
ytickformat('percentage');
ylabel('%TROP2+ (Orion)');
xlabel('%TROP2+ (IHC)');
title('TROP2+ (IHC verus Orion)','FontSize',14);

%% TROP2 (pCR verus RD, Fig S5A);

sum2 = sumAll;

figure('units','normalized','outerposition',[0.5 0.4 0.25 0.6]);

subplot(2,1,1);
myboxplot2dj(sum2.mean_TROP2p*100,sum2.Response);
ytickformat('percentage');
title('TROP2+ (Orion)');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

subplot(2,1,2);
myboxplot2dj(sum2.TROP2_percent,sum2.Response);
ytickformat('percentage');
title('TROP2+ (IHC)');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

%% TROP2 IHC (Fig S5B)

alldata.TROP2z = zscore(alldata.TROP2);
alltumor = alldata(alldata.E_cadherinp,:);

sumAll = sortrows(sumAll,"slideName","ascend");

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.45]);

boxplot(zscore(alltumor.TROP2z),alltumor.slideName,'symbol','');
h = findobj(gca, 'Tag', 'Box');  % Find all box objects
set(gca,'xticklabels',sumAll.PatientID)

for j = 1:length(h)
    if ismember(sumAll.Response(j),'R')
        patch(get(h(j), 'XData'), get(h(j), 'YData'), [1,0,0], 'FaceAlpha', 0.5);
    else
        patch(get(h(j), 'XData'), get(h(j), 'YData'), [0,0,1], 'FaceAlpha', 0.5);
    end
end

title('TROP2 ITH');
ylabel('zscore(TROP2)');

%% Cell counts & density (Fig S5C)

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.45]);

subplot(1,2,1);
myboxplot2dj(sum2.GroupCount/1000,sum2.Response);
title('Raw cell counts','FontSize',14);
ytickformat('%g K');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

subplot(1,2,2);
myboxplot2dj(sum2.GroupCount./arrayArea/1000,sum2.Response);
title('Cell Density','FontSize',14);
ylabel('cells/mm^2');
ytickformat('%g K');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

%% Plot size versus Ki67+ ratio boxplot (FigS5E)
sum2 = sumLog2TROP2pID;
sum2 = sum2(ismember(sum2.Response,'NR'),:);
figure,boxplot(sum2.mean_Ki_67p*100,sum2.log2_TROP2pID);
ytickformat('percentage');
xlabel('log2(Cluster Size)');
title('Cluster size versus Ki67+ ratio(RD)','FontSize',14);

%% Plot size versus Ki67+ ratio scatter (FigS5E)
sum2 = sumLog2TROP2pID;

colors = repmat([1,0,0],size(sum2,1),1);
for i = 1:size(sum2,1)
    if ismember(sum2.Response(i),'NR')
        colors(i,:) = [0,0,1];
    end
end

figure,gscatter(sum2.log2_TROP2pID,sum2.mean_Ki_67p*100,30,colors,'fill');
ytickformat('percentage');
xlabel('log2(Cluster Size)');
title('Cluster size versus Ki67+ ratio','FontSize',14);

%% Plot size versus Ki67+ ratio scatter (FigS5E updated 2025/11/16)
sum2 = sumLog2TROP2pID;
sum2 = sum2(ismember(sum2.Response,'R'),:);

figure,scatter(sum2.log2_TROP2pID,sum2.mean_Ki_67p*100,30,'r','fill');
[r1,p1] = corr(sum2.log2_TROP2pID,sum2.mean_Ki_67p*100);
hold on;

sum2 = sumLog2TROP2pID;
sum2 = sum2(ismember(sum2.Response,'NR'),:);
scatter(sum2.log2_TROP2pID,sum2.mean_Ki_67p*100,30,'b','fill');
[r2,p2] = corr(sum2.log2_TROP2pID,sum2.mean_Ki_67p*100);
h=lsline;
set(h(1),'color','b');
set(h(1),'LineWidth',1);
set(h(2),'color','r');
set(h(2),'LineWidth',1);
text(2.5, 85,strcat('r=',num2str(r1,'%0.2f')),'Color','r');
text(2.5, 81,strcat('p=',num2str(p1,'%0.5f')),'Color','r');
text(13.5, 14,strcat('r=',num2str(r2,'%0.2f')),'Color','b');
text(13.5, 10,strcat('p=',num2str(p2,'%0.5f')),'Color','b');
hold on;

ytickformat('percentage');
xlabel('log2(Cluster Size)');
ylabel('% Ki67+ cells');
title('Cluster size versus Ki67+ ratio','FontSize',14);

%% Plot relationship between TROP2 cluster size and Ki67+ (Fig S5E)

colors = repmat([1,0,0],size(sumAll,1),1);

for i = 1:size(sumAll,1)
    if ismember(sumAll.Response(i),'NR')
        colors(i,:) = [0,0,1];
    end
end

figure,scatter(sumAll.mean_Ki_67p*100,sumTROP2pID.mean_GroupCount,40,colors,'fill');
xtickformat('percentage');
xlabel('Ki67+');

ylabel('TROP2+ cluster size');

%% TLS identification example (Fig S7B)

figure,CycIF_tumorview(dataLSP22301,'CD20pID',9);
title('Algorithmatic TLS identification');
daspect([1 1 1]);

%% Stacked bar plot (LDA topics FigS7C)

sum2 = sumAll;
sum2 = sum2(ismember(sum2.Response,'R'),:);
figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);
subplot(1,3,1);
colors = jet(16);
% ---- legend----
subplot(1,20,1);
h1=bar(1,ones(16,1), 'stacked');
ylabel('topic');
set(gca,'ytick',1:16);

for k = 1:numel(h1)
    h1(k).FaceColor = colors(k,:);
    %h1(k).EdgeColor = [1 1 1];
end
xlim([0.2 1.8]);
ylim([0 16]);
%set(gca,'ytick',0.5:1:10.5);
%set(gca,'yticklabels',labelCT);
%set(gca,'ytick',[]);
set(gca,'XColor','none');
%set(gca,'YColor','none');

% ---- Actual plot (R)---
subplot(1,20,2:7);

h2= bar(sum2{:,72:87},'stacked');
ylim([0 1]);

for k = 1:numel(h2)
    h2(k).FaceColor = colors(k,:);
    %h2(k).EdgeColor = [1 1 1];
end
ylim([0 1]);
set(gca,'ytick',[]);
xlim([0.5 14.5]);
set(gca,'xtick',1:14);
set(gca,'xticklabels',strrep(sum2.PatientID,'#','P'));
title('pCR');



% ---- Actual plot (R)---
sum2 = sumAll;
sum2 = sum2(ismember(sum2.Response,'NR'),:);

subplot(1,20,8:20);

h2= bar(sum2{:,72:87},'stacked');
ylim([0 1]);

for k = 1:numel(h2)
    h2(k).FaceColor = colors(k,:);
    %h2(k).EdgeColor = [1 1 1];
end
ylim([0 1]);
set(gca,'ytick',[]);
xlim([0.5 34.5]);
set(gca,'xtick',1:34);
set(gca,'xticklabels',strrep(sum2.PatientID,'#','P'));
title('RD');
set(gcf,'color','w');

%% Select Topics and Ki67+  (Fig S7D)

sumT2 = sumT1(ismember(sumT1.topics,[11,13,9,10,8]),:);

sumT2.topicgroup = zeros(size(sumT2,1),1);
sumT2.topicgroup(ismember(sumT2.topics,[11,13]))=1;
sumT2.topicgroup(ismember(sumT2.topics,[9,10]))=2;
sumT2.topicgroup(ismember(sumT2.topics,[8]))=3;

figure,myboxplot3(sumT2.mean_TROP2p*100,sumT2.topicgroup);
ytickformat('percentage');
set(gca,'xticklabels',{'Topic 11/13','Topic 9/10','Topic 8'});
title('TROP2+ cells','FontSize',14);

%%  Total CD31 clusters (Fig S9B)

sum1 = varfun(@mean,tableVesM,'GroupingVariables','slideName');
sum1 = join(sum1,slideINFO,'keys','slideName');

sum1 = sortrows(sum1,'GroupCount',"ascend");

figure('units','normalized','outerposition',[0.5 0 0.5 1]);

b= barh(sum1.GroupCount,'FaceColor',[0,0,1],'EdgeColor',[1,1,1],'LineWidth',0.1);
hold on;
b.FaceColor = 'flat';

for i = 1:size(sum1,1)
    b.CData(i,:) = [0 0 1];
end

for i = 1:size(sum1,1)
    if ismember(sum1.Response(i),'R')
        b.CData(i,:) = [1 0 0];
    end
end

h = zeros(2, 1);
h(1) = barh(NaN,'Red');
h(2) = barh(NaN,'Blue');
legend(h, {'pCR','RD'},'Location','southeast','FontSize',16);

title('Total CD31+ Clusters','FontSize',14,'Interpreter','none');
set(gcf,'color','w');
set(gca,'ytick',1:48);
set(gca,'yticklabels',sum1.PatientID);
xlabel('Counts');

%%  Total CD31+ SMA+ clusters (Fig S9B)

sum1 = varfun(@mean,tableVesM(tableVesM.VesM==1,:),'GroupingVariables','slideName');
sum1 = join(sum1,slideINFO,'keys','slideName');

sum1 = sortrows(sum1,'GroupCount',"ascend");

figure('units','normalized','outerposition',[0.5 0 0.5 1]);

b= barh(sum1.GroupCount,'FaceColor',[0,0,1],'EdgeColor',[1,1,1],'LineWidth',0.1);
hold on;
b.FaceColor = 'flat';

for i = 1:size(sum1,1)
    b.CData(i,:) = [0 0 1];
end

for i = 1:size(sum1,1)
    if ismember(sum1.Response(i),'R')
        b.CData(i,:) = [1 0 0];
    end
end

h = zeros(2, 1);
h(1) = barh(NaN,'Red');
h(2) = barh(NaN,'Blue');
legend(h, {'pCR','RD'},'Location','southeast','FontSize',16);

title('CD31+/SMA+ Clusters','FontSize',14,'Interpreter','none');
set(gcf,'color','w');
set(gca,'ytick',1:48);
set(gca,'yticklabels',sum1.PatientID);
xlabel('Counts');


%%  Total CD31 * CD31+SMA+ clusters (Fig S9B merged)

sum1 = varfun(@mean,tableVesM,'GroupingVariables','slideName');
sum1 = join(sum1,slideINFO,'keys','slideName');

sum1 = sortrows(sum1,'GroupCount',"ascend");

figure('units','normalized','outerposition',[0.5 0 0.5 1]);

b= barh(sum1.GroupCount,'FaceColor',[0,0,1],'EdgeColor',[1,1,1],'LineWidth',0.1);
hold on;
b.FaceColor = 'flat';

for i = 1:size(sum1,1)
    b.CData(i,:) = [0 0 1];
    b.FaceAlpha = 0.4;
end

for i = 1:size(sum1,1)
    if ismember(sum1.Response(i),'R')
        b.CData(i,:) = [1 0 0];
    end
end

title('Total CD31+/SMA+ Clusters','FontSize',14,'Interpreter','none');
set(gcf,'color','w');
set(gca,'ytick',1:48);
set(gca,'yticklabels',sum1.PatientID);
xlabel('Counts');
hold on;

%----
sum2 = varfun(@mean,tableVesM(tableVesM.VesM==1,:),'GroupingVariables','slideName');
sum2 = join(sum2,slideINFO,'keys','slideName');
sum2.Temp1 = sum2.GroupCount;
sum2 = sum2(:,{'slideName','Temp1'});

sum1 = join(sum1,sum2,'keys','slideName');

b2= barh(sum1.Temp1,'FaceColor',[0,0,1],'EdgeColor',[1,1,1],'LineWidth',0.1);
b2.FaceColor = 'flat';

for i = 1:size(sum1,1)
    b2.CData(i,:) = [0 0 1];
end

for i = 1:size(sum1,1)
    if ismember(sum1.Response(i),'R')
        b2.CData(i,:) = [1 0 0];
    end
end

h = zeros(4, 1);
h(1) = barh(NaN,'Red');
h(2) = barh(NaN,'Blue');
h(3) = barh(NaN,'Red');
h(4) = barh(NaN,'Blue');

legend(h, {'pCR','RD','pCR','RD'},'Location','southeast','FontSize',16);

%% Density versus Ki67 (FigS5X)

sum2 = sumAll_2025;
sum2.Density = sum2.GroupCount./arrayArea;

colors = repmat([1,0,0],size(sum2,1),1);
for i = 1:size(sum2,1)
    if ismember(sum2.Response(i),'NR')
        colors(i,:) = [0,0,1];
    end
end

figure,scatter(sum2.Density,sum2.mean_Ki_67p*100,30,colors,'fill');
ytickformat('percentage');
xlabel('Cells/mm^2');
title('Cell Density versus Ki67+ ratio','FontSize',14);

%% Density versus CD45+ (FigS5X)

sum2 = sumAll_2025;
sum2.Density = sum2.GroupCount./arrayArea;

colors = repmat([1,0,0],size(sum2,1),1);
for i = 1:size(sum2,1)
    if ismember(sum2.Response(i),'NR')
        colors(i,:) = [0,0,1];
    end
end

figure,scatter(sum2.Density,sum2.mean_Ki_67p*100,30,colors,'fill');
h1=lsline;
h1.Color = 'k';
h1.LineWidth = 0.5;

[r, p]= corr(sum2.Density,sum2.mean_Ki_67p*100);
ylims = ylim;
xlims = xlim;
text (0.8*xlims(2),0.9*ylims(2),strcat('r=',num2str(r,'%0.3f')));
text (0.8*xlims(2),0.8*ylims(2),strcat('p=',num2str(p,'%0.3f')));


ytickformat('percentage');
xlabel('Cells/mm^2');
title('Cell Density versus Ki67+ ratio','FontSize',14);

%% Read TROP2 memembrane/Cytosol ratio from imageJ results

% dir1 = dir("*.csv");
% dir1 = struct2table(dir1);
% dir1 = sortrows(dir1,"name","ascend");
% listResults = dir1.name;

arrayTROP2ratios = zeros(length(slideName),2);  % col1 = memebrane, col2 = cytosol

for i = 1:length(listResults)
    table1 = readtable(listResults{i});
    arrayTROP2ratios(i,1)=table1.Mean(1);
    arrayTROP2ratios(i,2)=table1.Mean(2);
end
disp('Done');

%% Figure XX adjusted TROP2 ratios


figure,myboxplot2dj(sumAll_2025.TROP2ratios .* sumAll_2025.mean_TROP2p .* sumAll_2025.mean_E_cadherinp,sumAll_2025.Response);
title('TROP2 cytosol/membrane ratio','FontSize',14);
ylabel('Adjusted ratios');
set(gca,'xticklabels',{'pCR','RD'},'FontSize',14);

%% Figure XX TROP2 ratios in four groups (label)

sum1 = sumAll_2025;
sum2 = sum1(sum1.Recurr,:);
figure,myboxplot4d(sum1.TROP2ratios .* sum1.mean_TROP2p .* sum1.mean_E_cadherinp ,sum1.GroupOrder);
hold on;
scatter(sum2.GroupOrder,sum2.TROP2ratios .* sum2.mean_TROP2p .* sum2.mean_E_cadherinp ,50,'green','fill');
x = sum2.GroupOrder;
y = sum2.TROP2ratios .* sum2.mean_TROP2p .* sum2.mean_E_cadherinp;
txt = sum2.PatientID;

for i = 1:size(sum2,1)
    text(x(i)+0.1,y(i),txt(i));
end

set(gca,'xticklabels',tableGroupOrder.Group);
ylabel('Adjusted Ratios');
title('TROP2 cytosol/membrane ratio','FontSize',14);


%% Figure S5 (Scatter plot)

colors = repmat([0,0,1],48,1);
colors(ismember(sumAll.Response,'R'),:) = repmat([1,0,0],14,1);
figure,scatter(sumAll_2025.TROP2_percent,sumAll_2025.mean_TROP2p*100,40,colors,'fill');
lsline;
ytickformat('percentage');
xtickformat('percentage');
ylabel('Orion');
xlabel('IHC scores');
title('TROP2 quantification','FontSize',14);

temp1 = jonckheere_terpstra(sumAll_2025.mean_TROP2p, sumAll_2025.TROP2_percent);
text (10,70,strcat('JT z=',num2str(temp1.z)));
text (10,65,strcat('JT p=',num2str(temp1.p)));


%% Figure XXX TROP2 cluster size and cell density plot


figure,scatter(sumAll.GroupCount./arrayArea,sumTROP2pID.mean_GroupCount,40,colors,'fill');
xlabel('Cell Density (cells/mm^2)');
[r,p]= corr(sumAll.GroupCount./arrayArea,sumTROP2pID.mean_GroupCount);

ylabel('TROP2+ cluster size(cells)');
text(10000,2200,strcat('r=',num2str(r)));
text(10000,2000,strcat('p=',num2str(p)));