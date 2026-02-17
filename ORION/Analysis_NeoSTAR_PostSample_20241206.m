%% Analysis NeoSTAR post (Orion)
%  Jerry Lin 2024/12/06

%% Import data

dir1 = dir;
dir1 = struct2table(dir1);
dir1 = dir1.name;

filelist = dir1;
filelist = filelist(3:end);
slideName = cellfun(@(X) X(1:8),filelist,'UniformOutput',false);

tic;
for i = 1:length(slideName)
    filename = filelist{i};
    disp(strcat('Processing:',slideName{i}));
    data1 = CycIF_importMcMicro(filename,0.325);
    data1 = CycIF_removezero(data1);
    data1{:,2:end} = uint16(data1{:,2:end});
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Change channel names;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',num2str(i),':',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));
    data1.Properties.VariableNames(2) = "DNA";
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end


%% Reorder channels

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',num2str(i),':',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));
    data1 = data1(:,label1);
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end
%% test

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

%% Boxplot (all gated markers)

sum2 = sumTumor;
title1 = 'Tumor region';

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

%% Boxplot All double gates

sum3 = sumTumor;
title1 = 'Tumor region';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(doubleGates,1)
    marker1 = strcat('mean_',doubleGates{i,1},'p',doubleGates{i,2},'p');
    subplot(3,8,i);
    myboxplot2d(sum3{:,marker1}*100,sum3.Response);
    ytickformat('percentage');
    title(strcat(doubleGates{i,1},'+',doubleGates{i,2},'+'),'Interpreter','none');
end
sgtitle(title1);
set(gcf,'color','w');

%% Boxplot (all topics)

sum2 = sumAll_post;
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

%% Merge summary tables

temp1 = sumAll_pre(:,temp3);
temp1.Treat = repmat({'Pre'},size(temp1,1),1);
temp2 = sumAll_post(:,temp3);
temp2.Treat = repmat({'Post'},size(temp2,1),1);
sumAll_prepost = vertcat(temp1,temp2);


temp1 = sumTumor_pre(:,temp3);
temp1.Treat = repmat({'Pre'},size(temp1,1),1);
temp2 = sumTumor_post(:,temp3);
temp2.Treat = repmat({'Post'},size(temp2,1),1);
sumTumor_prepost = vertcat(temp1,temp2);

temp1 = sumBorder_pre(:,temp3);
temp1.Treat = repmat({'Pre'},size(temp1,1),1);
temp2 = sumBorder_post(:,temp3);
temp2.Treat = repmat({'Post'},size(temp2,1),1);
sumBorder_prepost = vertcat(temp1,temp2);


%% Generate paired summary

temp1 = sumAll_preS(:,temp3);
temp1.Treat = repmat({'Pre'},size(temp1,1),1);
temp2 = sumAll_postS(:,temp3);
temp2.Treat = repmat({'Post'},size(temp2,1),1);
sumAll_prepostS = vertcat(temp1,temp2);
%% Boxplot (all gated markers)

sum2 = sumAll_prepost;
title1 = 'Whole sample';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    subplot(3,6,i);
    myboxplot2d(sum2{:,marker1}*100,sum2.Treat);
    ytickformat('percentage');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);

%% Boxplot (all gated markers only pre

sum2 = sumAll_prepost;
sum2 = sum2(ismember(sum2.Response,'R'),:);
title1 = 'Only R';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    marker1 = strcat('mean_',labelp2{i});
    subplot(3,6,i);
    myboxplot2d(sum2{:,marker1}*100,sum2.Treat);
    ytickformat('percentage');
    title(labelp3{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);


%% Boxplot (all topics)

sum2 = sumAll_prepost;
sum2 = sum2(ismember(sum2.Response,'R'),:);
title1 = 'Only R';

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:maxT
    marker1 = strcat('mean_topic',num2str(i));
    subplot(3,6,i);
    myboxplot2d(sum2{:,marker1}*100,sum2.Treat);
    ytickformat('percentage');
    title(strcat('Topic',num2str(i)));
end

set(gcf,'color','w');
sgtitle(title1);

%% Generate VESp

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',num2str(i),':',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));
    data1.VESp = data1.CD31p | data1.SMAp;
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end


%% Paired Boxplot (pre-post R-NR)

figure,myboxplot2p(reshape(sumAll_prepostS.mean_TROP2p,2,30)',groupcolor);
ytickformat('percentage');
title('LDA topic','FontSize',18);
set(gca,'xticklabels',{'Pre','Post'});

%% Paired Boxplot (pre-post R-NR)
name1 = 'Ki_67';
title1 = 'Ki67+';
marker1 = strcat('mean_',name1,'p');

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);

subplot(1,2,1);
sum2 = sumAll_prepostS(ismember(sumAll_prepostS.Response,'NR'),:);
groupc = repmat({'b'},size(sum2,1)/2,1);
myboxplot2p(reshape(sum2{:,marker1}*100,2,size(sum2,1)/2)',groupc);
ytickformat('percentage');
title('RD','FontSize',16);
set(gca,'xticklabels',{'Pre','Post'});
subplot(1,2,2);
sum2 = sumAll_prepostS(ismember(sumAll_prepostS.Response,'R'),:);
groupc = repmat({'r'},size(sum2,1)/2,1);
myboxplot2p(reshape(sum2{:,marker1}*100,2,size(sum2,1)/2)',groupc);
ytickformat('percentage');
title('pCR','FontSize',16);
sgtitle(title1,'FontSize',16,'Interpreter','none');
set(gca,'xticklabels',{'Pre','Post'});

set(gcf,'color','w');

%% Paired Boxplot (pre-post R-NR)
name1 = 'CD8apKi_67';
title1 = 'CD8+Ki67+';
marker1 = strcat('mean_',name1,'p');

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);

subplot(1,2,1);
sum2 = sumTumor_prepostS(ismember(sumTumor_prepostS.Response,'NR'),:);
groupc = repmat({'b'},size(sum2,1)/2,1);
myboxplot2p(reshape(sum2{:,marker1}*100,2,size(sum2,1)/2)',groupc);
ytickformat('percentage');
title('RD','FontSize',16);
set(gca,'xticklabels',{'Pre','Post'});
subplot(1,2,2);
sum2 = sumTumor_prepostS(ismember(sumTumor_prepostS.Response,'R'),:);
groupc = repmat({'r'},size(sum2,1)/2,1);
myboxplot2p(reshape(sum2{:,marker1}*100,2,size(sum2,1)/2)',groupc);
ytickformat('percentage');
title('pCR','FontSize',16);
sgtitle(title1,'FontSize',16,'Interpreter','none');
set(gca,'xticklabels',{'Pre','Post'});

set(gcf,'color','w');

%% Assign to Celltype all slides

labelCT = {'Tumor1','Tumor2','Tumor3','Tumor4','Immune1','Immune2','Immune3','Immune4','Immune5','Stroma1','Stroma2'};
labelCTname = {'Tumor1(Ki67+)','Tumor2(Ki67-)','Tumor3(TROP2-)','Tumor4(MIX)','Immune1(CD4 T)','Immune2(CD8 T)','Immune3(Myeloid)','Immune4(CD163+)','Immune5(B cells)','Stroma1(CD31+)','Stroma2(Other)'};

labelCC = {'Tumor','Immune','Stroma'};
Catlist = [1 1 1 1 2 2 2 2 2 3 3];

tic;
for s = 1:length(slideName)
    disp(strcat('Processing:',slideName{s},'...',num2str(s)));
    data1 = eval(strcat('data',slideName{s}));
    data1.CellType=trainedModel_tree2.predictFcn(data1);
    for i = 1:11
        data1.CellCat(data1.CellType==i)= Catlist(i);
    end
    for i = 1:length(labelCT)
        data1{:,labelCT{i}} = data1.CellType ==i;
    end
    for j = 1:length(labelCC)
        data1{:,labelCC{j}} = data1.CellCat == j;
    end
    eval(strcat('data',slideName{s},'=data1;'));
    toc;
end

%% only RD sample

name1 = 'CD4pPD_1';
title1 = 'CD4+PD1+ cells';
marker1 = strcat('mean_',name1,'p');

figure;


sum2 = sumAll_prepostS(ismember(sumAll_prepostS.Response,'NR'),:);
groupc = repmat({'b'},size(sum2,1)/2,1);

myboxplot2p(reshape(sum2{:,marker1}*100,2,size(sum2,1)/2)',groupc);
ytickformat('percentage');
title(title1,'FontSize',16);
set(gca,'xticklabels',{'Pre','Post'});

%% PrePost analysis (cell-cell interaction)
% calculate all cluster score

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

%% PrePost analysis (cell-cell interaction)
% check cluster score per sample

i = 2;

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

%% PrePost analysis (cell-cell interaction)
% Test all pvalue (cluster score)

arrayClusterPvalues = ones(length(labelp),length(labelp));
arrayClusterFolds = zeros(length(labelp), length(labelp));
tic;
for i = 1:length(labelp)
    for j = 1:length(labelp)
        [arrayClusterPvalues(i,j),arrayClusterFolds(i,j)]=myttest2(squeeze(arrayClusterScore(i,j,:)),sumAll_prepostS.Treat);
    end
end
toc;

%% Plot interaction (pre/post)

Xpos = repmat((1:length(labelp))',1,length(labelp));
Ypos = repmat(1:length(labelp),length(labelp),1);

figure,scatter(Xpos(:),Ypos(:),-log(arrayClusterPvalues(:))*5,arrayClusterFolds(:),'fill');
xlim([0.5 length(labelp)+2.5]);
ylim([0.5 length(labelp)+0.5]);
hold on;
scatter(length(labelp)+1,2,-log(0.01)*5,'k');
text(length(labelp)+0.5,1,'p=0.01','FontSize',6);
scatter(length(labelp)+1,4,-log(0.001)*5,'k');
text(length(labelp)+0.5,3,'p=0.001','FontSize',6);
scatter(length(labelp)+1,6,-log(0.0001)*5,'k');
text(length(labelp)+0.5,5,'p=0.0001','FontSize',6);

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
colormap(redbluecmap);
colorbar;
caxis([0 2]);

%% PrePost analysis (cell-cell interaction)
% check 

name1 = 'CD8a';
name2 = 'SMA';

i = find(ismember(labelp,name1));
j = find(ismember(labelp,name2));

figure,myboxplot2d(squeeze(arrayClusterScore(i,j,:)),sumAll.Response);
title(strcat(name1,'-',name2),'FontSize',18);

%% remove extra data set

list1 = slideName_all(~ismember(slideName_all,slideName));
for i = 1:length(list1)
    eval(strcat('clear data',list1{i}));
end


%% PrePost analysis (cell-cell interaction, only RD)
% Test all pvalue (cluster score)

arrayClusterPvaluesRD = ones(length(labelp),length(labelp));
arrayClusterFoldsRD = zeros(length(labelp), length(labelp));
tic;
for i = 1:length(labelp)
    for j = 1:length(labelp)
        [arrayClusterPvaluesRD(i,j),arrayClusterFoldsRD(i,j)]=myttest2(squeeze(arrayClusterScore(i,j,ismember(sumAll_prepostS.Response,'NR'))),sumAll_prepostS.Treat(ismember(sumAll_prepostS.Response,'NR')));
    end
end
toc;

%% Plot interaction (Pre/post only RD)

Xpos = repmat((1:length(labelp))',1,length(labelp));
Ypos = repmat(1:length(labelp),length(labelp),1);

figure,scatter(Xpos(:),Ypos(:),-log(arrayClusterPvaluesRD(:))*5,arrayClusterFoldsRD(:),'fill');
xlim([0.5 length(labelp)+2.5]);
ylim([0.5 length(labelp)+0.5]);
hold on;
scatter(length(labelp)+1,2,-log(0.01)*5,'k');
text(length(labelp)+0.5,1,'p=0.01','FontSize',6);
scatter(length(labelp)+1,4,-log(0.001)*5,'k');
text(length(labelp)+0.5,3,'p=0.001','FontSize',6);
scatter(length(labelp)+1,6,-log(0.0001)*5,'k');
text(length(labelp)+0.5,5,'p=0.0001','FontSize',6);

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
colormap(redbluecmap);
colorbar;
caxis([0 2]);

%% Paired Boxplot (pre-post R-NR)
%name1 = 'CD20';
%title1 = 'CD20+';
%marker1 = strcat('mean_',name1,'p');

figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);

subplot(1,2,1);
sum2 = sumAll_prepostS(ismember(sumAll_prepostS.Response,'NR'),:);
groupc = repmat({'b'},size(sum2,1)/2,1);
myboxplot2p(reshape(sum2.mean_CD8ap./sum2.mean_FOXP3p,2,size(sum2,1)/2)',groupc);
%ytickformat('percentage');
title('RD','FontSize',16);
set(gca,'xticklabels',{'Pre','Post'});
subplot(1,2,2);
sum2 = sumAll_prepostS(ismember(sumAll_prepostS.Response,'R'),:);
groupc = repmat({'r'},size(sum2,1)/2,1);
myboxplot2p(reshape(sum2.mean_CD8ap./sum2.mean_FOXP3p,2,size(sum2,1)/2)',groupc);
%ytickformat('percentage');
title('pCR','FontSize',16);
sgtitle('CD8+/FOXP3+','FontSize',16,'Interpreter','none');
set(gca,'xticklabels',{'Pre','Post'});

set(gcf,'color','w');


%% Display LDA topics

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

%% paired boxplot for all gated markers (Fig12B)

sum2 = sumAll_prepostS(ismember(sumAll_prepostS.Response,'NR'),:);
sum2 = sortrows(sum2,"slideName","ascend");
sum2 = sortrows(sum2,"PatientID","ascend");
groupc = repmat({'b'},size(sum2,1)/2,1);


figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(labelp)
    subplot(3,6,i);
    marker1 = strcat('mean_',labelp2{i});
    myboxplot2p(reshape(sum2{:,marker1}*100,2,size(sum2,1)/2)',groupc);
    ytickformat('percentage');
    title(labelp3{i},'Interpreter','none');
    set(gca,'xticklabels',{'Pre','Post'});
end

set(gcf,'color','w');

%% Assign TLS counts (FigS12X)

list1 = sumAll_prepostS.slideName;
TLScount = zeros(length(list1),1);
for i = 1:length(list1)
    data1 =eval(strcat('data',list1{i}));
    TLScount(i) = max(data1.TLS);
end
sumAll_prepostS.TLScount = TLScount;

sum2 = sumAll_prepostS;
groupc = repmat({'b'},size(sum2,1)/2,1);
figure,myboxplot2p(reshape(sum2.TLScount,2,size(sum2,1)/2)',groupc);
set(gca,'xticklabels',{'Pre','Post'});
title('TLS count','FontSize',14);