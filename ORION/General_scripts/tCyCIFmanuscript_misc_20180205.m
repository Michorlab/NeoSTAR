%% t-CYCIF manuscript, figures for 16cycles
%  Jerry Lin 2018/02/05


%% Side by side Ab comparision

dataA = sampleM36459A;
dataB = sampleM36459B;

figure;

subplot(2,2,1);
CycIF_tumorview(dataA,'Cy5_4',1);
xlim([0 25000]);
ylim([0 16000]);


subplot(2,2,2);
CycIF_tumorview(dataB,'Cy5_14',1);
xlim([0 25000]);
ylim([0 16000]);

subplot(2,2,3);
findcutoff(log(dataA.FITC_4),2,0.05);
xlim([6 11]);

subplot(2,2,4);
findcutoff(log(dataB.FITC_14),2,0.05);
xlim([6 11]);

%% 2D scatter between Cycles (within same sample)

data1 =sampleTONSILB;

markers = {'Cy3_3','Cy3_7','Cy3_12','Cy3_16'};
cycles = {'Cycle 3','Cycle 7','Cycle 12','Cycle 16'};
labelname = 'Vimentin';

figure;
fig = 1;
for i = 1:4
    for j=i+1:4
       subplot(2,3,fig);
       fig= fig+1;
       dscatter(log(data1{:,markers{i}}),log(data1{:,markers{j}}));colormap(jet);
       ratio = mean(log(data1{:,markers{i}})./ log(data1{:,markers{j}}));
       corr1 = corr(log(data1{:,markers{i}}),log(data1{:,markers{j}}));
       xlabel(strcat(labelname,{' in '},cycles(i),{' (log scale)'}));
       ylabel(strcat(labelname,{' in '},cycles(j),{' (log scale)'}));
       
       %title(strcat(cycles{i},{' versus '},cycles{j}));
       title(['Ratio=',num2str(ratio,'%0.2f'),' rho=',num2str(corr1,'%0.2f')])
    end
end

%% Compare between Sample A&B

dataname = 'TONSIL';

dataA = eval(strcat('sample',dataname,'A'));
dataB = eval(strcat('sample',dataname,'B'));

ch = 'Cy3';
cyA = 5;
cyB = 15;
chA = strcat(ch,'_',num2str(cyA));
chB = strcat(ch,'_',num2str(cyB));

labelname = Abtable{cyA,strcat(ch,'_A')};

figure;

data1 = log(dataA{:,chA});
data2 = log(dataB{:,chB});
min1 = prctile([data1; data2],0.5);
max1 = prctile([data1; data2],99.5);
nbins = min1:(max1-min1)/200:max1;

h1 = smooth(histc(data1,nbins),7);
plot(nbins,h1,'g');
hold on;
h2 = smooth(histc(data2,nbins),7);
plot(nbins,h2,'r');
hold off;

Hcomb = max(h1,h2);
Hover = min(h1,h2);

overarea = trapz(nbins,Hover)/trapz(nbins,Hcomb);


% h1=histfit(log(dataA{:,chA}),100,'kernel');
% delete(h1(1));
% h1(2).Color = 'g';
% hold on;
% 
% h2=histfit(log(dataB{:,chB}),100,'kernel');
% delete(h2(1));
% h2(2).Color = 'r';
% 
legend(['Cycle ',num2str(cyA)],['Cycle ',num2str(cyB)]);
xlabel(strcat(labelname,{' intensity (log scale)'}));
ylabel('Frequency');
title(strcat(dataname,'(overlap=',num2str(overarea,'%0.2f'),')'));











