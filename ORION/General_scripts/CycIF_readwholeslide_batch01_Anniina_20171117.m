%% Batch Read slides (from ImageJ data)
%  required the variables :: labels
% 
%  Jerry Lin 2017/08/24

%% Initialization

slideDir = {'C:\IMAGETEMP\Parental_T2\',...
'C:\IMAGETEMP\Parental_T3\',...
'C:\IMAGETEMP\Parental_T4\',...
'C:\IMAGETEMP\SCP01B_T1\',...
'C:\IMAGETEMP\SCP01B_T2\',...
'C:\IMAGETEMP\SCP01B_T3\',...
'C:\IMAGETEMP\SCP3B_T1\',...
'C:\IMAGETEMP\SCP3B_T2\',...
'C:\IMAGETEMP\SCP3B_T3\',...
'C:\IMAGETEMP\SCP17B_T1\',...
'C:\IMAGETEMP\SCP17B_T2\',...
'C:\IMAGETEMP\SCP17B_T3\',...
'C:\IMAGETEMP\SCP29_T1\',...
'C:\IMAGETEMP\SCP29_T2\',...
'C:\IMAGETEMP\SCP29_T3\',...
'C:\IMAGETEMP\SCP32_T1\',...
'C:\IMAGETEMP\SCP32_T2\',...
'C:\IMAGETEMP\scp32_T3\'};

%input slide name
slideName = {'Parental_T2 ',...
'Parental_T3',...
'Parental_T4',...
'SCP01B_T1',...
'SCP01B_T2',...
'SCP01B_T3',...
'SCP3B_T1',...
'SCP3B_T2',...
'SCP3B_T3',...
'SCP17B_T1',...
'SCP17B_T2',...
'SCP17B_T3',...
'SCP29_T1',...
'SCP29_T2',...
'SCP29_T3',...
'SCP32_T1',...
'SCP32_T2',...
'scp32_T3'};

%input columns & rows of each slides
colarray= [11,11,11,8,9,14,12,13,15,23,13,11,15,9,7,12,10,10,22,9,7,10,6,9,27,18,9,10,8];
rowarray = [8,7,9,7,9,6,10,12,8,9,6,8,12,9,13,13,11,16];

xlims = 416;    % dimension for 20x image on RareCyte
ylims = 351;    % dimension for 20x image on RareCyte

samples = 10000;
chs = 40;
int_cut = 2000;


%alldata = cell(rows,cols);
allsample = table;

%% Reading each slide

for slide = 1:length(slideName)
    
myDIR = slideDir{slide};

myName = slideName{slide}; %input('Please input file name:','s');
cols = colarray(slide); %input('Columns=');
rows = rowarray(slide); %input('Rows=');

alldata =table;
totalframe = rows*cols; 

for i=1:totalframe
  filename = strcat(myDIR,'Results-',myName,'-',num2str(i),'.csv');
  if exist(filename,'file')
    temp1 = array2table(CycIF_readtable03(chs,filename),'VariableNames',labels);
    temp1 = CycIF_filterbyhoechst02(temp1,1:(chs/4-2),int_cut);
    temp1.frame = repmat(i,length(temp1.X),1);
    
    r = floor((i-1)/cols)+1;
    c = i - (r-1)*cols;
    
    temp1.COL = repmat(c,length(temp1.X),1);
    temp1.ROW = repmat(r,length(temp1.X),1);
    temp1.Xt = temp1.X + (c-1)* xlims;
    
    temp1.Yt = temp1.Y + (r-1)* ylims;
    
    
    if isempty(alldata)
        alldata = temp1;
        %%eachdata{i} = temp1;
    else
        alldata = vertcat(alldata,temp1);
        %%eachdata{i} = temp1;
    end
  end  
  display(['Processing:',filename]);
end

myName = strrep(myName,'-','');
sample1 = datasample(alldata,samples);
eval(strcat('data',myName,'=alldata;'));
eval(strcat('sample',myName,'=sample1;'));

sample1.slidename = repmat({myName},length(sample1.X),1);
if(isempty(allsample))
    allsample = sample1;
else
    allsample = vertcat(allsample,sample1);
end

clear alldata sample1;

end
clear samples rows cols slide temp1 totalframe r i c ch chs 