%% Read RPCNL QUANT (for Caitlin xenograft)
%  required the variables :: labels
%  v2: require user inputs (name, cols, rows,chs); rename alldata;
%  Jerry Lin 2016/12/14


%% Initialization


myDIR = uigetdir('F:\Caitlin-20170217');
rows = 15;
cols = 16;
xlim = 1664;
ylim = 1404;
samples = 20000;


%% start reading csv filess

myName = input('Please input file name:','s');
%%eachdata = cell(totalframe);

filename = strcat(myDIR,'\Results-',myName,'-DAPI.csv');
dataDAPI = importcsv_rpcnl(filename,2,inf);

filename = strcat(myDIR,'\Results-',myName,'-FITC.csv');
dataFITC = importcsv_rpcnl(filename,2,inf);

filename = strcat(myDIR,'\Results-',myName,'-Cy3.csv');
dataCy3 = importcsv_rpcnl(filename,2,inf);

filename = strcat(myDIR,'\Results-',myName,'-Cy5.csv');
dataCy5 = importcsv_rpcnl(filename,2,inf);

%% Assign other parameters

temp1 = dataDAPI;
temp1.DAPI = dataDAPI.Mean;
temp1.FITC = dataFITC.Mean;
temp1.Cy3 = dataCy3.Mean;
temp1.Cy5 = dataCy5.Mean;
temp1.Row = ceil(temp1.Slice/cols);
temp1.Col = temp1.Slice - (temp1.Row-1)*cols;
temp1.Xt = temp1.X + (temp1.Col-1)* xlim;
temp1.Yt = temp1.Y + (temp1.Row-1)* ylim;
sample1 = datasample(temp1,samples);

eval(strcat('data_',myName,' = temp1;'));
eval(strcat('sample_',myName,' = sample1;'));

clear dataDAPI dataFITC dataCy3 dataCy5 temp1 sample1;
clear xlim ylim;

