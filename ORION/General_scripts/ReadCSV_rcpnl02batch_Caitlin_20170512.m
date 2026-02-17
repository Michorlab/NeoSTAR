%% Read RPCNL QUANT (for Caitlin xenograft)
%  required the variables :: labels
%  v2: require user inputs (name, cols, rows,chs); rename alldata;
%  Jerry Lin 2016/12/14


%% Initialization


myDIR = uigetdir('C:\IMAGETEMP\CAITLIN-20170512');

slideName = {'CMD1N','CMD2N','CMD3N','CMD4N','CMD5N','CMD6N','CMD7N','CMD8N','CMD9N','CMD10N','CMD11N','CMD12N'};

colarray = [7,7,5,3,8,6,8,5,7,9,6,7];
rowarray = [4,5,5,6,9,6,8,3,6,7,6,7];

xlims = 1664;    % dimension for 10x image on RareCyte
ylims = 1404;    % dimension for 10x image on RareCyte

samples = 20000;
allsample = table;

%% Reading each slide

for slide = 1:length(slideName)
    
rows = rowarray(slide);
cols = colarray(slide);

myName = slideName{slide};

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
temp1.Rows = ceil(temp1.Slice/cols);
temp1.Cols = temp1.Slice - (temp1.Rows-1)*cols;
temp1.Xt = temp1.X + (temp1.Cols-1)* xlims;
temp1.Yt = temp1.Y + (temp1.Rows-1)* ylims;
sample1 = datasample(temp1,samples);

eval(strcat('data_',myName,' = temp1;'));
eval(strcat('sample_',myName,' = sample1;'));

sample1.sname = repmat({myName},size(sample1,1),1);

if(isempty(allsample))
    allsample = sample1;
else
    allsample = vertcat(allsample,sample1);
end

clear dataDAPI dataFITC dataCy3 dataCy5 temp1 sample1;

end  %for slide