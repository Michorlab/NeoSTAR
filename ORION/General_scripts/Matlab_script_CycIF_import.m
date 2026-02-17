%% Import CycIF data from text files.
% Script for importing data from the following text file:
%
%   Results-[Row][Col][Filed].txt
%
%   Jerry Lin 2015/09/21
%

%% Initialize variables.
rows = {'A','B','C','D','E','F','G','H'};
cols = {'01','02','03','04','05','06','07','08','09','10','11','12'};
flds = {'fld01','fld02','fld03','fld04','fld05','fld06','fld07','fld08','fld09','fld10','fld11','fld12'};

mypath = uigetdir;
sitedata = cell(8,12,6);
welldata = cell(8,12);
wellsum = cell(8,12);

for r=2:7;
    for c=2:11;
        for f=1:12;
            site = strcat(rows(r),cols(c),flds(f));
            well = strcat(rows(r),cols(c));
            filename = strjoin(strcat(mypath,'\Results-',site,'.txt'));
            

disp(strcat('Processing:',filename));
%%filename = 'C:\CycIF\Results-B02fld01.txt';
delimiter = '\t';
startRow = 3;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: text (%s)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%s%f%f%f%f%f%f%f%f%f%f%f%s%s%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
temp = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','Label','Area','Mean','StdDev','Min','Max','X','Y','Perim','Circ','IntDen','Median','VarName14','RawIntDen','Slice','AR','Round','Solidity'});
if(f<10)
    temp.Channel = cellfun(@(x) x(20:length(x)),temp.Label,'UniformOutput',false);
else
    temp.Channel = cellfun(@(x) x(21:length(x)),temp.Label,'UniformOutput',false);
end
temp.Well = repmat(well,length(temp.Label),1);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
sitedata{r,c,f}=temp;
if(istable(welldata{r,c}))
    welldata{r,c} = [welldata{r,c};temp];
else
    welldata{r,c} = temp;
end


        end
        %% processing well summary
        temp1 = varfun(@mean,welldata{r,c},'GroupingVariables','Channel','InputVariables','Mean');
        channels = temp1.Channel;
        sizes = temp1.GroupCount;
        
        allchs = zeros(sizes(1),length(channels));
        
        for i=1:length(channels);
            allchs(:,i) = welldata{r,c}.Mean(strcmp(welldata{r,c}.Channel,channels{i}));
        end
        channels(i+1) = {'Area'};
        allchs(:,i+1) = welldata{r,c}.Area(strcmp(welldata{r,c}.Channel,channels{1}));
        DAPI = welldata{r,c}.Mean(strcmp(welldata{r,c}.Channel,'DAPI-0001'));
        AREA = allchs(:,i+1);
        channels(i+2) = {'IntDAPI'};
        allchs(:,i+2) = DAPI .* AREA;
        wellsum{r,c}=allchs;
        %tableTemp = array2table(allchs,'VariableNames',channels);
        clearvars allchs DAPI AREA temp1 sizes;
     end
end
