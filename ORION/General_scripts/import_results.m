function PIP_results = import_results(filename, startRow, endRow)
%% import_results files as table (from ImageJ dataset)
% Jerry 20160401


%% Initialize variables.
%mypath = uigetdir('F:\20160309-PIP-10x');


delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%filename = strcat(mypath,'\Result-',filename,'-CFP.xls');
display(filename);



formatSpec = '%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
PIP_results = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','Label','Area','Mean','StdDev','Min','Max','X','Y','XM','YM','Perim','BX','BY','Width','Height','Major','Minor','Angle','Circ','Feret','IntDen','Median','Skew','Kurt','VarName26','RawIntDen','Slice','FeretX','FeretY','FeretAngle','MinFeret','AR','Round','Solidity'});

