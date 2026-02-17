% Define input file
inputFile = uigetfile('*.txt');

% Read all lines
fileLines = readlines(inputFile);

% Regular expression to find LSP IDs
expr = 'LSP\d{5}';

% Initialize cell arrays
ids = {};
paths = {};

% Extract IDs and file paths
for i = 1:length(fileLines)
    line = fileLines(i);
    tokens = regexp(line, expr, 'match');
    if ~isempty(tokens)
        ids{end+1,1} = tokens{1};
        paths{end+1,1} = line;
    end
end

% Create table
T = table(ids, paths, 'VariableNames', {'ID', 'FileLocation'});

% Sort table by ID
T = sortrows(T, 'ID');

% Compute duplicate flags
[uniqueIDs, ~, idx] = unique(T.ID, 'stable');
duplicateCounts = accumarray(idx, 1);
isDuplicate = duplicateCounts(idx) > 1;

% Add column
T.IsDuplicate = isDuplicate;

% Preview
disp(head(T));
outfilename = strcat('parsed_',inputFile,'.csv');
% Save result
writetable(T, outfilename);