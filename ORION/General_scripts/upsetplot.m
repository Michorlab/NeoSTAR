%% upset_from_tables.m
% UpSet plot generator for MATLAB
% Inputs: one or more tables
%   Format A (wide/binary): columns are sets (0/1 or logical). Optional ID column.
%   Format B (long): two columns: ID and Set (categorical/string/char).
%
% Example (wide):
%   T = table((rand(30,1)>0.5),(rand(30,1)>0.5),(rand(30,1)>0.5), ...
%             'VariableNames',{'SetA','SetB','SetC'});
%   upsetplot(T);
%
% Example (long):
%   TL = table([1;1;2;3;3],[ "SetA";"SetB";"SetA";"SetB";"SetC"], ...
%              'VariableNames',{'ID','Set'});
%   upsetplot(TL,'InputFormat','long');
%
% Author: ChatGPT
% -------------------------------------------------------------------------

function out = upsetplot(varargin)
% out = upsetplot(T1, T2, ..., 'Name',Value,...)
%
% Name-Value options:
%   'InputFormat'  : 'auto' (default) | 'wide' | 'long'
%   'IDVar'         : string/char for ID column (wide only), default '' (auto)
%   'LongIDVar'     : default 'ID'   (long only)
%   'LongSetVar'    : default 'Set'  (long only)
%   'SetVars'       : cellstr/string array of set column names (wide only) (default auto)
%   'SetOrder'      : 'desc' (by set size) | 'asc' | 'input' | string array of names
%   'TopN'          : number of intersections to show (default 20)
%   'MinIntersection': minimum intersection size to keep (default 1)
%   'IncludeEmpty'  : false (default) | true (include empty intersections)
%   'ShowSetSizes'  : true (default) | false (left bars)
%   'SortIntersections' : 'desc' (default) | 'asc' | 'none'
%   'Title'         : char/string (default '')
%
% Returns struct:
%   out.Sets, out.Intersections, out.MembershipMatrix, out.Fig, out.Axes, out.Table

% ---------------- Parse inputs ----------------
tables = {};
k = 1;
while k <= nargin && istable(varargin{k})
    tables{end+1} = varargin{k}; %#ok<AGROW>
    k = k + 1;
end
if isempty(tables)
    error('Provide at least one input table.');
end

p = inputParser;
p.addParameter('InputFormat','auto',@(x)ischar(x)||isstring(x));
p.addParameter('IDVar','',@(x)ischar(x)||isstring(x));
p.addParameter('LongIDVar','ID',@(x)ischar(x)||isstring(x));
p.addParameter('LongSetVar','Set',@(x)ischar(x)||isstring(x));
p.addParameter('SetVars',[],@(x)iscellstr(x)||isstring(x)||isempty(x));
p.addParameter('SetOrder','desc',@(x)ischar(x)||isstring(x)||iscellstr(x)||isstring(x));
p.addParameter('TopN',20,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('MinIntersection',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('IncludeEmpty',false,@islogical);
p.addParameter('ShowSetSizes',true,@islogical);
p.addParameter('SortIntersections','desc',@(x)ischar(x)||isstring(x));
p.addParameter('Title','',@(x)ischar(x)||isstring(x));
p.parse(varargin{k:end});
opt = p.Results;

% ---------------- Normalize to membership matrix ----------------
fmt = lower(string(opt.InputFormat));
if fmt=="auto"
    fmt = inferFormat(tables{1}, string(opt.LongIDVar), string(opt.LongSetVar));
end

switch fmt
    case "wide"
        [M, setNames, rowIDs] = wideTablesToMatrix(tables, string(opt.IDVar), opt.SetVars);
    case "long"
        [M, setNames, rowIDs] = longTablesToMatrix(tables, string(opt.LongIDVar), string(opt.LongSetVar));
    otherwise
        error('InputFormat must be auto|wide|long.');
end

% Ensure logical
M = M ~= 0;

% ---------------- Order sets ----------------
setSizes = sum(M,1);
setNames = string(setNames);

if isstring(opt.SetOrder) || ischar(opt.SetOrder)
    so = lower(string(opt.SetOrder));
    if so=="desc"
        [~,idx] = sort(setSizes,'descend');
    elseif so=="asc"
        [~,idx] = sort(setSizes,'ascend');
    elseif so=="input"
        idx = 1:numel(setNames);
    else
        % user provided a single string not matching keywords -> treat as error
        error('Unknown SetOrder. Use desc|asc|input or provide a list of set names.');
    end
else
    % list of set names
    desired = string(opt.SetOrder);
    [tf,idx] = ismember(desired,setNames);
    if ~all(tf)
        missing = desired(~tf);
        error('SetOrder includes unknown set(s): %s', strjoin(missing,", "));
    end
end
M = M(:,idx);
setNames = setNames(idx);
setSizes = setSizes(idx);

% ---------------- Compute intersections ----------------
[nRows,nSets] = size(M);
% represent each row membership as binary code
codes = uint64(0);
for j = 1:nSets
    codes = codes + uint64(M(:,j)) .* bitshift(uint64(1), j-1);
end

% Optionally remove empty intersection (code==0)
if ~opt.IncludeEmpty
    keepRows = codes ~= 0;
    codesK = codes(keepRows);
else
    codesK = codes;
end

% Count unique intersections
[uCodes,~,ic] = unique(codesK,'stable');
counts = accumarray(ic, 1);

% Filter by min size
keepI = counts >= opt.MinIntersection;
uCodes = uCodes(keepI);
counts = counts(keepI);

% Sort intersections
sortMode = lower(string(opt.SortIntersections));
if sortMode=="desc"
    [counts,ord] = sort(counts,'descend');
    uCodes = uCodes(ord);
elseif sortMode=="asc"
    [counts,ord] = sort(counts,'ascend');
    uCodes = uCodes(ord);
elseif sortMode=="none"
    % keep as-is
else
    error('SortIntersections must be desc|asc|none.');
end

% Take TopN
topN = min(opt.TopN, numel(uCodes));
uCodes = uCodes(1:topN);
counts = counts(1:topN);

% Build membership matrix for intersections (topN x nSets)
IM = false(topN, nSets);
for i = 1:topN
    for j = 1:nSets
        IM(i,j) = bitand(uCodes(i), bitshift(uint64(1), j-1)) ~= 0;
    end
end

% ---------------- Plot ----------------
fig = figure('Color','w','Name','UpSet Plot');
tlo = tiledlayout(fig, 2, 2, 'TileSpacing','compact','Padding','compact');

% Axes layout:
% (1,1) optional set-size bars (left) spanning rows 2? We'll use 2x2:
% Top-right: intersection bars
% Bottom-right: membership matrix
% Bottom-left: set-size bars (horizontal)
axTop   = nexttile(tlo, 2);   % (row1,col2)
axMat   = nexttile(tlo, 4);   % (row2,col2)
axLeft  = nexttile(tlo, 3);   % (row2,col1)
axDummy = nexttile(tlo, 1);   % (row1,col1) placeholder for title etc.
axis(axDummy,'off');

% Title
if strlength(string(opt.Title))>0
    title(axDummy, string(opt.Title), 'FontWeight','bold', 'FontSize',12);
else
    title(axDummy, 'UpSet plot', 'FontWeight','bold', 'FontSize',12);
end

% Intersection bars (top)
bar(axTop, 1:topN, counts, 'FaceAlpha',1);
axTop.XLim = [0.5, topN+0.5];
axTop.XTick = [];
ylabel(axTop,'Intersection size');
grid(axTop,'on');
box(axTop,'off');

% Membership matrix (bottom-right)
% Coordinates: x = intersection index, y = set index
axes(axMat); %#ok<LAXES>
cla(axMat);
hold(axMat,'on');

% background faint dots (all positions)
[xg,yg] = meshgrid(1:topN, 1:nSets);
scatter(axMat, xg(:), yg(:), 20, 'filled', 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0);

% filled dots for membership
[idxX, idxY] = find(IM'); % IM' is nSets x topN
scatter(axMat, idxY, idxX, 32, 'filled');

% connect vertical lines for each intersection where >=2 sets
for i = 1:topN
    ys = find(IM(i,:));
    if numel(ys) >= 2
        plot(axMat, [i i], [min(ys) max(ys)], 'LineWidth', 1.5);
    end
end

axMat.YLim = [0.5, nSets+0.5];
axMat.XLim = [0.5, topN+0.5];
axMat.YDir = 'reverse';
axMat.YTick = 1:nSets;
axMat.YTickLabel = setNames;
axMat.XTick = 1:topN;
axMat.XTickLabel = repmat({''},1,topN);
xlabel(axMat,'Intersections (top N)');
grid(axMat,'off');
box(axMat,'off');

% Set size bars (bottom-left)
if opt.ShowSetSizes
    barh(axLeft, 1:nSets, setSizes);
    axLeft.YDir = 'reverse';
    axLeft.YTick = 1:nSets;
    axLeft.YTickLabel = setNames;
    xlabel(axLeft,'Set size');
    grid(axLeft,'on');
    box(axLeft,'off');
else
    axis(axLeft,'off');
end

% Align Y tick labels between left and matrix if showing set sizes
if opt.ShowSetSizes
    axMat.YTickLabel = setNames;
end

% ---------------- Outputs ----------------
out = struct();
out.Sets = table(setNames(:), setSizes(:), 'VariableNames', {'Set','Size'});
out.Intersections = table(uCodes(:), counts(:), 'VariableNames', {'Code','Size'});
out.MembershipMatrix = IM;
out.RowIDs = rowIDs;
out.Fig = fig;
out.Axes = struct('Top',axTop,'Matrix',axMat,'SetSizes',axLeft);
out.Table = tables;

end

%% -------- helper: infer format --------
function fmt = inferFormat(T, longID, longSet)
vars = string(T.Properties.VariableNames);
if any(vars==longID) && any(vars==longSet)
    fmt = "long";
else
    % if many logical/numeric 0/1 columns, assume wide
    fmt = "wide";
end
end

%% -------- helper: wide tables to matrix --------
function [M, setNames, rowIDs] = wideTablesToMatrix(tables, idVar, setVarsOpt)
% Concatenate rows; sets are union of set columns across all tables
allVars = string.empty;
for i = 1:numel(tables)
    allVars = union(allVars, string(tables{i}.Properties.VariableNames), 'stable');
end

% Determine ID column (optional)
if strlength(idVar)==0
    % auto-detect common ID-ish names
    candidates = ["ID","Id","Sample","Cell","Name"];
    idVar = "";
    for c = candidates
        if any(allVars==c), idVar = c; break; end
    end
end

% Determine set variables
if isempty(setVarsOpt)
    % choose columns that are logical or numeric and not the ID var
    setNames = string.empty;
    for i = 1:numel(tables)
        T = tables{i};
        v = string(T.Properties.VariableNames);
        for j = 1:numel(v)
            if v(j)==idVar, continue; end
            col = T.(v(j));
            if islogical(col) || (isnumeric(col) && all(ismember(unique(col(~isnan(col))), [0 1])))
                setNames(end+1) = v(j); %#ok<AGROW>
            end
        end
    end
    setNames = unique(setNames,'stable');
else
    setNames = string(setVarsOpt);
end

% Build combined membership matrix
M = false(0, numel(setNames));
rowIDs = [];
for i = 1:numel(tables)
    T = tables{i};
    n = height(T);
    Mi = false(n, numel(setNames));
    for j = 1:numel(setNames)
        nm = setNames(j);
        if any(string(T.Properties.VariableNames)==nm)
            col = T.(nm);
            Mi(:,j) = col ~= 0;
        else
            Mi(:,j) = false(n,1);
        end
    end
    M = [M; Mi]; %#ok<AGROW>
    if strlength(idVar)>0 && any(string(T.Properties.VariableNames)==idVar)
        rowIDs = [rowIDs; T.(idVar)]; %#ok<AGROW>
    else
        rowIDs = [rowIDs; (1:n)' + (size(M,1)-n)]; %#ok<AGROW>
    end
end
end

%% -------- helper: long tables to matrix --------
function [M, setNames, rowIDs] = longTablesToMatrix(tables, idVar, setVar)
% Combine into single long table
TL = vertcat(tables{:});
vars = string(TL.Properties.VariableNames);
if ~any(vars==idVar) || ~any(vars==setVar)
    error('Long format requires columns %s and %s.', idVar, setVar);
end

ids = TL.(idVar);
sets = TL.(setVar);

% Normalize types
if iscell(ids), ids = string(ids); end
if iscell(sets), sets = string(sets); end
ids = string(ids);
sets = string(sets);

rowIDs = unique(ids,'stable');
setNames = unique(sets,'stable');

nRows = numel(rowIDs);
nSets = numel(setNames);
M = false(nRows, nSets);

[~,rid] = ismember(ids, rowIDs);
[~,sid] = ismember(sets, setNames);

for k = 1:numel(rid)
    if rid(k)>0 && sid(k)>0
        M(rid(k), sid(k)) = true;
    end
end
end
