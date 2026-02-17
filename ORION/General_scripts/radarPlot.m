function ax = radarPlot(X, labels, varargin)
%RADARPLOT  Radar (spider) plot, toolbox-free.
%
%   ax = radarPlot(X, labels)
%   ax = radarPlot(X, labels, 'Name', Value, ...)
%
% Inputs
%   X      : nSeries x nAxes numeric (or 1 x nAxes vector)
%   labels : 1 x nAxes cellstr or string array
%
% Name-Value options
%   'Parent'        : axes handle (default: gca)
%   'Limits'        : [min max] or 2 x nAxes (default: auto from data)
%   'GridRings'     : number of grid rings (default: 5)
%   'ShowGrid'      : true/false (default: true)
%   'ShowSpokes'    : true/false (default: true)
%   'LabelRadius'   : radius for labels (default: 1.12)
%   'LineWidth'     : line width for series (default: 2)
%   'Marker'        : marker char e.g. 'o' (default: 'o')
%   'Fill'          : true/false fill polygons (default: false)
%   'FillAlpha'     : 0..1 (default: 0.10)
%   'ThetaZeroTop'  : true/false set 0° at top (default: true)
%   'Clockwise'     : true/false (default: true)
%
% Notes
% - Colors are not set explicitly (MATLAB default ColorOrder).
% - If you need per-axis limits, pass 'Limits' as 2 x nAxes.

% ---------- Parse inputs ----------
if isvector(X)
    X = X(:).'; % force row vector
end
if ~isnumeric(X) || isempty(X)
    error('radarPlot:BadX','X must be a non-empty numeric array.');
end

nSeries = size(X,1);
nAxes   = size(X,2);

if isstring(labels), labels = cellstr(labels); end
if ischar(labels),   labels = cellstr(labels); end %#ok<CHARTEN>
if ~iscell(labels) || numel(labels) ~= nAxes
    error('radarPlot:LabelMismatch','labels length must equal size(X,2).');
end

p = inputParser;
p.FunctionName = 'radarPlot';

addParameter(p,'Parent',[],@(h) isempty(h) || isgraphics(h,'axes'));
addParameter(p,'Limits',[],@(v) isempty(v) || (isnumeric(v) && (isequal(size(v),[1 2]) || isequal(size(v),[2 nAxes]))));
addParameter(p,'GridRings',5,@(v) isnumeric(v)&&isscalar(v)&&v>=1);
addParameter(p,'ShowGrid',true,@(v) islogical(v)&&isscalar(v));
addParameter(p,'ShowSpokes',true,@(v) islogical(v)&&isscalar(v));
addParameter(p,'LabelRadius',1.12,@(v) isnumeric(v)&&isscalar(v)&&v>0);
addParameter(p,'LineWidth',2,@(v) isnumeric(v)&&isscalar(v)&&v>0);
addParameter(p,'Marker','o',@(s) ischar(s) || isstring(s));
addParameter(p,'Fill',false,@(v) islogical(v)&&isscalar(v));
addParameter(p,'FillAlpha',0.10,@(v) isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
addParameter(p,'ThetaZeroTop',true,@(v) islogical(v)&&isscalar(v));
addParameter(p,'Clockwise',true,@(v) islogical(v)&&isscalar(v));

parse(p,varargin{:});
opt = p.Results;

% ---------- Axes ----------
if isempty(opt.Parent)
    ax = gca;
else
    ax = opt.Parent;
end
axes(ax); %#ok<LAXES>
holdState = ishold(ax);
cla(ax);
hold(ax,'on');

axis(ax,'equal');
axis(ax,'off');

% ---------- Angles ----------
theta = linspace(0, 2*pi, nAxes+1);
theta(end) = []; % nAxes points

% Rotate so first axis is at top if requested
if opt.ThetaZeroTop
    theta = theta + pi/2;
end
% Make clockwise if requested (flip direction)
if opt.Clockwise
    theta = fliplr(theta);
end

thetaC = [theta theta(1)]; % closed

% ---------- Limits & scaling ----------
lims = opt.Limits;
if isempty(lims)
    % auto limits from data (ignore NaNs)
    mn = min(X(:), [], 'omitnan');
    mx = max(X(:), [], 'omitnan');
    if ~isfinite(mn) || ~isfinite(mx) || mn == mx
        mn = 0; mx = 1;
    end
    lims = [mn mx];
end

if isequal(size(lims), [1 2])
    % global limits
    minV = lims(1);
    maxV = lims(2);
    if maxV == minV, maxV = minV + 1; end
    scale = @(v) (v - minV) ./ (maxV - minV);
    Xs = scale(X);
else
    % per-axis limits: 2 x nAxes
    minV = lims(1,:);
    maxV = lims(2,:);
    denom = maxV - minV;
    denom(denom==0) = 1;
    Xs = (X - minV) ./ denom;
end
Xs = max(0, min(1, Xs)); % clamp to [0,1]

% ---------- Grid ----------
if opt.ShowGrid
    rings = opt.GridRings;
    for r = (1:rings)/rings
        xg = r*cos(thetaC);
        yg = r*sin(thetaC);
        plot(ax, xg, yg, ':', 'HandleVisibility','off');
    end
end

% ---------- Spokes ----------
if opt.ShowSpokes
    for k = 1:nAxes
        plot(ax, [0 cos(theta(k))], [0 sin(theta(k))], ':', 'HandleVisibility','off');
    end
end

% ---------- Labels ----------
for k = 1:nAxes
    tx = opt.LabelRadius*cos(theta(k));
    ty = opt.LabelRadius*sin(theta(k));
    text(ax, tx, ty, labels{k}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',7);
end

% ---------- Plot series ----------
for i = 1:nSeries
    rho = Xs(i,:);
    rhoC = [rho rho(1)];

    x = rhoC .* cos(thetaC);
    y = rhoC .* sin(thetaC);

    h = plot(ax, x, y, ...
        'LineWidth', opt.LineWidth, ...
        'Marker', char(opt.Marker));

    if opt.Fill
        % Use the line color for fill (without specifying a custom color)
        fc = h.Color;
        patch(ax, x, y, fc, 'FaceAlpha', opt.FillAlpha, 'EdgeColor','none', ...
            'HandleVisibility','off');
        % Keep line on top
        uistack(h,'top');
    end
end

% Restore hold state
if ~holdState
    hold(ax,'off');
end

end
