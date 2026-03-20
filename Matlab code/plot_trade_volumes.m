%------------------------------------------------------------
% Script:   plot_trade_volumes.m
% Purpose:  Plots export and import volumes by country
%           (emerging solid, advanced dashed) and prints
%           the values to console.
%
% Requires: Run linkage_model.m first — uses Trade_country_year,
%           N, T, S, year_labels, country_list from workspace.
%
% Output:   Volume Plot and Tables
%------------------------------------------------------------

%% === Figure: Trade Volumes - emerging (solid) vs advanced (dashed) ===

% --- Global style ---
set(groot, 'defaultAxesFontName',             'Times New Roman');
set(groot, 'defaultTextFontName',             'Times New Roman');
set(groot, 'defaultAxesFontSize',             12);
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

% --- Country groupings ---
emerging_list = {'BRA','CHN','IND','MEX','RUS'};
advanced_list    = {'CAN','DEU','ESP','FRA','GBR','ITA','JPN','KOR','USA'};
all_countries = [emerging_list, advanced_list];

% --- Colors ---
cmap = [
    0.00, 0.45, 0.70;
    0.00, 0.69, 0.31;
    0.85, 0.33, 0.10;
    0.49, 0.18, 0.56;
    0.93, 0.69, 0.13;
    0.30, 0.75, 0.93;
    0.64, 0.08, 0.18;
    0.47, 0.67, 0.19;
    1.00, 0.60, 0.78;
    0.00, 0.00, 0.55;
    0.75, 0.75, 0.00;
    0.50, 0.50, 0.50;
    0.87, 0.49, 0.00;
    0.00, 0.50, 0.50;
];

% --- Compute Exports and Imports ---
Exports = zeros(N, T);
Imports = zeros(N, T);

for t = 1:T
    X = Trade_country_year{t};
    for s = 1:S
        Xs = X(:,:,s);
        Xs(1:N+1:end) = 0;
        Exports(:,t) = Exports(:,t) + sum(Xs, 2);
        Imports(:,t) = Imports(:,t) + sum(Xs, 1)';
    end
end

% --- Figure with extra space at bottom for legend ---
fig = figure('Units','inches','Position',[1 1 13 6.5]);

subplot_data  = {Exports, Imports};
subplot_ylbl = {'Export volume ($10^5$ million USD)', 'Import volume ($10^6$ million USD)'};
subplot_title = {'Exports by country', 'Imports by country'};

h_lines = gobjects(numel(all_countries), 1);

for sp = 1:2
    ax = subplot(1, 2, sp);
    pos = get(ax, 'Position');
    set(ax, 'Position', [pos(1), pos(2)+0.18, pos(3), pos(4)-0.18]);
    hold(ax, 'on');

    for c = 1:numel(all_countries)
        cname = all_countries{c};
        idx   = find(strcmpi(country_list, cname));
        if isempty(idx), continue; end

        if ismember(cname, emerging_list)
            ls = '-';  mk = 'o';
        else
            ls = '--'; mk = 's';
        end

        h = plot(ax, year_labels, subplot_data{sp}(idx,:), ...
            'LineStyle',  ls,        ...
            'Marker',     mk,        ...
            'Color',      cmap(c,:), ...
            'LineWidth',  1.4,       ...
            'MarkerSize', 5,         ...
            'DisplayName', cname);

        if sp == 1
            h_lines(c) = h;
        end
    end

    xlabel(ax, 'Year');
    ylabel(ax, subplot_ylbl{sp});
    title(ax,  subplot_title{sp});
    xticks(ax, year_labels);
    xtickangle(ax, 45);
    grid(ax, 'on');
    hold(ax, 'off');
end

% --- Single shared legend: centered and closer to plots ---
lgd = legend(h_lines, all_countries, ...
    'Orientation', 'horizontal',     ...
    'NumColumns',   7,               ...
    'Units',        'normalized',    ...
    'FontSize',      9);

drawnow;  % force MATLAB to compute legend size before reading Position
lgd_width  = lgd.Position(3);
lgd_height = lgd.Position(4);
lgd.Position = [(1 - lgd_width) / 2, 0.06, lgd_width, lgd_height];

%% === Print plotted values to console ===

fprintf('\n=== EXPORT VOLUMES ===\n');
fprintf('%-6s', 'Year');
for c = 1:numel(all_countries)
    cname = all_countries{c};
    idx   = find(strcmpi(country_list, cname));
    if isempty(idx), continue; end
    fprintf('%12s', cname);
end
fprintf('\n');

for t = 1:T
    fprintf('%-6d', year_labels(t));
    for c = 1:numel(all_countries)
        cname = all_countries{c};
        idx   = find(strcmpi(country_list, cname));
        if isempty(idx), continue; end
        fprintf('%12.4f', Exports(idx, t));
    end
    fprintf('\n');
end

fprintf('\n=== IMPORT VOLUMES ===\n');
fprintf('%-6s', 'Year');
for c = 1:numel(all_countries)
    cname = all_countries{c};
    idx   = find(strcmpi(country_list, cname));
    if isempty(idx), continue; end
    fprintf('%12s', cname);
end
fprintf('\n');

for t = 1:T
    fprintf('%-6d', year_labels(t));
    for c = 1:numel(all_countries)
        cname = all_countries{c};
        idx   = find(strcmpi(country_list, cname));
        if isempty(idx), continue; end
        fprintf('%12.4f', Imports(idx, t));
    end
    fprintf('\n');
end