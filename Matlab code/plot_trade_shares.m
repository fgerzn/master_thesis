%------------------------------------------------------------
% Script:   plot_trade_shares.m
% Purpose:  Computes and plots intra-group vs inter-group trade
%           shares for emerging and advanced economies,
%           averaged over all years and excluding domestic trade.
%
% Requires: Run linkage_model.m first — uses Trade_country_year,
%           N, T, country_list from workspace.
%
% Output:   Trade shares plot and table
%------------------------------------------------------------

%% === Intra-group trade shares (excl. domestic, ROW) ===

% --- Global style ---
set(groot, 'defaultAxesFontName',             'Times New Roman');
set(groot, 'defaultTextFontName',             'Times New Roman');
set(groot, 'defaultAxesFontSize',             12);
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

% --- Average trade volume over years (sum over sectors) ---
trade_vol = zeros(N, N);
for t = 1:T
    X = Trade_country_year{t};
    trade_vol = trade_vol + sum(X, 3);
end
trade_vol = trade_vol / T;

% --- Exclude ROW ---
is_ROW       = strcmp(country_list, 'ROW');
keep         = ~is_ROW;
vol_filt     = trade_vol(keep, keep);
country_filt = country_list(keep);
N_filt       = sum(keep);

% --- Exclude domestic trade ---
vol_filt_nodiag = vol_filt;
vol_filt_nodiag(logical(eye(N_filt))) = 0;

% --- Group definitions ---
advanced    = {'CAN','FRA','ITA','JPN','GBR','USA','DEU','KOR','ESP'};
emerging = {'BRA','IND','CHN','RUS','MEX'};

is_advanced    = ismember(country_filt, advanced);
is_emerging = ismember(country_filt, emerging);

groups      = {is_emerging, is_advanced};
group_names = {'Emerging Economies', 'Advanced Economies'};
nG          = numel(groups);

% --- Compute shares ---
intra = zeros(nG, 1);
extra = zeros(nG, 1);

for g = 1:nG
    mask_g   = groups{g};
    total    = sum(vol_filt_nodiag(mask_g, :), 'all');
    within   = sum(vol_filt_nodiag(mask_g, mask_g), 'all');
    intra(g) = 100 * within / total;
    extra(g) = 100 - intra(g);
end

% --- Plot ---
fig = figure('Units', 'inches', 'Position', [1 1 6 3]);
ax  = axes(fig);

b = bar(ax, 1:nG, [intra, extra], 'stacked', 'BarWidth', 0.3);

b(1).FaceColor = [0.00, 0.45, 0.70];   % within-group  - blue
b(2).FaceColor = [0.80, 0.80, 0.80];   % outside-group - light grey

set(ax, ...
    'XTick',      1:nG,        ...
    'XTickLabel', group_names, ...
    'YLim',       [0, 100],    ...
    'Box',        'off');

ylabel(ax, 'Share of total trade (\%)');
grid(ax, 'on');

% --- Centered legend below plot ---
lgd = legend(ax, {'Within group', 'With other group'}, ...
    'Orientation', 'horizontal', ...
    'Units',       'normalized', ...
    'FontSize',     10);

drawnow;
lgd_width    = lgd.Position(3);
lgd_height   = lgd.Position(4);
lgd.Position = [(1 - lgd_width) / 2, 0.05, lgd_width, lgd_height];

% Shift axes up slightly to make room for legend
pos = get(ax, 'Position');
set(ax, 'Position', [pos(1), pos(2)+0.11, pos(3), pos(4)-0.11]);

exportgraphics(fig, 'intra_group_trade_shares.pdf', 'ContentType', 'vector');
fprintf('Saved intra_group_trade_shares.pdf\n');

% --- Print plotted values ---
fprintf('\n=== Intra-group trade shares (excl. domestic, avg. over years) ===\n');
fprintf('%-22s  %12s  %12s\n', 'Group', 'Within (%)', 'Outside (%)');
fprintf('%s\n', repmat('-', 1, 50));
for g = 1:nG
    fprintf('%-22s  %11.2f%%  %11.2f%%\n', group_names{g}, intra(g), extra(g));
end