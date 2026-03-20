%------------------------------------------------------------
% Script:   linkage_model.m
% Purpose:  Full Armington model with sectors and intermediate
%           input linkages. Solves for wage hats and price index
%           hats year-by-year using a rolling baseline, and plots
%           the emerging/advanced GDP ratio (data vs model).
%           Also includes the option to for scenario analysis
%
% Inputs:   DATA_*.mat              — trade data (one per year)
%           PARAMETERS_YY_AGG.mat  — structural parameters (one per year,
%                                    produced by clean_aggregate_WIOD.m)
%           country_list.csv        — country labels
%
% Output:   gdp_ratio_data_vs_model.pdf
%           Table printed to console: explained convergence (%)
%
% Requires: Optimization Toolbox (fsolve)
%------------------------------------------------------------

%% === clear all
clear; clc; close all;

%% === Settings
%------------------------------------------------------------
emerging_codes = {'BRA','IND','CHN','RUS','MEX'};
advanced_codes = {'CAN','FRA','ITA','JPN','GBR','USA','DEU','KOR','ESP'};
%------------------------------------------------------------

%% === Load DATA_*.mat in correct chronological order ===
files_raw = dir('DATA_*.mat');
nFiles = numel(files_raw);

years_true = zeros(nFiles,1);
for k = 1:nFiles
    yy = str2double(files_raw(k).name(6:7));
    years_true(k) = (yy >= 90) * (1900+yy) + (yy < 90) * (2000+yy);
end
[year_labels, sort_idx] = sort(years_true);
files = files_raw(sort_idx);

Trade_country_year = cell(nFiles,1);
for k = 1:nFiles
    data = load(files(k).name);
    Trade_country_year{k} = data.TRADEs;  % keep full sector structure
    fprintf("Loaded %s (%d)\n", files(k).name, year_labels(k));
end

country_list = readcell('country_list.csv');
N = numel(country_list);

%% === Parameters
load("PARAMETERS_95_AGG.mat")
DELTA = PARAMETERS.DELTA;        
beta  = PARAMETERS.BETA;           
GAMMA = PARAMETERS.GAMMA; 

%% === scenario option for a proportionally changed sigma
sigma = PARAMETERS.SIGMA(:);      % S x 1
%sigma = sigma * (8/ mean(sigma));

%% === Add tiny epsilon to all trade flows to avoid zeros ===
epsilon = 1e-6;  
for k = 1:nFiles
    Trade_country_year{k} = Trade_country_year{k} + epsilon;
end
disp('Added small epsilon to all trade flows to avoid zeros.');

%% === Trade cost changes               

N = size(Trade_country_year{1},1);
S = size(Trade_country_year{1},3);
T = nFiles;

tau = NaN(N,N,S,T);

for t = 1:T
    X = Trade_country_year{t};    % N x N x S
    for s = 1:S
        expo = 1/(2*(sigma(s)-1));
        Xii = diag(X(:,:,s));
        for i = 1:N
            for j = 1:N
                if X(i,j,s)>0 && X(j,i,s)>0 && Xii(i)>0 && Xii(j)>0
                    tau(i,j,s,t) = ((Xii(i)*Xii(j)) / ...
                                     (X(i,j,s)*X(j,i,s)))^expo;
                end
            end
        end
    end
end

% tau is N x N x S x T
tau_hat = NaN(size(tau,1), size(tau,2), size(tau,3), size(tau,4));

for t = 2:size(tau,4)
    tau_hat(:,:,:,t) = tau(:,:,:,t) ./ tau(:,:,:,t-1);
end

%% === Scenario Analysis for certain countries/groups/placebo
% placebo case
tau_hat(:,:,:,:)=1;

%   china case
% === Keep tau_hat only for links involving country 3; set all others to 1 ===
% === Set tau_hat to 1 only for links involving country 1 ===
%{
country_idx = 3;

mask = false(N,N);
mask(country_idx,:) = true;      % row 1
mask(:,country_idx) = true;      % column 1

for t = 2:size(tau_hat,4)
    for s = 1:size(tau_hat,3)
        tmp = tau_hat(:,:,s,t);
        tmp(mask) = 1;
        tau_hat(:,:,s,t) = tmp;
    end
end
%}
% (t=1 stays NaN by construction)         

%% === SOLVER SETUP (all years t=2..T) ===
% Solves for: w_hat (N x 1), P_hatI (N x S), and implied lambda_hat (N x N x S)
% Numéraire: mean(w_hat)=1 (implemented as sum(log w_hat)=0)

param_files = cell(nFiles,1);
for k = 1:nFiles
    yy = files(k).name(6:7);   % same suffix as DATA_YY...
    param_files{k} = ['PARAMETERS_' yy '_AGG.mat'];
end                       

% --- baseline lambda_ijs for each year from trade tensor X (exports i -> importer j) ---
lambda_base = NaN(N,N,S,T);
for t = 1:T
    X = Trade_country_year{t};
    for s = 1:S
        Ej = sum(X(:,:,s),1);                 % 1 x N (total expenditure of importer j)
        lambda_base(:,:,s,t) = X(:,:,s) ./ Ej; % divide each column j by Ej(j)
    end
end

% --- precompute exponents for price index equation (log form) ---
% exponent(i,t,s) = sum_a DELTA(i,t,a)*GAMMA(i,s,a)/(sigma(t)-1)
Expo = zeros(N,S,S);               % (i,s,t)
for i = 1:N
    for s = 1:S
        g = reshape(GAMMA(i,s,:),[],1);   % a x 1
        for tt = 1:S
            d = reshape(DELTA(i,tt,:),1,[]); % 1 x a
            Expo(i,s,tt) = (d * g) / (sigma(tt)-1);
        end
    end
end

% --- storage ---
w_hat_full  = NaN(N,T);      % full model (observed tau_hat)
P_hatI_all  = NaN(N,S,T);
lam_hat_all = NaN(N,N,S,T);

w_hat_full(:,1) = ones(N,1);


% --- initial guesses ---
x0_full = zeros((N-1) + N*S,1);  
opts = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-10,'StepTolerance',1e-10);

for t = 2:T
    P = load(param_files{t-1});
    phi = P.PARAMETERS.PHI;
    % === Trade-cost change t-1 -> t (FULL model)
    tauh = tau_hat(:,:,:,t);          % N x N x S

    % === Moving baseline shares (from actual t-1 data)
    L0 = lambda_base(:,:,:,t-1);      % N x N x S

    % === Moving baseline data (actual t-1)
    X_base = Trade_country_year{t-1}; % N x N x S

    % Exporter-sector totals
    Eexp_base = squeeze(sum(X_base,2));   % N x S

    % Baseline income (actual t-1)
    Y0 = sum(beta .* Eexp_base, 2);       % N x 1

    % Baseline gross revenue
    REV0 = squeeze(sum(sum(X_base,2),3)); % N x 1

    % Baseline absorption
    Eabs_base = squeeze(sum(sum(X_base,1),3))';  % N x 1

    % Trade imbalance
    Tprime = Eabs_base - REV0;            % N x 1

    f_full = @(x) equilibrium_residuals( ...
        x,N,S,beta,sigma,phi,L0,tauh,Y0,REV0,Tprime,Expo);

    xsol_full = fsolve(f_full,x0_full,opts);
    x0_full   = xsol_full;   

    [w_hat_t, P_hatI_t, lam_hat_t] = unpack_and_imply( ...
        xsol_full,N,S,beta,sigma,L0,tauh,Expo);

    w_hat_full(:,t)      = w_hat_t;
    P_hatI_all(:,:,t)    = P_hatI_t;
    lam_hat_all(:,:,:,t) = lam_hat_t;

end

%% === Figure 1: emerging vs advanced GDP Ratio (Data vs Model) ===

% --- wage hats: set first year = 1 and cumulate ---
w_cum = ones(N,T);
for t = 2:T
    w_cum(:,t) = w_cum(:,t-1) .* w_hat_full(:,t);
end

% --- "GDP" using picture formula: Y_i = sum_s beta_is * sum_j X_ijs ---
GDP_data  = zeros(N,T);
GDP_model = zeros(N,T);

% baseline year-1 Y_init from formula (not simple sum of trade flows)
X1 = Trade_country_year{1};               % N x N x S
E1 = squeeze(sum(X1,2));                  % N x S   (sum over importers j)
Y_init = sum(beta .* E1, 2);              % N x 1

for t = 1:T
    Xt = Trade_country_year{t};           % N x N x S
    Et = squeeze(sum(Xt,2));              % N x S
    GDP_data(:,t) = sum(beta .* Et, 2);   % N x 1  (picture formula)

    GDP_model(:,t) = w_cum(:,t) .* Y_init;
end

is_emerging = ismember(country_list, emerging_codes);
is_advanced    = ismember(country_list, advanced_codes);

ratio_data  = sum(GDP_data(is_emerging,:),1) ./ sum(GDP_data(is_advanced,:),1);
ratio_model = sum(GDP_model(is_emerging,:),1) ./ sum(GDP_model(is_advanced,:),1);

% --- Global style ---
set(groot, 'defaultAxesFontName',             'Times New Roman');
set(groot, 'defaultTextFontName',             'Times New Roman');
set(groot, 'defaultAxesFontSize',             12);
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

fig = figure('Units', 'inches', 'Position', [1 1 6 4.5]);
ax  = axes(fig);
pos = get(ax, 'Position');
set(ax, 'Position', [pos(1), pos(2)+0.10, pos(3), pos(4)-0.10]);
hold(ax, 'on');
plot(ax, year_labels, ratio_data,  'o-',  'LineWidth', 1.4, 'Color', [0.00, 0.45, 0.70], 'DisplayName', 'Data');
plot(ax, year_labels, ratio_model, 's--', 'LineWidth', 1.4, 'Color', [0.85, 0.33, 0.10], 'DisplayName', 'Linkage Model');
hold(ax, 'off');
xlabel(ax, 'Year'); ylabel(ax, 'Emerging / Advanced Economies GDP Ratio');
xticks(ax, year_labels); xtickangle(ax, 45); grid(ax, 'on');
lgd = legend(ax, 'Orientation', 'horizontal', 'Units', 'normalized', 'FontSize', 10);
drawnow;
lgd_width = lgd.Position(3); lgd_height = lgd.Position(4);
lgd.Position = [(1 - lgd_width)/2, 0.01, lgd_width, lgd_height];

%% === Table: % of DATA change explained by MODEL (emerging/advanced ratio) ===
% Explained(t) = 100 * (ratio_model(t)-ratio_model(1)) / (ratio_data(t)-ratio_data(1))

dData  = ratio_data  - ratio_data(1);
dModel = ratio_model - ratio_model(1);

explained_pct = 100 * (dModel ./ dData);

% baseline year is 0/0 -> undefined
explained_pct(1) = NaN;

% also guard against any year where data change is zero
explained_pct(dData == 0) = NaN;

% Put years as columns
yearCols = matlab.lang.makeValidName(string(year_labels));

explainedTable = array2table(explained_pct, ...
    'RowNames', "Explained % (emerging/advanced ratio change)", ...
    'VariableNames', cellstr(yearCols));

disp("=== Explained-by-model (percent of data change explained) ===");
disp(explainedTable);

%% === local functions ===
function F = equilibrium_residuals(x,N,S,beta,sigma,phi,L0,tauh,Y,REV0,Tprime,Expo)
    [w_hat, P_hatI, lam_hat] = unpack_and_imply(x,N,S,beta,sigma,L0,tauh,Expo);

    % (1) Wage / income equation: N equations, but drop one due to numéraire
    RHS = zeros(N,1);
    for s = 1:S
        % denom for each importer j
        denom = zeros(1,N);
        for j = 1:N
            tmp = 0;
            for m = 1:N
                cm = (w_hat(m)^beta(m,s)) * (P_hatI(m,s)^(1-beta(m,s)));
                tmp = tmp + L0(m,j,s) * (tauh(m,j,s)*cm)^(1-sigma(s));
            end
            denom(j) = tmp;
        end

        for i = 1:N
            for j = 1:N
                ci = (w_hat(i)^beta(i,s)) * (P_hatI(i,s)^(1-beta(i,s)));
                lamhat_ijs = (tauh(i,j,s)*ci)^(1-sigma(s)) / denom(j);
                RHS(i) = RHS(i) + beta(i,s) * (L0(i,j,s)*lamhat_ijs) * phi(j,s) * (w_hat(j)*REV0(j) + Tprime(j));
            end
        end
    end

    wage_res = w_hat .* Y - RHS;   % Nx1
    wage_res = wage_res(1:N-1);    % enforce mean(w_hat)=1 via parametrization

    % (2) Price index equation in logs: log P_hatI(i,s) - [log w_hat(i) + sum_t Expo(i,s,t)*log hatlambda_dom(i,t)]
    hat_dom = zeros(N,S);
    for tt = 1:S
        hat_dom(:,tt) = diag(lam_hat(:,:,tt)); % hat lambda_{ii,tt}
    end
    logP_rhs = log(w_hat) + squeeze(sum(Expo .* reshape(log(hat_dom),N,1,S), 3)); % NxS
    price_res = log(P_hatI) - logP_rhs;  % NxS

    F = [wage_res; price_res(:)];
end

function [w_hat, P_hatI, lam_hat] = unpack_and_imply(x,N,S,beta,sigma,L0,tauh,Expo)
    % unpack with geometric-mean numéraire: sum(log w_hat) = 0
lw = zeros(N,1);
lw(1:N-1) = x(1:N-1);
lw(N) = -sum(lw(1:N-1));
w_hat = exp(lw);
    lp = x((N-1)+1:end);
    P_hatI = reshape(exp(lp),N,S);

    % implied lambda_hat from equation (3)
    lam_hat = NaN(N,N,S);
    for s = 1:S
        denom = zeros(1,N);
        for j = 1:N
            tmp = 0;
            for m = 1:N
                cm = (w_hat(m)^beta(m,s)) * (P_hatI(m,s)^(1-beta(m,s)));
                tmp = tmp + L0(m,j,s) * (tauh(m,j,s)*cm)^(1-sigma(s));
            end
            denom(j) = tmp;
        end
        for i = 1:N
            ci = (w_hat(i)^beta(i,s)) * (P_hatI(i,s)^(1-beta(i,s)));
            for j = 1:N
                lam_hat(i,j,s) = (tauh(i,j,s)*ci)^(1-sigma(s)) / denom(j);
            end
        end
    end
end



