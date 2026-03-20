%------------------------------------------------------------
% Script:   baseline_model.m
% Purpose:  One-sector, no-intermediates Armington model.
%           Aggregates trade across sectors, computes tau and
%           hat algebra wage changes year-by-year, and plots
%           the emerging/advanced GDP ratio (data vs model).
%
% Inputs:   DATA_*.mat         — trade data (one file per year,
%                                produced by clean_aggregate_WIOD.m)
%           country_list.csv   — country labels
%
% Output:   Figure: emerging/advanced GDP ratio (data vs model)
%           Table printed to console: explained convergence (%)
%           
% Requires: Optimization Toolbox (fsolve)
%------------------------------------------------------------
%% === clear all
clear; clc; close all

%%
% SETTINGS
%------------------------------------------------------------
sigma     = 3.94;        % uniform trade elasticity (trade-weighted avg of sigma_s)
base_year = 1995;        % baseline year
emerging_codes = {'BRA','IND','CHN','RUS','MEX'};
advanced_codes    = {'CAN','FRA','ITA','JPN','GBR','USA','DEU','KOR','ESP'};
%%------------------------------------------------------------

%% === Chronological ordering: read DATA_*.mat sorted as 1995–1999, 2000–2007 ===

files_raw = dir('DATA_*.mat');
nFiles    = numel(files_raw);

years_true = zeros(nFiles,1);

for k = 1:nFiles
    fname = files_raw(k).name;
    yy    = str2double(fname(6:7));

    if yy >= 90
        years_true(k) = 1900 + yy;   % 95→1995, 99→1999
    else
        years_true(k) = 2000 + yy;   % 00→2000 ... 07→2007
    end
end

[year_labels, sort_idx] = sort(years_true);
files = files_raw(sort_idx);

%% === Load trade matrices in sorted order ===
Trade_country_year = cell(nFiles,1);

for k = 1:nFiles
    fname = files(k).name;
    data  = load(fname);
    TRADEs = data.TRADEs;

    % Aggregate over sectors (1-sector world)
    Trade_country_year{k} = sum(TRADEs,3);

    fprintf('Loaded %s (year = %d)\n', fname, year_labels(k));
end

country_list = readcell('country_list.csv');
n = numel(country_list);


%% === Baseline from data ===
base_idx  = find(year_labels == base_year,1);

if isempty(base_idx)
    error('Baseline year 1995 not found.');
end

X_base = Trade_country_year{base_idx};

Y_base_95 = sum(X_base,2);            % GDP

%% === Compute τ(t) for every year ===
tau_year = cell(nFiles,1);

for k = 1:nFiles
    Xk = Trade_country_year{k};
    tau = zeros(n);

    for i = 1:n
        for j = 1:n
            tau(i,j) = ((Xk(i,i)*Xk(j,j))/(Xk(i,j)*Xk(j,i)))^(1/(2*(sigma-1)));
        end
    end

    tau_year{k} = tau;
end

%% === Containers for step hats ===
w_hat_all      = cell(nFiles,1);
lambda_hat_all = cell(nFiles,1);
tau_hat_all    = cell(nFiles,1);
welfare_hat_all = cell(nFiles,1);
w_cum = nan(n, nFiles);
w_cum(:, base_idx) = ones(n,1);

%% === Dynamic hat algebra loop ===
for k = 1:nFiles

    if k == base_idx
        % Baseline year: hats = 1
        tau_hat_all{k}     = ones(n);
        w_hat_all{k}       = ones(n,1);
        lambda_hat_all{k}  = ones(n);
        welfare_hat_all{k} = ones(n,1);

        fprintf('Year %d: baseline (hats = 1)\n', year_labels(k));
        continue;
    end

    prev = k - 1;

    % DATA baseline at t-1 
X_prev = Trade_country_year{prev};

Y_base = sum(X_prev,2);          % baseline GDP/revenue from data (Nx1)
E_prev = sum(X_prev,1);          % baseline expenditure from data (1xN)
lambda_base = X_prev ./ E_prev;  % baseline trade shares from data (NxN)

T_base_step = E_prev' - Y_base;  % baseline transfers from data (Nx1)
    % τ̂_t = τ(t)/τ(t-1) (data-based)
    tau_hat = tau_year{k} ./ tau_year{prev};
    tau_hat_all{k} = tau_hat;

    % Solve for step wages ŵ_t = w(t)/w(t-1)
    w_hat_step = solve_w_hat_mean_numeraire( ...
                        lambda_base, tau_hat, sigma, Y_base, T_base_step);
    w_hat_all{k} = w_hat_step;
    w_cum(:,k) = w_cum(:,prev) .* w_hat_step;
    % Compute step λ̂_t
    C   = (tau_hat .* (w_hat_step .* ones(1,n))).^(1 - sigma);
    den = sum(lambda_base .* C, 1);
    lambda_hat_step = C ./ den;
    lambda_hat_all{k} = lambda_hat_step;

    fprintf('Processed year %d using DATA baseline.\n', year_labels(k));
end

%% === GDP (data) ===
GDP_data = nan(n, nFiles);
for k = 1:nFiles
    GDP_data(:,k) = sum(Trade_country_year{k},2);
end

%% === GDP (model) ===
% Counterfactual GDP: chain cumulated wage hats from 1995 baseline
GDP_model = w_cum .* Y_base_95;   

%% === emerging vs advanced GDP ratio ===

is_emerging = ismember(country_list, emerging_codes);
is_advanced    = ismember(country_list, advanced_codes);

ratioGDP_data  = nan(nFiles,1);
ratioGDP_model = nan(nFiles,1);

for k = 1:nFiles
    ratioGDP_data(k)  = sum(GDP_data(is_emerging,k))  / sum(GDP_data(is_advanced,k));
    ratioGDP_model(k) = sum(GDP_model(is_emerging,k)) / sum(GDP_model(is_advanced,k));
end

%% === Plot: emerging/advanced GDP ratio ===
[years_sorted, idx] = sort(year_labels);

figure;
plot(years_sorted, ratioGDP_data(idx), 'o-','LineWidth',1.4); hold on;
plot(years_sorted, ratioGDP_model(idx),'s--','LineWidth',1.4);
xlabel('Year');
ylabel('emerging GDP / advanced GDP');
title('emerging vs advanced: GDP Ratio (Data vs Model)');
legend('Data','Model');
grid on;

%% === Explained convergence (equation 4 in paper) ===
% Explained(t) = 100 * (ratio_model(t) - ratio_model(1)) / 
%                      (ratio_data(t)  - ratio_data(1))
% A value of x% means the model accounts for x% of the observed
% change in the emerging/advanced GDP ratio relative to 1995.

dData  = ratioGDP_data  - ratioGDP_data(base_idx);
dModel = ratioGDP_model - ratioGDP_model(base_idx);

explained_pct = 100 * (dModel ./ dData);

% Baseline year is 0/0 by construction
explained_pct(base_idx) = NaN;

% Guard against years where data change is zero
explained_pct(dData == 0) = NaN;

% Display as a table
yearCols = matlab.lang.makeValidName(string(year_labels'));
explainedTable = array2table(explained_pct', ...
    'VariableNames', cellstr(yearCols), ...
    'RowNames',      {'Explained (%)'});

disp('=== Explained convergence: baseline model ===');
disp(explainedTable);
%% === Helper functions (mean-wage numéraire + analytical Jacobian) ===
function w_hat = solve_w_hat_mean_numeraire(lambda, tau_hat, sigma, Y, T_base_step)
    % Numéraire: mean(w_hat) = 1 (geometric-mean normalisation)
    % Returns N×1 vector w_hat. Uses analytical Jacobian for speed/stability.
 
    N   = length(Y);
    tol = 1e-12;
 
    opts = optimset('Display','off', ...
                    'TolFun',tol, 'TolX',tol, ...
                    'Algorithm','levenberg-marquardt', ...
                    'Jacobian','on');
 
    % Initial guess: all wages unchanged (w_hat = 1)
    x0 = ones(N,1);
 
    % Solve N unknowns with N equations: (N-1) wage eqs (skip i=1) + 1 numéraire
    [w_hat, ~, exitflag] = fsolve(@(w) wage_eq_mean_numeraire_withJ(w, lambda, tau_hat, sigma, Y, T_base_step), x0, opts);
 
    if exitflag <= 0
        warning('solve_w_hat_mean_numeraire: fsolve did not report success (exitflag=%d).', exitflag);
    end
end
 
 
function [F, J] = wage_eq_mean_numeraire_withJ(w_hat, lambda, tau_hat, sigma, Y, T)
    % Builds residuals F and analytical Jacobian J
    %   rows 1..N-1 : wage equations for i = 2..N
    %   row N       : numéraire 1 - mean(w) = 0
 
    N = length(w_hat);
 
    % ----- Precompute common objects -----
    % C_ij = (tau_ij * w_i)^(1 - sigma)
    C = (tau_hat .* (w_hat .* ones(1,N))).^(1 - sigma);   % N x N
 
    % denom_j = sum_m lambda_mj * C_mj
    denom = sum(lambda .* C, 1);                           % 1 x N
 
    % S_ij = (lambda_ij * C_ij) / denom_j
    S = (lambda .* C) ./ (ones(N,1) * denom);              % N x N
 
    % Z_j = w_j * Y_j + T_j
    Z = w_hat .* Y + T;                                    % N x 1
 
    % R_i = sum_j S_ij * Z_j  (the RHS term per i)
    R = S * Z;                                             % N x 1
 
    % H_{ik} = sum_j S_ij * Z_j * S_kj  = [ S * diag(Z) * S' ]_{ik}
    H = S * (diag(Z) * S');                                % N x N
 
    % ----- Residuals -----
    % Wage eqs for i = 2..N: F_i = w_i * Y_i - R_i
    F = zeros(N,1);
    F(1:N-1) = w_hat(2:N) .* Y(2:N) - R(2:N);
 
    % Numéraire: 1 - mean(w) = 0
    F(N) = 1 - mean(w_hat);
 
    % ----- Analytical Jacobian -----
    one_minus_sigma = 1 - sigma;
    J_core = diag(Y) ...
           - one_minus_sigma * ( diag(R ./ w_hat) - H * diag(1 ./ w_hat) ) ...
           - S * diag(Y);
 
    % Arrange rows to match F:
    % rows 1..N-1 correspond to i=2..N; last row is numéraire derivative
    J = zeros(N,N);
    J(1:N-1, :) = J_core(2:N, :);
    J(N, :)     = -ones(1,N) / N;   % d(1 - mean(w))/dw = -1/N each
 
end