%------------------------------------------------------------
% Script:   multisector_model.m
% Purpose:  Multi-sector, no-intermediates Armington model.
%           Computes tau and hat algebra wage changes year-by-year
%           using sector-specific elasticities, and plots the
%           emerging/advanced GDP ratio (data vs model).
%
% Inputs:   DATA_*.mat         — trade data (one file per year)
%                        
%           country_list.csv   — country labels
%           SIGMA.csv          — sector-specific elasticities (S x 1)
%
% Output:   Figure: emerging/advanced GDP ratio (data vs model)
%           Table printed to console: explained convergence (%)
%
% Requires: Optimization Toolbox (fsolve)
%------------------------------------------------------------

%% === clear all
clear; clc; close all;

%% === SETTINGS
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
sigma_s = readmatrix('SIGMA.csv');  
S = length(sigma_s);

%% === Initial-year GDP ===
X_init = Trade_country_year{1};
Y_init = sum(sum(X_init,3),2);  % N×1

%% === Storage ===
w_cum_all = cell(nFiles,1);
% Initialize baseline
w_cum_all{1} = ones(N,1);

%% === Moving Baseline Loop ===
for a = 1:(nFiles-1)
    X_base = Trade_country_year{a};

    Y_base = sum(sum(X_base,3),2);
    E_base = squeeze(sum(sum(X_base,1),3))';
    if size(E_base,1)==1, E_base = E_base'; end
    T_base = E_base - Y_base;

    E_js = squeeze(sum(X_base,1)); % N×S
    % Total expenditure shares across sectors — equal to final demand shares
    % since beta_is = 1 for all i,s in this no-intermediates model
    alpha_js = E_js ./ (E_base * ones(1,S)); 

    lambda_s = zeros(N,N,S);
    for s = 1:S
        denom = sum(X_base(:,:,s),1);
        lambda_s(:,:,s) = X_base(:,:,s) ./ (ones(N,1)*denom);
    end

    % --- ε-safe Novy τ at baseline (avoid dividing by 0) --- 
    tau_base = zeros(N,N,S);
    for s = 1:S
        Xs = X_base(:,:,s);
        pos = Xs(Xs > 0);
        eps_s = 1e-6 * min(pos);
        Xs_eps = Xs;
        Xs_eps(Xs_eps == 0) = eps_s;

        sig = sigma_s(s);
        for i = 1:N
            for j = 1:N
                tau_base(i,j,s) = ((Xs_eps(i,i)*Xs_eps(j,j))/(Xs_eps(i,j)*Xs_eps(j,i)))^(1/(2*(sig-1)));
            end
        end
    end

    % --- Next year ---
    k = a+1;
    X = Trade_country_year{k};

    tau_s = zeros(N,N,S);
    for s = 1:S
        Xs = X(:,:,s);
        pos = Xs(Xs > 0);
        eps_s = 1e-6 * min(pos);
        Xs_eps = Xs;
        Xs_eps(Xs_eps == 0) = eps_s;

        sig = sigma_s(s);
        for i = 1:N
            for j = 1:N
                tau_s(i,j,s) = ((Xs_eps(i,i)*Xs_eps(j,j))/(Xs_eps(i,j)*Xs_eps(j,i)))^(1/(2*(sig-1)));
            end
        end
    end

    tau_hat_s = tau_s ./ tau_base;

    % Solve for wages
    w_hat = solve_w_hat_multisector_mean_numeraire(lambda_s, tau_hat_s, sigma_s, alpha_js, Y_base, T_base);
    

    % Cumulate wages
    w_cum_all{k} = w_cum_all{k-1} .* w_hat;

    fprintf("Processed: baseline %d → %d\n", year_labels(a), year_labels(k));
end

%% === GDP: Data vs Model ===
GDP_data = zeros(N,nFiles);
GDP_model = zeros(N,nFiles);

for k = 1:nFiles
    GDP_data(:,k) = sum(sum(Trade_country_year{k},3),2);
end

for k = 1:nFiles
    GDP_model(:,k) = w_cum_all{k} .* Y_init;
end

%% === emerging vs advanced ratio ===
is_emerging = ismember(country_list, emerging_codes);
is_advanced = ismember(country_list, advanced_codes);

ratio_data = sum(GDP_data(is_emerging,:),1) ./ sum(GDP_data(is_advanced,:),1);
ratio_model = sum(GDP_model(is_emerging,:),1) ./ sum(GDP_model(is_advanced,:),1);

figure; hold on;
plot(year_labels, ratio_data, 'o-', 'LineWidth',1.4, 'DisplayName',"Data");
plot(year_labels, ratio_model, 's--', 'LineWidth',1.4, 'DisplayName',"Model");
xlabel('Year'); ylabel('emerging/advanced GDP Ratio');
title('emerging vs advanced (Data vs GE Model)');
legend('Location','NorthWest'); grid on;

%% === Explained convergence table
dData  = ratio_data  - ratio_data(1);
dModel = ratio_model - ratio_model(1);
explained_pct = 100 * (dModel ./ dData);
explained_pct(1) = NaN;
explained_pct(dData == 0) = NaN;
yearCols = matlab.lang.makeValidName(string(year_labels'));
explainedTable = array2table(explained_pct, ...
    'VariableNames', cellstr(yearCols), ...
    'RowNames',      {'Explained (%)'});
disp('=== Explained convergence: multisector model ===');
disp(explainedTable);

%% ==== Multisector solver: mean-wage numéraire + analytical Jacobian ====
function w_hat = solve_w_hat_multisector_mean_numeraire(lambda_s, tau_hat_s, sigma_s, alpha_js, Y, T)
    % Numéraire: mean(w_hat) = 1
    % Unknowns: all N wages 
    % Equations: N-1 wage eqs (i=2..N) + 1 numéraire
    N   = length(Y);
    tol = 1e-12;
 
    opts = optimset('Display','off', ...
                    'TolFun',tol, 'TolX',tol, ...
                    'Algorithm','levenberg-marquardt', ...
                    'Jacobian','on');
 
    x0 = ones(N,1);
    [w_hat, ~, exitflag] = fsolve(@(w) wage_eq_multisector_mean_numeraire_withJ( ...
        w, lambda_s, tau_hat_s, sigma_s, alpha_js, Y, T), x0, opts);
 
    if exitflag <= 0
        warning('solve_w_hat_multisector_mean_numeraire: fsolve did not report success (exitflag=%d).', exitflag);
    end
end
 
function [F, J] = wage_eq_multisector_mean_numeraire_withJ(w_hat, lambda_s, tau_hat_s, sigma_s, alpha_js, Y, T)
    % Multisector wage system with mean-wage numéraire.
    % Equation ordering:
    %   rows 1..N-1 : wage equations for i = 2..N
    %   row N       : numéraire 1 - mean(w) = 0
 
    N = length(w_hat);
    S = length(sigma_s);
 
    % Common objects
    Z = w_hat .* Y + T;        % N x 1
 
    % Accumulators over sectors
    R  = zeros(N,1);           % RHS per i (sum_s R_s)
    Jc = zeros(N,N);           % core Jacobian contribution to be subtracted from diag(Y)
 
    for s = 1:S
        sig  = sigma_s(s);
        one_minus_sig = 1 - sig;
 
        % C_s(i,j) = (tau_hat_s(i,j,s) * w_i)^(1 - sigma_s)
        C_s = (tau_hat_s(:,:,s) .* (w_hat * ones(1,N))).^(one_minus_sig);   % N x N
 
        % denom_s(j) = sum_m lambda_s(m,j,s) * C_s(m,j)
        denom_s = sum(lambda_s(:,:,s) .* C_s, 1);                           % 1 x N
 
        % S_s(i,j) = [lambda_s(i,j,s) * C_s(i,j)] / denom_s(j)
        S_s = (lambda_s(:,:,s) .* C_s) ./ (ones(N,1) * denom_s);            % N x N
 
        % Sector-weighted expenditure: Zs = α_js(:,s) .* Z
        Zs = alpha_js(:,s) .* Z;                                            % N x 1
 
        % R_s = S_s * Zs
        R_s = S_s * Zs;                                                     % N x 1
        R   = R + R_s;
 
        % H_s = S_s * diag(Zs) * S_s'
        H_s = S_s * (diag(Zs) * S_s');                                      % N x N
 
        % Contribution to core Jacobian:
        J_core_s = one_minus_sig * ( diag(R_s ./ w_hat) - H_s * diag(1 ./ w_hat) ) ...
                 + S_s * diag(alpha_js(:,s) .* Y);
 
        Jc = Jc + J_core_s;
    end
 
    % Residuals:
    % Wage eqs i=2..N: F_i = w_i*Y_i - R_i
    F = zeros(N,1);
    F(1:N-1) = w_hat(2:N) .* Y(2:N) - R(2:N);
 
    % Numéraire: mean(w) = 1
    F(N) = 1 - mean(w_hat);
 
    % Jacobian:
    % J_core = diag(Y) - Jc
    J_core = diag(Y) - Jc;
 
    % Arrange rows to match F’s ordering
    J = zeros(N,N);
    J(1:N-1, :) = J_core(2:N, :);
    J(N, :)     = -ones(1,N)/N;   % d(1 - mean(w))/dw = -1/N each
end
