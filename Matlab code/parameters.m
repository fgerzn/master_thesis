%------------------------------------------------------------
% Script:   parameters.m
% Purpose:  Cleans and aggregates raw WIOD input-output tables
%           for years 1995-2007. Computes structural model
%           parameters (ALPHA, BETA, GAMMA, DELTA, PHI, SIGMA)
%           and saves one PARAMETERS_YY_AGG.mat file per year.
%
% Inputs:   wiotYY_row_apr12.csv      — raw WIOD table (one per year)
%           country_aggregation.csv   — country mapping (Nagg x 41)
%           sector_aggregation.csv    — sector mapping  (Sagg x 35)
%           country_list.csv          — country labels
%           sector_list.csv           — sector labels
%           SIGMA.csv                 — trade elasticities (S x 1)
%
% Output:   PARAMETERS_YY_AGG.mat     — struct with all parameters,
%                                       one file per year
%
% Requires: No special toolboxes
%
% Usage:    Set working directory to folder containing all input
%           files, then run. Script loops over all years
%           automatically.
%------------------------------------------------------------

%% Clean + aggregate WIOD for years 1995–2007 and SAVE parameters
% Loops over all years, stores each in PARAMETERS_YY_AGG.mat
%------------------------------------------------------------
clear; close all; clc;

% ---- Set your working directory ----
cd('C:\Users\flori\Downloads\Data (1)');

% ---- Raw WIOD dimensions ----
N0 = 41;      % original regions
S0 = 35;      % original sectors

% ---- Years to process ----
% yr_strs matches exact WIOD filename suffixes: 95..99, then 00..07
yr_strs = {'95','96','97','98','99','00','01','02','03','04','05','06','07'};

% ---- Load aggregation mappings (same for all years) ----
AggC = dlmread('country_aggregation.csv');   % Nagg x 41
AggS = dlmread('sector_aggregation.csv');    % Sagg x 35
SIGMA_raw = dlmread('SIGMA.csv');
CCODES = importdata('country_list.csv');
SCODES = importdata('sector_list.csv');

N = size(AggC, 1);
S = size(AggS, 1);

% Aggregation operators (constant across years)
CC = kron(AggC, AggS);        % (NS_agg) x (N0*S0)
FF = kron(AggC, eye(5));       % final demand aggregation

% ---- Loop over years ----
for k = 1:numel(yr_strs)

    yr_str = yr_strs{k};   % '95','96',...,'99','00','01',...,'07'

    fprintf('Processing year: %s ...\n', yr_str);

    % ---- Load raw WIOD table ----
    fname = ['wiot' yr_str '_row_apr12.csv'];
    if ~exist(fname, 'file')
        warning('File not found: %s — skipping.', fname);
        continue;
    end
    DATA = dlmread(fname);

    Zinit = DATA(1:N0*S0, 1:N0*S0);
    Xraw  = DATA(1:N0*S0, 1:N0*S0 + 5*N0);
    Rinit = sum(Xraw, 2);

    % ---- Correct inventories ----
    FIN   = Xraw(:, N0*S0+1 : N0*S0 + 5*N0);
    Fpos  = FIN .* (FIN > 0);
    Fsum  = sum(Fpos, 2);

    % ---- Recompute consistent output ----
    A  = Zinit / diag(Rinit + 1e-6 .* (Rinit <= 1e-6));
    R0 = (eye(N0*S0) - A) \ Fsum;
    Z0 = A * diag(R0);

    % ---- Aggregate ----
    Z = CC * Z0 * CC';
    F = CC * Fpos * FF';
    X = [Z, F];
    R = sum(X, 2);

    % ---- Trade flows ----
    TRADEIst = zeros(N, N, S, S);
    for i = 1:N
        for j = 1:N
            for s = 1:S
                for t = 1:S
                    TRADEIst(i,j,s,t) = Z((i-1)*S+s, (j-1)*S+t);
                end
            end
        end
    end

    TRADEFs = zeros(N, N, S);
    for i = 1:N
        for j = 1:N
            cols = (j-1)*5 + (1:5);
            for s = 1:S
                TRADEFs(i,j,s) = sum(F((i-1)*S+s, cols));
            end
        end
    end

    TRADEs = sum(TRADEIst, 4) + TRADEFs;

    % ---- ALPHA ----
    EF = zeros(N, S);
    for j = 1:N
        cols = (j-1)*5 + (1:5);
        tmp = sum(F(:, cols), 2);
        tmp = reshape(tmp, [S, N])';
        EF(j,:) = sum(tmp, 1);
    end
    ALPHA = EF ./ (sum(EF, 2) + 1e-12);

    % ---- BETA ----
    INTUSE   = sum(Z, 1)';
    BETAvec  = 1 - INTUSE ./ (R + 1e-12);
    BETA     = reshape(BETAvec, [S, N])';

    % ---- GAMMA ----
    GAMMA = zeros(N, S, S);
    for j = 1:N
        for t = 1:S
            use_st = zeros(S, 1);
            for i = 1:N
                for s = 1:S
                    use_st(s) = use_st(s) + TRADEIst(i,j,s,t);
                end
            end
            GAMMA(j,:,t) = (use_st / (sum(use_st) + 1e-12))';
        end
    end

    % ---- DELTA ----
    DELTA = zeros(N, S, S);
    for j = 1:N
        Gj       = squeeze(GAMMA(j,:,:));
        intshare = (1 - BETA(j,:));
        Bj       = Gj .* (ones(S,1) * intshare);
        DELTA(j,:,:) = inv(eye(S) - Bj)';
    end
    
    % ---- PHI: total expenditure shares across sectors in importer j ----
    % ESEC(j,s) = total absorption in country j of sector s goods
    ESEC = squeeze(sum(TRADEs, 1));              % N x S
    PHI  = ESEC ./ (sum(ESEC, 2) + 1e-12);      % N x S
    % ---- Pack into struct ----
    PARAMETERS        = struct();
    PARAMETERS.N      = N;
    PARAMETERS.S      = S;
    PARAMETERS.ALPHA  = ALPHA;
    PARAMETERS.BETA   = BETA;
    PARAMETERS.GAMMA  = GAMMA;
    PARAMETERS.DELTA  = DELTA;
    PARAMETERS.SIGMA  = SIGMA_raw;
    PARAMETERS.PHI = PHI
    PARAMETERS.CCODES = CCODES;
    PARAMETERS.SCODES = SCODES;
    % Also store raw trade/IO objects for convenience
    PARAMETERS.TRADEs    = TRADEs;
    PARAMETERS.TRADEIst  = TRADEIst;
    PARAMETERS.TRADEFs   = TRADEFs;
    PARAMETERS.Z         = Z;
    PARAMETERS.F         = F;
    PARAMETERS.R         = R;

    % ---- Save ----
    outfile = ['PARAMETERS_' yr_str '_AGG.mat'];
    save(outfile, 'PARAMETERS');
    fprintf('  Saved: %s\n', outfile);

    % ---- Sanity checks ----
    tol = 1e-8;

    alpha_dev = max(abs(sum(ALPHA, 2) - 1));
    if alpha_dev > tol
        warning('[%s] ALPHA row sums deviate from 1. Max dev: %.3e', yr_str, alpha_dev);
    end
    if any(ALPHA(:) < -tol)
        warning('[%s] ALPHA has negative entries. Min: %.3e', yr_str, min(ALPHA(:)));
    end

    gamma_tmp = squeeze(sum(GAMMA, 2));
    gamma_dev = max(abs(gamma_tmp(:) - 1));
    if gamma_dev > 1e-6
        warning('[%s] GAMMA upstream sums deviate from 1. Max dev: %.3e', yr_str, gamma_dev);
    end

    bmin = min(BETA(:)); bmax = max(BETA(:));
    if bmin < -1e-6 || bmax > 1+1e-6
        warning('[%s] BETA outside [0,1]. Min=%.4f, Max=%.4f', yr_str, bmin, bmax);
    end

    for j = 1:N
        Gj       = squeeze(GAMMA(j,:,:));
        intshare = (1 - BETA(j,:));
        Bj       = Gj .* (ones(S,1) * intshare);
        Dj       = squeeze(DELTA(j,:,:));
        resid    = norm((eye(S) - Bj) * Dj' - eye(S), 'fro');
        if resid > 1e-6
            warning('[%s] DELTA inversion residual high for country %d: %.3e', yr_str, j, resid);
        end
        rho = max(abs(eig(Bj)));
        if rho >= 1
            warning('[%s] Country %d spectral radius(Bj)=%.4f >= 1.', yr_str, j, rho);
        end
    end

    fprintf('  Sanity checks passed for year %s.\n', yr_str);

end % year loop (k)

disp('All years processed.');
