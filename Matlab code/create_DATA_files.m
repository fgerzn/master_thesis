%------------------------------------------------------------
% Script:   create_DATA_files.m
% Purpose:  Cleans and aggregates raw WIOD input-output tables
%           for years 1995-2007 and saves TRADEs for each year.
%           Follows Costinot and Rodriguez-Clare (2014).
%
% Inputs:   wiotYY_row_apr12.csv    — raw WIOD table (one per year)
%           country_aggregation.csv — country mapping (Nagg x 41)
%           sector_aggregation.csv  — sector mapping  (Sagg x 35)
%           country_list.csv        — country labels
%           sector_list.csv         — sector labels
%
% Output:   DATA_YY.mat             — contains TRADEs, CCODES, SCODES
%                                     one file per year
%
% Requires: No special toolboxes
%
% Usage:    Set working directory to folder containing all input
%           files, then run. Script loops over all years automatically.
%------------------------------------------------------------
clear; close all; clc;

% ---- Set your working directory ----
cd('C:\Users\flori\Downloads\Data (1)');

% ---- Raw WIOD dimensions ----
N0 = 41;   % original number of regions
S0 = 35;   % original number of sectors

% ---- Years to process ----
yr_strs = {'95','96','97','98','99','00','01','02','03','04','05','06','07'};

% ---- Load aggregation mappings (same for all years) ----
AggC = dlmread('country_aggregation.csv');   % Nagg x 41
AggS = dlmread('sector_aggregation.csv');    % Sagg x 35
CCODES = importdata('country_list.csv');
SCODES = importdata('sector_list.csv');

N = size(AggC, 1);
S = size(AggS, 1);

% Aggregation operators (constant across years)
CC = kron(AggC, AggS);     % aggregation operator for intermediate flows
FF = kron(AggC, eye(5));   % aggregation operator for final demand

% ---- Loop over years ----
for k = 1:numel(yr_strs)

    yr_str = yr_strs{k};
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

    % ---- Correct inventories: set negative final demand entries to 0 ----
    FIN  = Xraw(:, N0*S0+1 : N0*S0 + 5*N0);
    Fpos = FIN .* (FIN > 0);
    Fsum = sum(Fpos, 2);

    % ---- Recompute consistent output under fixed coefficients ----
    A  = Zinit / diag(Rinit + 1e-6 .* (Rinit <= 1e-6));
    R0 = (eye(N0*S0) - A) \ Fsum;
    Z0 = A * diag(R0);

    % ---- Aggregate countries and sectors ----
    Z = CC * Z0 * CC';
    F = CC * Fpos * FF';
    X = [Z, F];
    R = sum(X, 2);

    % ---- Intermediate trade flows: i->j, upstream s -> downstream t ----
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

    % ---- Final trade flows: i->j, sector s ----
    TRADEFs = zeros(N, N, S);
    for i = 1:N
        for j = 1:N
            cols = (j-1)*5 + (1:5);
            for s = 1:S
                TRADEFs(i,j,s) = sum(F((i-1)*S+s, cols));
            end
        end
    end

    % ---- Total trade flows: intermediate + final ----
    TRADEs = sum(TRADEIst, 4) + TRADEFs;

    % ---- Save ----
    outfile = ['DATA_' yr_str '.mat'];
    save(outfile, 'TRADEs', 'CCODES', 'SCODES');
    fprintf('  Saved: %s\n', outfile);

end

disp('All years processed.');