%------------------------------------------------------------------------
%
% run_mod.m
%
% January 2026
%
% This file:
%   1. Loads model parameters and options from model_parameters.m
%   2. Calls Dynare with appropriate macro definitions
%   3. Runs the model simulation
%
% Usage:
%   Simply run this script after setting options in model_parameters.m
%   buildfile.m needs to be run first to create needed mod files
%
% Replication codes: " Firm Heterogeneity, Capital Misallocation and 
% Optimal Monetary Policy"
%
% (c) 2026 Beatriz Gonzalez, Galo Nuño and Dominik Thaler
%
%
%------------------------------------------------------------------------



%% Load Parameters and Options

model_parameters

%% Run dynare

% Call Dynare with arguments (set in model_parameters_noK)
eval(['dynare ', modelfilename, ' nostrict savemacro -DN=' num2str(N),...
    ' -DOMP=' num2str(OMP), ...
    ' -DIc=' num2str(Ic),' -DIf=' num2str(If),' -DIs=' num2str(Is),...
    ' -Dmodelfilename=["',modelfilename,'"]']); 
