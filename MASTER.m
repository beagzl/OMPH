%% MASTER Script for HANK/RANK Model Simulations
% This is the main control script for running HANK and RANK model simulations.
%
% Workflow:
%   1. Generate model files using Buildfile.m (run once per model)
%   2. Configure parameters in model_parameters.m
%   3. Run simulations for HANK and/or RANK models
%
% Author: Beatriz Gonzalez, Galo Nuño, Dominik Thaler
% Date: January 2026


% modelfilename: 'RANK' for RANK model; 'HANK' for HANK model; 


%% ========================================================================
% SECTION 1: GENERATE MODEL FILES
% ========================================================================
% Run this section once to create Dynare .mod files for each model type.
% Only needs to be rerun if you modify model equations in Buildfile.m

clear, clc;
modelfilename = 'RANK';      
buildfile;


clear, clc;
modelfilename = 'HANK';      
buildfile;


%% ========================================================================
% SECTION 2: RUN SIMULATIONS
% ========================================================================
% Configure simulation parameters in model_parameters.m before running.
%
% Example configuration for comparing RANK vs HANK Ramsey policy in 
% response to a time-preference shock:
% -------------------------------------------------
% OMP = 1;              % Ramsey optimal monetary policy
% timeless = 1;         % Timeless perspective
% shocktype = 2;        % Time preference shock
% shocksize = -100;     % -100 basis points (-1%)

clear; clc;

modelfilename = 'HANK';      
run_mod;

modelfilename = 'RANK';      
run_mod;



