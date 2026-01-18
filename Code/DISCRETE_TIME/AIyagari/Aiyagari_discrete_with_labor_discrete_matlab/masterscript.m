%==========================================================================
% Code to solve a heterogeneous agent neoclassical model a la 
% Aiyagari (1994,QJE) with idiosyncratic productivity and labor supply.
%
% Michael Dobrew, michael.dobrew@gmail.com, December 2020
%==========================================================================

% _________________________________________________________________________
%% Initialize workspace and load directories
% _________________________________________________________________________
clc
clear
close all

addpath(genpath('Aiyagari_with_labor_discrete\functions'))
casename = 'SS_Aiyagari_labor';
dir = pwd;

% _________________________________________________________________________
%% Initialize parameters and grid
% _________________________________________________________________________
disp('Initializing parameters and grid')
setParameters
makegrids

% _________________________________________________________________________
%% Solve for steady state capital, policy functions and distribution
% _________________________________________________________________________
disp('Solving Steady State by EGM')
[pol, JD, SS] = main_steadystate(par, mpar, grid, mesh, Pz);

% _________________________________________________________________________
%% Calculate steady state capital and further statistics
% _________________________________________________________________________

SS_stats

%% Save
filename=casename;
save(filename)
