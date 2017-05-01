clear all; close all; clc;
addpath('bin');                      % add binary files
addpath('bin\FLORISSE_M');           % add FLORIS model files
addpath('bin\FLORISSE_M\functions'); % add FLORIS model functions

% Show figures
plotResults = true;

% Import DEL look-up table
DEL_table = load('./LUT_database/lut_example.mat'); % Load LUT of choice

% Load model, turbine and topology settings
modelStruct = floris_param_model('default');    % Load default FLORIS model settings
turbType    = floris_param_turbine('nrel5mw');  % Load NREL 5MW turbine properties

siteStruct.LocIF =   [300,    100.0,  turbType.hub_height % 9 turbine scenario
                      300,    300.0,  turbType.hub_height
                      300,    500.0,  turbType.hub_height
                      1000,   100.0,  turbType.hub_height
                      1000,   300.0,  turbType.hub_height
                      1000,   500.0,  turbType.hub_height
                      1600,   100.0,  turbType.hub_height
                      1600,   300.0,  turbType.hub_height
                      1600,   500.0,  turbType.hub_height];

% Atmospheric settings
siteStruct.uInfIf   = 8;        % x-direction flow speed inertial frame (m/s)
siteStruct.vInfIf   = 0;        % y-direction flow speed inertial frame (m/s)
siteStruct.rho      = 1.1716;   % Atmospheric air density (kg/m3)

% Setup optimization settings
optimStruct.optConst        = 1;                            % Weighting factor. Power only = 1, Loads only = 0.
optimStruct.iterations      = 1000;                          % Optimization iterations  [-]
optimStruct.maxYaw          = +30;                          % Largest  yaw angle [radians]
optimStruct.minYaw          = -30;                          % Smallest yaw angle [radians]
optimStruct.axInd           = 1/3*ones(size(siteStruct.LocIF,1)); % Axial induction factors
optimStruct.windUncertainty = [-12:4:12];                   % Additional wind disturbance range (symmetric)

Pref                        = 10000*10^6;                   % Reference Power [W]

% Run optimization
[yaw_opt,J_Pws_opt,J_DEL_opt,J_sum_opt] = optimizeL4FLORIS(modelStruct,turbType,siteStruct,optimStruct,plotResults,Pref);