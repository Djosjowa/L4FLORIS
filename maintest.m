clear all; close all; clc;
addpath('bin');                      % add binary files
addpath('bin\FLORISSE_M');           % add FLORIS model files
addpath('bin\FLORISSE_M\functions'); % add FLORIS model functions

% Plot results
plotResults = false;

% Import DEL look-up table
DEL_table = load('./LUT_database/LUT_Ben.mat'); % Load LUT of choice

% Load model, turbine and topology settings
modelStruct = floris_param_model('default');    % Load default FLORIS model settings
modelStruct.axialIndProvided = true;
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
optimStruct.optConst        = 0.5;                          % Weighting factor. Power only = 1, Loads only = 0.
optimStruct.iterations      = 100;                          % Optimization iterations  [-]
optimStruct.maxYaw          = +30;                          % Largest  yaw angle [radians]
optimStruct.minYaw          = -30;                          % Smallest yaw angle [radians]
optimStruct.windUncertainty = [-12:4:12];                   % Additional wind disturbance range (symmetric)
optimStruct.minA            = 0;
optimStruct.maxA            = 0.31;
optimStruct.Pref            = 10*10^6;     % Reference power [W]
optimStruct.Pbandwidth      = 0.05*10^6;

%optimization
input0 = [.25.*ones(1,9),zeros(1,9)];
yawmin = -30;
yawmax = 30;
amin   = 0;
amax   = 0.31;
lb     = [amin*ones(1,9),yawmin*ones(1,9)];
ub     = [amax.*ones(1,9),yawmax*ones(1,9)];
A = [];
b = [];

output = @(input)optimfunc(input,optimStruct,siteStruct,modelStruct,turbType,DEL_table,plotResults);

[input,fval] = fmincon(output,input0,A,b,lb,ub)






