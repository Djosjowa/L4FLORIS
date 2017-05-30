clear; close all; clc;
addpath('bin');                      % add binary files
addpath('bin\FLORISSE_M');           % add FLORIS model files
addpath('bin\FLORISSE_M\functions'); % add FLORIS model functions

% Plot results
plotResults = true;

% Import DEL look-up table
load('./LUT_database/4D_LUT_complete_cleaned.mat'); % Load LUT of choice

% Load model, turbine and topology settings
modelStruct = floris_param_model('default');    % Load default FLORIS model settings
modelStruct.axialIndProvided = true;            % Value set to true so FLORIS uses user defined axial induction factors

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
optimStruct.iterations      = 100;                            % Optimization iterations  [-]
optimStruct.minYaw          = -30;                          % Smallest yaw angle [radians]
optimStruct.maxYaw          = +30;                          % Largest  yaw angle [radians]
optimStruct.initYaw         = 0;                            % Yaw angle used for first iteration [radians]
optimStruct.minA            = 0;                            % Smallest axial induction factor
optimStruct.maxA            = 1/3;                         % Largest axial induction factor
optimStruct.initA           = 0.08;                          % Axial induction factor used for first iteration
optimStruct.windUncertainty = [-12:4:12];                   % Additional wind disturbance range (symmetric)
optimStruct.Pref            = 8*10^6;     % Reference power [W]
optimStruct.Pbandwidth      = 0.05*10^6;

% Run optimization
[a_opt,a_tries,yaw_opt,yaw_tries,J_Pws_opt,J_DEL_opt,J_sum_opt] = optimizeL4FLORIS(modelStruct,turbType,siteStruct,optimStruct,LUT,plotResults);

for i = 1:optimStruct.iterations-1
    a_diff(i,:) = a_tries(i+1,:)-a_tries(i,:);
    yaw_diff(i,:) = yaw_tries(i+1,:)-yaw_tries(i,:);
end
avgstepA = mean(abs(a_diff(a_diff~=0)))    
avgstepYaw = mean(abs(yaw_diff(yaw_diff~=0))) 