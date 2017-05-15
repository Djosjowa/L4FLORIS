clear all; close all; clc; tic;
addpath bin % Add binary files (external functions)

% LUT settings: overwrite default parameters
inputData.parameters       = {'yWake'}; %{'yWake','u_mean'};    % Variable names to cover
inputData.ranges           = {-20:10:20}; %{-20:20:20,6:1:8};     % Ranges to cover for LUT.parameters

% Set up settings
inputData.T                 = 200;          % Simulation duration [seconds]
inputData.dt                = 0.250;        % Timestep for FAST simulations
inputData.plotProfile       = false;        % Plot the inflow profiles (y,z)
inputData.HH                = 90.0;         % Hub height
inputData.zWake             = inputData.HH; % Height of wake center
inputData.yWake             = 0;            % Lateral C2C displacement between rotor and wake
inputData.yaw               = 0;           % yaw angle of turbine [degrees]
inputData.destinationFolder = 'C2C_only';   % Destination folder name

% Static settings: vertical grid
inputData.y               = 300:-5:-300;   % lateral dimension (NOTE: MUST BE POSITIVE TO NEGATIVE)
inputData.z               = 0:5:200;       % vertical dimension

% FAST settings and LUT coverage
inputData.u_fs            = 8.0;          % Freestream wind speed [m/s]
inputData.Gaussian_A      = 5;            % Wake deficit in center [m/s]
inputData.Gaussian_omegay = 50;           % Spread for Gaussian wake shape
inputData.Gaussian_omegaz = 50;           % Spread for Gaussian wake shape
% https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

% Setup destination folder
foldername = ['inflowProfiles/' inputData.destinationFolder];
% if exist(foldername) % Delete existing folder
%     user_response = input(['Are you sure you want to delete folder ' foldername '? (y/n)  '],'s');
%     if lower(user_response(1)) == 'y'
%         rmdir(foldername ,'s');
%     else
%         error('Please rename the old folder or change the destination folder name.');
%     end;
% end;
if exist(foldername)
    rmdir(foldername ,'s');
end


mkdir(foldername); % Create folder

% Fill up N-D output filenames matrix & generate inflow profiles
[filesMat,~] = nested_generateInflows(inputData.parameters,inputData.ranges,inputData);

% Save workspace for future use
save(['inflowProfiles/' inputData.destinationFolder '/workspace.mat']);
disp('Finished generating inflow profiles.');
toc;