clear all; close all; clc; tic;
addpath bin % Add binary files (external functions)

% LUT settings: overwrite default parameters
inputData.parameters       = {'yWake','Dwake','yaw','U_fs'};   % Variable names to cover
inputData.ranges           = {-20:5:20,200:100:300,-20:2:20,[6 8]};  % Ranges to cover for LUT.parameters

% Set up settings
inputData.T                 = 200;          % Simulation duration [seconds]
inputData.dt                = 0.250;        % Timestep for FAST simulations
inputData.plotProfile       = false;        % Plot the inflow profiles (y,z)
inputData.HH                = 90.0;         % Hub height
inputData.zWake             = inputData.HH; % Height of wake center
inputData.yWake             = 0;            % Lateral C2C displacement between rotor and wake
inputData.yaw               = 0;           % yaw angle of turbine [degrees]
inputData.destinationFolder = 'C2C_only';   % Destination folder name
inputData.doWindShear       = false;         % Take wind shear into account
inputData.shear_const       = 0.1;          % Constant for shear effect on smooth surfaces(oceans, sand, etc.)

% Static settings: vertical grid
inputData.y               = 300:-5:-300;   % lateral dimension (NOTE: MUST BE POSITIVE TO NEGATIVE)
inputData.z               = 0:5:200;       % vertical dimension

% FAST settings and LUT coverage
inputData.U_fs            = 8.0;          % Freestream wind speed [m/s]
inputData.Gaussian_A      = 5;            % Wake deficit in center [m/s]
inputData.Gaussian_omegay = 50;           % Spread for Gaussian wake shape
inputData.Gaussian_omegaz = 50;           % Spread for Gaussian wake shape
inputData.Drotor          = 126.4;        % Diameter rotor
inputData.Dwake           = 176.4;        % Diameter outerwake
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
maxLength = 0;
% find the longest parameter vector
for i = 1:length(inputData.ranges)
    if length(inputData.ranges{i}) > maxLength
        longestRangeIdx = i;
        maxLength = inputData.ranges{i}; 
    end
end
longestRange = inputData.ranges{longestRangeIdx};
% make newParameters and newRanges with the longest parameter vector put in
% first
if longestRangeIdx == 1
    newParameters = inputData.parameters;
    newRanges = inputData.ranges(2:end);
elseif longestRangeIdx == length(inputData.ranges)
    newParameters = [inputData.parameters(end) , inputData.parameters(1:end-1)];
    newRanges = inputData.ranges(1:end-1);
else
    newParameters = [inputData.parameters(longestRangeIdx) , inputData.parameters(1:longestRangeIdx-1) , inputData.parameters(longestRangeIdx+1:end)];
    newRanges = [inputData.ranges(1:longestRangeIdx-1) , inputData.ranges(longestRangeIdx+1:end)];
end
parfor i = 1:length(inputData.ranges{1})
    nested_generateInflows(newParameters,[longestRange(i) , newRanges],inputData);
end

% Save workspace for future use
save(['inflowProfiles/' inputData.destinationFolder '/workspace.mat']);
disp('Finished generating inflow profiles.');
toc;