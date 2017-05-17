


% Load FAST path
addpath(genpath('../FAST_module'))   %Load Path needed for execution of this code

% Import inflow profiles
T  = WS.inputData.T;  % simulation time in seconds
dt = WS.inputData.dt; % sampling time in seconds

% Import windfiles based on LUTsettings.m
delete('wind.*'); % Delete existing wind files and overwrite new
copyfile(['..\' filename '.wnd'],'wind.wnd');
copyfile(['..\' filename '.sum'],'wind.sum');
outputDirectory = [cd '\MLife_out\' filename];

if ~exist('outputDirectory', 'dir')
  mkdir('outputDirectory');
end

% Simulate Simulink
sim('NREL5MW_Baseline_loads_2014b')  % Simulate FAST

% Set up settings
mlifeSettings   = 'calculate_DEL_v2.mlif';  % Mlife settings file

% Simulate MLife
[fatigue, statistics] = mlife(mlifeSettings, cd, outputDirectory); 

% Process data from MLife
[~,~,rawXlsx] = xlsread([outputDirectory '\inter_result_Short-term.xlsx'],2);
DEL           = mean(cell2mat(rawXlsx(9,[4:6]))); % Take mean or? Which components DEL do we pick?
