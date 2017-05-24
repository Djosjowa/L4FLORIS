clear;
addpath bin
fillinExistingLUT = true;   % set this to true if you want to extend an existing lut, set it false if you want to make a new one
existingLUT = '..\LUT_database\LUT.mat';

%% Load parameters
inputData.DELSetName = 'dw180yaw0';   % name of the DEL set that is used
if fillinExistingLUT
    load(existingLUT);
    inputData.parameters = LUT.parameters;
    for i = 1:length(LUT.parameters)
        inputData.ranges{i} = LUT.(LUT.parameters{i});
    end
else
    load([cd '\DEL_files\' inputData.DELSetName '\parameters.mat']);   % parameters and ranges are loaded in
    inputData.parameters = parameters;
    inputData.ranges = ranges;
end

%% Build the LUT
[LUT.table,~,DELsAdded] = nested_buildLUT(inputData);
% If the LUT doesn't exist yet, generate the parameter and range fields
if ~fillinExistingLUT
    LUT.parameters = parameters;
    for i = 1:length(ranges)
        LUT.(parameters{i}) = ranges{i};
    end
end

%% Save LUT

save(['..\LUT_database\' inputData.DELSetName '.mat'],'LUT');
disp(['LUT building complete ' DELsAdded ' DEL values were added']);