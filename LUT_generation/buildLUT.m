clear;
addpath bin

%% Load parameters
inputData.DELSetName = '4_parameters';   % name of the DEL set that is used
load([cd '\DEL_files\' inputData.DELSetName '\parameters.mat']);   % parameters and ranges are loaded in
inputData.parameters = parameters;
inputData.ranges = ranges;

%% Build the LUT
[LUT.table,~] = nested_buildLUT(inputData);
LUT.parameters = parameters;
for i = 1:length(ranges)
    LUT.(parameters{i}) = ranges{i};
end
save(['..\LUT_database\' inputData.DELSetName '.mat'],'LUT');