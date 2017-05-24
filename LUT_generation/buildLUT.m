clear;
addpath bin
fillinExistingLUT = true;   % set this to true if you want to extend an existing lut, set it false if you want to make a new one

%% Load parameters
existingLUT = '..\LUT_database\LUT.mat';    % path of existing LUT that will be extended if fillinExistingLUT is true
inputData.DELSetName = 'Dwake255yaw-30';   % name of the DEL set that is used
if fillinExistingLUT
    load(existingLUT);
    inputData.parameters = LUT.parameters;
    for i = 1:length(LUT.parameters)
        inputData.ranges{i} = LUT.(LUT.parameters{i});
    end
    inputData.table = LUT.table;
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
disp(['LUT building complete ' num2str(DELsAdded) ' DEL values were added']);

% Check if user want to save this new LUT
user_response = input(['Are you sure you want to save this LUT (y/n)  '],'s');
if lower(user_response(1)) == 'y'
    if fillinExistingLUT
        save(existingLUT,'LUT');
    else
        save(['..\LUT_database\' inputData.DELSetName '.mat'],'LUT');
    end
else
    disp('LUT not saved');
end
