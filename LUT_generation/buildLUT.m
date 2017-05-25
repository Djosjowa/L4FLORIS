clear; clc;
addpath bin
fillinExistingLUT = true;   % set this to true if you want to extend an existing LUT, set it false if you want to make a new one
doOverwrite = false;        % do you want to overwrite existing DEL values in the LUT?

%% Load parameters
existingLUT = '..\LUT_database\4D_LUT_3yaws.mat';    % path of existing LUT that will be extended if fillinExistingLUT is true
inputData.DELSetName = '4_parameters';   % name of the DEL set that is used
inputData.doOverwrite = doOverwrite;
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
[LUT.table,~,counter] = nested_buildLUT(inputData);
% If the LUT doesn't exist yet, generate the parameter and range fields
if ~fillinExistingLUT
    LUT.parameters = parameters;
    for i = 1:length(ranges)
        LUT.(parameters{i}) = ranges{i};
    end
end

%% Save LUT
disp(['LUT building complete ' num2str(counter.DELsAdded) ' DEL values were added, and ' num2str(counter.DELsOverwritten) ' were overwritten']);

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
