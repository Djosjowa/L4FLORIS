function [ nametemp ] = nested_filenamer( inputData )
% Create name for inflow parameter set
elementNames = inputData.parameters;
nametemp = [];
for j = 1:length(elementNames)
    nametemp = [nametemp elementNames{j} '=' num2str(inputData.(elementNames{j}))];
end
end