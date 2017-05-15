function [saveMat,saveMatIdx] = nestedloop_generateInflows( paramsLoop,rangeLoop,inputData,saveMat,saveMatIdx,Ld )
if nargin <= 5 % Set up initial settings
    % Create N-D output matrices
    Nd            = [];
    saveMat       = {};
    saveMatIdx    = {};
    for j = 1:length(rangeLoop)
        Nd            = [Nd length(rangeLoop{j})];
        saveMatIdx{j} = 1;
    end;
    Ld = 1;
end;
yWake = inputData.(paramsLoop{1});
yWake_range = rangeLoop{1}(i);
%Nested loop for N-dimensional LUT generation
% Ld = loop depth
if length(paramsLoop) >= 1
    parfor i = 1:length(rangeLoop{1})
%         saveMatIdx{Ld} = i;
        yWake = yWake_range(i); % Update corresponding parameter
%         [saveMat,saveMatIdx] = nested_generateInflows({paramsLoop{2:end}},{rangeLoop{2:end}},inputData,saveMat,saveMatIdx,Ld+1);
        
        % Save data
        if length(paramsLoop) == 1  % Save data at lowest level
            filename = ['inflowProfiles\' inputData.destinationFolder '\' nested_filenamer( inputData )];
            fun_generateInflow(inputData,filename);
            saveMat{1}  = filename; % Save file names corresponding to entries
        end;
    end;
end

