function [saveMat,saveMatIdx] = nestedloop_generateInflows( paramsLoop,rangeLoop,inputData,saveMat,saveMatIdx,Ld )
if nargin <= 5 % Set up initial settings
    % Create N-D output matrices
    Nd            = [];
    saveMat       = {};
    saveMatIdx    = {};
    for i = 1:length(rangeLoop)
        Nd            = [Nd length(rangeLoop{i})];
        saveMatIdx{i} = 1;
    end
    Ld = 1;
end
%Nested loop for N-dimensional LUT generation
% Ld = loop depth
if length(paramsLoop) >= 1
    for i = 1:length(rangeLoop{1})
        saveMatIdx{Ld} = i;
        inputData.(paramsLoop{1}) = rangeLoop{1}(i); % Update corresponding parameter
        [saveMat,saveMatIdx] = nested_generateInflows({paramsLoop{2:end}},{rangeLoop{2:end}},inputData,saveMat,saveMatIdx,Ld+1);
        
        % Save data
        if length(paramsLoop) == 1  % Save data at lowest level
            filename = ['inflowProfiles\' inputData.destinationFolder '\' nested_filenamer( inputData )];
            fun_generateInflow(inputData,filename);
            saveMat{saveMatIdx{:}}  = filename; % Save file names corresponding to entries
        end
%     [saveMat,saveMatIdx] = nested_generateInflows_loop( i,paramsLoop,rangeLoop,Ld,saveMat,saveMatIdx,inputData);
    end
end

