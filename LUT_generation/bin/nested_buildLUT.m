function [saveMat,saveMatIdx,counter] = nested_buildLUT( inputData,saveMat,saveMatIdx,paramsLoop,rangeLoop,counter,Ld )
if nargin <= 5 % Set up initial settings
    % Setup inputData parameters
    paramsLoop = inputData.parameters;
    rangeLoop  = inputData.ranges;
    % Create N-D output matrices
    Nd            = [];
    saveMatIdx    = {};
    for i = 1:length(rangeLoop)
        Nd            = [Nd length(rangeLoop{i})];
        saveMatIdx{i} = 1;
    end
    if isfield(inputData,'table')
        saveMat = inputData.table;
    elseif length(Nd) == 1
        saveMat = zeros(Nd,1);
    else
        saveMat = zeros(Nd);
    end
    Ld = 1;
    counter.DELsAdded = 0;
    counter.DELsOverwritten = 0;
end
%Nested loop for N-dimensional LUT table building
% Ld = loop depth
if length(paramsLoop) >= 1
    for i = 1:length(rangeLoop{1})
        saveMatIdx{Ld} = i;
        inputData.(paramsLoop{1}) = rangeLoop{1}(i); % Update corresponding parameter
        [saveMat,saveMatIdx,counter] = nested_buildLUT(inputData,saveMat,saveMatIdx,{paramsLoop{2:end}},{rangeLoop{2:end}},counter,Ld+1);
        
        % Read data from txt files
        if length(paramsLoop) == 1  % Save data at lowest level
            DELfile = [cd '\DEL_files\' inputData.DELSetName '\' nested_filenamer(inputData) '.txt'];
            if exist(DELfile,'file')
                DEL = dlmread(DELfile);
                if saveMat(saveMatIdx{:}) == 0
                    saveMat(saveMatIdx{:})  = DEL; % Save file names corresponding to entries
                    counter.DELsAdded = counter.DELsAdded+1;
                elseif inputData.doOverwrite
                    saveMat(saveMatIdx{:})  = DEL; % Save file names corresponding to entries
                    counter.DELsOverwritten = counter.DELsOverwritten+1;
                end
                
            end
        end
    end
end



