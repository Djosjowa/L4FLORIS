function [saveMat,saveMatIdx] = nested_calculateDELs( WS,paramsLoop,rangeLoop,saveMat,saveMatIdx,Ld )
if nargin <= 5 % Set up initial settings
    % Setup inputData parameters
    paramsLoop = WS.inputData.parameters;
    rangeLoop  = WS.inputData.ranges;
    
    % Create N-D output matrices
    Nd            = [];
    saveMatIdx    = {};
    for i = 1:length(rangeLoop)
        Nd            = [Nd length(rangeLoop{i})];
        saveMatIdx{i} = 1;
    end
    if length(Nd) == 1
        saveMat = zeros(Nd,1);
    else
        saveMat = zeros(Nd);
    end
    Ld = 1;
end
%Nested loop for N-dimensional LUT generation
% Ld = loop depth
if length(paramsLoop) >= 1
    for i = 1:length(rangeLoop{1})
        saveMatIdx{Ld} = i;
        WS.inputData.(paramsLoop{1}) = rangeLoop{1}(i); % Update corresponding parameter
        [saveMat,saveMatIdx] = nested_calculateDELs(WS,{paramsLoop{2:end}},{rangeLoop{2:end}},saveMat,saveMatIdx,Ld+1);
        
        % Save data from DEL simulations
        if length(paramsLoop) == 1  % Save data at lowest level
            paramCombiName = nested_filenamer( WS.inputData );
            filename = [WS.foldername '\' paramCombiName];
            cd FAST_module          % Go into FAST folder
            run('windField2LUT');   % Call in-folder file for processing
            cd ..                   % Return to initial working directory
            dlmwrite([WS.DELfoldername '\' paramCombiName '.txt'],DEL);
            saveMat(saveMatIdx{:})  = DEL; % Save file names corresponding to entries
            
        end
    end
end

