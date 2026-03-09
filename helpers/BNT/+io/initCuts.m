% Initialize information about cut files
%
% This function detects which system (Axona/MClust) was used for cutting and processes
% information about cut files. The processing can include cut file autodetection.
% This function is used to decouple recording system from cut software.
%
function sessionData = initCuts(sessionData)
    sessionData.cutApp = bntConstants.CutApp.none;
    isMClustCuts = false;
    for e = 1:length(bntConstants.MClustExtensions)
        curExt = bntConstants.MClustExtensions{e};
        if ~isempty(dir(fullfile(sessionData.path, sprintf('*.%s', curExt))))
            isMClustCuts = true;
            break;
        end
    end
    isAxonaCuts = size(dir(sprintf('%s*.cut', sessionData.sessions{1})), 1) > 0; % NB! some versions of MClust can produce .cut files
    isAxonaCuts = isAxonaCuts | ~isempty(dir(fullfile(sessionData.path, sprintf('*%s*.cut', sessionData.basename))));
    isAxonaCuts = isAxonaCuts | ~isempty(dir(fullfile(sessionData.path, '*.cut')));
    isDbMaker = ~isempty(dir(sprintf('%s_T*C*.mat', sessionData.sessions{1})));
    
    if isDbMaker
        matFile = strcat(sessionData.sessions{1}, '_pos.mat');
        if exist(matFile, 'file')
            t = load(matFile);
            if isfield(t, 'recSystem')
                sessionData.cuts = {'*_T%uC%u.mat'};
                sessionData.cutApp = bntConstants.CutApp.mclust;
                return;
            end
        end
    end

    if ~isMClustCuts && ~isAxonaCuts
        if isempty(sessionData.cuts)
            error('BNT:badCuts', 'There is no cut information in the input file, and I''ve failed to guess any cut files. Check in your input file that data paths are correct.');
        end
        % Run both
        sessionData = io.axona.detectAxonaCuts(sessionData);
        sessionData = io.detectMClustCuts(sessionData);
        if isempty(sessionData.cutHash)
            error('BNT:badCuts', 'There is cut information in the input file, but I''ve failed to load it (find any cut files). Check in your input file that data paths are correct.');
        end
    end

    try
        % there will be an error in case we try to load .cut files from MClust
        if isAxonaCuts
            sessionData = io.axona.detectAxonaCuts(sessionData);
            sessionData.cutApp = bntConstants.CutApp.tint;
        end
    catch
    end

    if isMClustCuts
        sessionData.cutApp = bntConstants.CutApp.mclust;
        sessionData = io.detectMClustCuts(sessionData);
    end
    
    switch sessionData.cutApp
        case bntConstants.CutApp.tint
            sessionData = io.axona.calculateCutHashes(sessionData);
            
        case bntConstants.CutApp.mclust
            % do nothing for now, cut hash is initialized in detectMClustCuts
    end
end