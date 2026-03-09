function sessionData = detectTetrodes(sessionData)
    if strcmpi(sessionData.system, bntConstants.RecSystem.Neuralynx) == false
        error('Invalid recording system');
    end
    
    possibleTetrodes = 1:16;
    tetrodes = [];
    if isdir(sessionData.sessions{1})
        workDir = sessionData.sessions{1};
    else
        workDir = sessionData.path;
    end
    
    for i = 1:length(possibleTetrodes)
        candidate = possibleTetrodes(i);
        tFile = sprintf('%s/TT%u.ntt', workDir, candidate);
        if exist(tFile, 'file') == 0
            continue;
        end
        tetrodes = cat(1, tetrodes, candidate);
    end
    if ~isempty(tetrodes)
        units = [tetrodes nan(length(tetrodes), 1)];
        sessionData.units = units;
    end    
end