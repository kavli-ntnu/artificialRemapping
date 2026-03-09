function sessionData = detectRecordingSystem(sessionData)    
    % Find out if this is Axona or NeuraLynx data
    % In case of multiple sessions, assume that combined sessions have been recorded using the same system.
    % Thus, check only the first session.
    setCandidate = sprintf('%s.set', sessionData.sessions{1});
    if ~isempty(dir(setCandidate))
        % check if this is a Virmen virtual reality, which comes alongside of Axona
        virmenCandidate = sprintf('%s*.vr', sessionData.sessions{1});
        isVirmen = ~isempty(dir(virmenCandidate));
        if isVirmen
            setVirmen();
        else
            setAxona();
        end
        return;
    end

    if ~isempty(dir(fullfile(sessionData.sessions{1}, '*.nvt')))
        setNeuraLynx;
    end
    
    %% check if it is in old dbMaker format
    matFile = sprintf('%s_pos.mat', sessionData.sessions{1});
    if exist(matFile, 'file') ~= 0
        tt = load(matFile);
        if isfield(tt, 'recSystem')
            if strcmpi(tt.recSystem, 'Axona')
                setAxona();
            else
                setNeuraLynx();
            end
        end
    end
    
    function setAxona()
        sessionData.system = bntConstants.RecSystem.Axona;
        sessionData.sampleTime = 0.02; % 50 Hz
        sessionData.videoSamplingRate = 50;
    end
        
    function setVirmen()
        sessionData.system = bntConstants.RecSystem.Axona;
        sessionData.sampleTime = 0.02; % 50 Hz
        sessionData.videoSamplingRate = 50;
    end

    function setNeuraLynx()
        sessionData.system = bntConstants.RecSystem.Neuralynx;
        sessionData.sampleTime = 0.04; % 25 Hz
        sessionData.videoSamplingRate = 25;
    end            
end
