function opto_stim(src, event, args)

hSI = src.hSI;

switch event.EventName
    
    case 'acqModeStart'
        
        disp('Starting loop')
        
        % Setup as acquisition starts
        hSI.extCustomProps.nFramesAcq = 0;
        hSI.hBeams.beamsOff();
        hSI.extCustomProps.frameCounts = {[]};
        hSI.extCustomProps.laserPower = {[]};
        hSI.hBeams.powers = 0;
       
        % Stim timing
        stimTimes = args; % [startTime, endTime]
        fps = hSI.hRoiManager.scanVolumeRate * hSI.hFastZ.numFramesPerVolume;
        hSI.extCustomProps.stimStartFrame = ceil(stimTimes(1) * fps);
        hSI.extCustomProps.stimEndFrame = ceil(stimTimes(2) * fps);        
        disp(['Stim on time: ', num2str(stimTimes(1)), ' sec'])
        disp(['Stim off time: ', num2str(stimTimes(2)), ' sec'])
        disp(['Stim on from frames ', ...
                num2str(hSI.extCustomProps.stimStartFrame), ' to ', ...
                num2str(hSI.extCustomProps.stimEndFrame)])
    case 'acqModeDone'
        
        disp('Block finished')
        
        % Trim the empty array off the end of the log data variables
        hSI.extCustomProps.frameCounts(end) = [];
        hSI.extCustomProps.laserPower(end) = [];
        
        % Save log file for the block that just finished
        optoStimInfo = hSI.extCustomProps;        
        baseName = hSI.hScan2D.logFileStem;
        saveDir = hSI.hScan2D.logFilePath;
        save(fullfile(saveDir, [baseName, '.mat']), 'optoStimInfo');
        
    case 'acqDone'
        
        % Reset variables for next trial
        hSI.extCustomProps.nFramesAcq = 0;
        hSI.extCustomProps.frameCounts{end + 1} = [];
        hSI.extCustomProps.laserPower{end + 1} = [];
        disp(['Finished with trial ' num2str(hSI.hScan2D.logFileCounter - 1)])
        
    case 'frameAcquired'
        
        % Get stim start and end frames
        stimStartFrame = hSI.extCustomProps.stimStartFrame;
        stimEndFrame = hSI.extCustomProps.stimEndFrame;

        % Increment frame counter
        fileCount = hSI.hScan2D.logFileCounter;
        hSI.extCustomProps.nFramesAcq = hSI.extCustomProps.nFramesAcq + 1;
        hSI.extCustomProps.frameCounts{fileCount}(end + 1) ...
                = hSI.extCustomProps.nFramesAcq;
        
        if hSI.extCustomProps.nFramesAcq == stimStartFrame
            % Turn laser power up at start of stim
            hSI.hBeams.powers = hSI.extCustomProps.stimPower;
            disp(['Setting laser power to ', num2str(hSI.hBeams.powers)])
        elseif hSI.extCustomProps.nFramesAcq == stimEndFrame
            % Turn laser power down at end of stim
            hSI.hBeams.powers = 0;
            disp(['Setting laser power to ', num2str(hSI.hBeams.powers)])
        end
        
        % Record laser power
        hSI.extCustomProps.laserPower{hSI.hScan2D.logFileCounter}(end + 1) ...
            = hSI.hBeams.powers;
        
end%case


end