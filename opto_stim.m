function opto_stim(src, event, args)

hSI = src.hSI;

switch event.EventName
    
    case 'acqModeStart'
        
        disp('Starting loop')
        
        % Get handles for stimulus and control ROIs
        hStimROI = hSI.hRoiManager.roiGroupMroi.rois(1);
        hControlROI = hSI.hRoiManager.roiGroupMroi.rois(2);
        
        % Setup as acquisition starts
        hSI.extCustomProps.nFramesAcq = 0;
        hSI.hBeams.beamsOff();
        hSI.extCustomProps.frameCounts = {[]};
        hSI.extCustomProps.stimROIPower = {[]};
        hSI.extCustomProps.controlROIPower = {[]};
        hSI.hRoiManager.roiGroupMroi.rois(1).powers = 0.3;
        hSI.hRoiManager.roiGroupMroi.rois(2).powers = hSI.extCustomProps.stimPower;
       
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
        disp('--------------------------------------------')
        disp('Starting trial 0...')
    case 'acqModeDone'
        
        disp('Block finished')
        
        % Trim the empty array off the end of the log data variables
        hSI.extCustomProps.frameCounts(end) = [];
        hSI.extCustomProps.stimROIPower(end) = [];
        hSI.extCustomProps.controlROIPower(end) = [];
        
        % Save log file for the block that just finished
        optoStimInfo = hSI.extCustomProps;        
        baseName = hSI.hScan2D.logFileStem;
        saveDir = hSI.hScan2D.logFilePath;
        save(fullfile(saveDir, [baseName, '.mat']), 'optoStimInfo');
        
    case 'acqDone'
        
        % Reset variables for next trial
        hSI.extCustomProps.nFramesAcq = 0;
        hSI.extCustomProps.frameCounts{end + 1} = [];
        hSI.extCustomProps.stimROIPower{end + 1} = [];
        hSI.extCustomProps.controlROIPower{end + 1} = [];
        disp(['Starting trial ' num2str(hSI.hScan2D.logFileCounter), '...'])
        
    case 'frameAcquired'
        
        % Get handles for stimulus and control ROIs
        hStimROI = hSI.hRoiManager.roiGroupMroi.rois(1);
        hControlROI = hSI.hRoiManager.roiGroupMroi.rois(2);
        
        % Get stim start and end frames
        stimStartFrame = hSI.extCustomProps.stimStartFrame;
        stimEndFrame = hSI.extCustomProps.stimEndFrame;

        % Increment frame counter
        fileCount = hSI.hScan2D.logFileCounter;
        hSI.extCustomProps.nFramesAcq = hSI.extCustomProps.nFramesAcq + 1;
        hSI.extCustomProps.frameCounts{fileCount}(end + 1) ...
                = hSI.extCustomProps.nFramesAcq;
        
        if hSI.extCustomProps.nFramesAcq == stimStartFrame
            % Switch laser power to stim ROI
            hSI.hRoiManager.roiGroupMroi.rois(1).powers = hSI.extCustomProps.stimPower;
            hSI.hRoiManager.roiGroupMroi.rois(2).powers = 0.3;
            disp('Setting laser to full power in stim ROI')
%             disp(['Setting stim ROI power to ', num2str(hStimROI.powers)])
%             disp(['Setting control ROI power to ', num2str(hControlROI.powers)])
            
            hControlROI.enable = 0;

        elseif hSI.extCustomProps.nFramesAcq == stimEndFrame
            % Switch laser power to control ROI
            hSI.hRoiManager.roiGroupMroi.rois(1).powers = 0.3;
            hSI.hRoiManager.roiGroupMroi.rois(2).powers = hSI.extCustomProps.stimPower;
            disp('Setting laser to full power in control ROI')
%             disp(['Setting stim ROI power to ', num2str(hStimROI.powers)])
%             disp(['Setting control ROI power to ', num2str(hControlROI.powers)])
            hControlROI.enable = 1;
        end
        
        % Record laser powers
        hSI.extCustomProps.stimROIPower{hSI.hScan2D.logFileCounter}(end + 1) ...
            = hStimROI.powers;
        hSI.extCustomProps.controlROIPower{hSI.hScan2D.logFileCounter}(end + 1) ...
            = hControlROI.powers;
        
end%case


end