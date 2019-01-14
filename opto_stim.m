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
%         hSI.extCustomProps.laserPower = {[]};
        hSI.extCustomProps.stimROIPower = {[]};
        hSI.extCustomProps.controlROIPower = {[]};
        hStimROI.powers = 0.3;
        hControlROI.powers = hSI.extCustomProps.stimPower;
        hSI.hBeams.updateBeamBufferAsync(true);
        %         hSI.hBeams.powers = 0.3;
        
        % Get stim timing
        stimTimes = args; % [startTime, endTime]
        fps = hSI.hRoiManager.scanVolumeRate * hSI.hFastZ.numFramesPerVolume;
        if ~isempty(stimTimes)
            hSI.extCustomProps.stimStartFrame = ceil(stimTimes(1) * fps);
            hSI.extCustomProps.stimEndFrame = ceil(stimTimes(2) * fps);
            disp(['Stim on time: ', num2str(stimTimes(1)), ' sec'])
            disp(['Stim off time: ', num2str(stimTimes(2)), ' sec'])
            disp(['Stim on from frames ', ...
                num2str(hSI.extCustomProps.stimStartFrame), ' to ', ...
                num2str(hSI.extCustomProps.stimEndFrame)])
        else
            hSI.extCustomProps.stimStartFrame = [];
            hSI.extCustomProps.stimEndFrame = [];
<<<<<<< HEAD
            disp('No photostimulation in this trial');
        end
        
=======
            disp('Not using photostim in this trial...')
        end
>>>>>>> 8dd7c0075e93887185104d55cca11123a3374754
        disp('--------------------------------------------')
        disp('Starting trial 0...')
        
    case 'acqModeDone'
        
        disp('Block finished')
        
        % Trim the empty array off the end of the log data variables
        hSI.extCustomProps.frameCounts(end) = [];
%         hSI.extCustomProps.laserPower(end) = [];
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
%         hSI.extCustomProps.laserPower{end + 1} = [];
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
<<<<<<< HEAD
        if ~isempty(stimStartFrame)
=======

        if ~isempty(stimStartFrame) && ~isempty(stimEndFrame)
>>>>>>> 8dd7c0075e93887185104d55cca11123a3374754
            if hSI.extCustomProps.nFramesAcq == stimStartFrame
                % Switch laser power to stim ROI
                hStimROI.powers = hSI.extCustomProps.stimPower;
                hControlROI.powers = 0.3;
<<<<<<< HEAD
                disp(['Setting laser to ', num2str(hSI.extCustomProps.stimPower), '% power in stim ROI'])
=======
                hSI.hBeams.updateBeamBufferAsync(true);
                disp(['Setting laser to ', num2str(hSI.extCustomProps.stimPower), '% power in stim ROI'])
                %               hSI.hBeams.powers = hSI.extCustomProps.stimPower;
                
>>>>>>> 8dd7c0075e93887185104d55cca11123a3374754
                
            elseif hSI.extCustomProps.nFramesAcq == stimEndFrame
                % Switch laser power to control ROI
                hStimROI.powers = 0.3;
                hControlROI.powers = hSI.extCustomProps.stimPower;
<<<<<<< HEAD
                disp(['Setting laser to ', num2str(hSI.extCustomProps.stimPower), '% power in control ROI'])
                
=======
                hSI.hBeams.updateBeamBufferAsync(true);
                disp(['Setting laser to ', num2str(hSI.extCustomProps.stimPower), '% power in control ROI'])
                %               hSI.hBeams.powers = 0.3;
>>>>>>> 8dd7c0075e93887185104d55cca11123a3374754
            end
        end
        
        % Record laser powers
        hSI.extCustomProps.stimROIPower{hSI.hScan2D.logFileCounter}(end + 1) ...
            = hStimROI.powers;
        hSI.extCustomProps.controlROIPower{hSI.hScan2D.logFileCounter}(end + 1) ...
            = hControlROI.powers;
%           hSI.extCustomProps.laserPower{hSI.hScan2D.logFileCounter}(end + 1)...
%               = hSI.hBeams.powers;
        
end%case


end