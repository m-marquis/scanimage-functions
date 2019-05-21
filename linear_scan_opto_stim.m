function linear_scan_opto_stim(src, event, varargin)
% Custom user function to allow online control of laser power in multiple linear scan ROIs at specific times in
% each trial by counting the acquired frames and updating ROI power accordingly.
%
% First argument is a 2-element vector specifying [startTime, endTime] in seconds.
%
% Second argument specifies what laser power % to use for the photostimulation (laser power at the 
%       beginning of the acquisition will be used as the imaging laser power).

hSI = src.hSI;

switch event.EventName
    
    case 'acqModeStart'
        
        disp('Starting grab')
        
        % Get handles for stimulus and control ROIs
        allRois = hSI.hRoiManager.roiGroupLineScan;
        scanRois = []; scanRoiNums = [];
        for iRoi = 1:numel(allRois.rois)
            currRoiName = allRois.rois(iRoi).scanfields.shortDescription;
            if ~strcmp(currRoiName(7:end), 'pause') && ...
                        ~strcmp(currRoiName(7:end), 'park')
                scanRoiNums(end + 1) = iRoi;
            end
        end
        scanRois = allRois.rois(scanRoiNums);
        hStimRoi = scanRois(1).scanfields;
        hControlRoi = scanRois(2).scanfields;
        hImageRois = scanRois(3:end);
        
        % Setup laser power as acquisition starts
        hSI.extCustomProps.nFramesAcq = 0;
        hSI.hBeams.beamsOff();
        hSI.extCustomProps.frameCounts = [];
        hSI.extCustomProps.stimROIPowerLog = [];
        hSI.extCustomProps.controlROIPowerLog = [];
        hSI.extCustomProps.interleaveTrials = 0;
        hStimRoi.powers = 0.1;
        hControlRoi.powers = hSI.extCustomProps.imagingPower;
        for iRoi = 1:numel(hImageRois)
            hImageRois(iRoi).scanfields.powers = hSI.extCustomProps.imagingPower;
        end
        hSI.hBeams.updateBeamBufferAsync(true);
        
        % Calculate trial start timing
        stimTimes = varargin{1}; % [startTime, endTime]
        hSI.extCustomProps.stimROIPower = varargin{2}; % Get stim timing power
        cps = hSI.hRoiManager.scanFrameRate;
        cpt = hSI.extCustomProps.cyclesPerTrial;
        trialStartCycles = [];
        for iTrial = 1:hSI.extCustomProps.nTrials
            if iTrial == 1
                trialStartCycles(iTrial) = 1;
            else
                trialStartCycles(iTrial) = 1 + (cpt * iTrial);
            end
        end
        hSI.extCustomProps.trialStartCycles = trialStartCycles;
        
        % Calculate stim timing
        if ~isempty(stimTimes)
            stimStartCycles = [];
            stimEndCycles = [];
            relStimStartCycle = ceil(stimTimes(1) * cps) - 1;
            relStimEndCycle = ceil(stimTimes(2) * cps) - 1;
            if relStimStartCycle < 1
                relStimStartCycle = 1;
            end
            if relStimEndCycle >= cpt - 1
                relStimEndCycle = cpt - 2;
            end
            for iTrial = 1:hSI.extCustomProps.nTrials
                stimStartCycles(iTrial) = trialStartCycles(iTrial) + relStimStartCycle;
                stimEndCycles(iTrial) = trialStartCycles(iTrial) + relStimEndCycle;
            end
            hSI.extCustomProps.stimStartCycles = stimStartCycles;
            hSI.extCustomProps.stimEndCycles = stimEndCycles;
            disp(['Photostim at ', num2str(hSI.extCustomProps.stimROIPower), '% power'])
            disp(['Stim on time: ', num2str(stimTimes(1)), ' sec'])
            disp(['Stim off time: ', num2str(stimTimes(2)), ' sec'])
            disp(['Stim on from cycles ', ...
                num2str(relStimStartCycle), ' to ', ...
                num2str(relStimEndCycle), ' of each trial'])
        else
            hSI.extCustomProps.stimStartCycles = [];
            hSI.extCustomProps.stimEndCycles = [];
            disp('No photostimulation in this block');
        end
        disp('--------------------------------------------')
        disp('Starting trial 1...')
        
    case 'acqModeDone'
        
        disp('Block finished')
%         
%         % Trim the empty array off the end of the log data variables
%         hSI.extCustomProps.frameCounts(end) = [];
%         hSI.extCustomProps.stimROIPowerLog(end) = [];
%         hSI.extCustomProps.controlROIPowerLog(end) = [];
        
        % Save log file for the block that just finished
        optoStimInfo = hSI.extCustomProps;        
        baseName = hSI.hScan2D.logFileStem;
        saveDir = hSI.hScan2D.logFilePath;
        save(fullfile(saveDir, [baseName, '.mat']), 'optoStimInfo');
        
        % Clear extCustomProps
        hSI.hBeams.powers = hSI.extCustomProps.imagingPower;
        hSI.extCustomProps = [];
        
    case 'frameAcquired'
        
        % Get handles for stimulus and control ROIs
        allRois = hSI.hRoiManager.roiGroupLineScan.rois;
        scanRois = []; scanRoiNums = [];
        for iRoi = 1:numel(allRois)
            currRoiName = allRois(iRoi).scanfields.shortDescription;
            if ~strcmp(currRoiName(7:end), 'pause') && ...
                        ~strcmp(currRoiName(7:end), 'park')
                scanRoiNums(end + 1) = iRoi;
            end
        end
        scanRois = allRois(scanRoiNums);
        hStimRoi = scanRois(1).scanfields;
        hControlRoi = scanRois(2).scanfields;
        
        % Get trial and stim timing
        trialStartCycles = hSI.extCustomProps.trialStartCycles;
        stimStartCycles = hSI.extCustomProps.stimStartCycles;
        stimEndCycles = hSI.extCustomProps.stimEndCycles;
        
        % Increment frame counter
        hSI.extCustomProps.nFramesAcq = hSI.extCustomProps.nFramesAcq + 1;
        hSI.extCustomProps.frameCounts(end + 1) = hSI.extCustomProps.nFramesAcq;
        disp(['frameCounterForDisplay count: ', num2str(hSI.frameCounterForDisplay)])
        disp(['frameAcquired event count: ', num2str(hSI.extCustomProps.frameCounts(end))])

        
        % Modify laser power if necessary
        if ~isempty(stimStartCycles)                
                if find(stimStartCycles == hSI.extCustomProps.nFramesAcq)
                    % Switch laser power to stim ROI
                    hStimRoi.powers = hSI.extCustomProps.stimROIPower;
                    hControlRoi.powers = 0.3;
                    hSI.hBeams.updateBeamBufferAsync(true);
                    disp(['Setting laser to ', num2str(hSI.extCustomProps.stimROIPower), '% power in stim ROI'])
                    
                elseif find(stimEndCycles == hSI.extCustomProps.nFramesAcq)
                    % Switch laser power to control ROI
                    hStimRoi.powers = 0.3;
                    hControlRoi.powers = hSI.extCustomProps.stimROIPower;
                    hSI.hBeams.updateBeamBufferAsync(true);
                    disp(['Setting laser to ', num2str(hSI.extCustomProps.stimROIPower), '% power in control ROI'])
                    
                end
        end
        
        % Notify user if starting new trial
        if find(trialStartCycles == hSI.extCustomProps.nFramesAcq)
           disp(['Starting trial #', num2str(find(trialStartCycles == hSI.extCustomProps.nFramesAcq))]); 
        end
        
        % Record laser powers
        hSI.extCustomProps.stimROIPowerLog(end + 1) ...
            = hStimRoi.powers;
        hSI.extCustomProps.controlROIPowerLog(end + 1) ...
            = hControlRoi.powers;
        
end%case


end