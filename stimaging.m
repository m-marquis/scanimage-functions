function stimaging(src, event, varargin)
% Allows simultaneous imaging and photostimulation at a specific time during a trial and only within
% a specific ROI by modulating the overall laser power on a frame-by-frame basis. Every other frame
% during the photostim period, the laser will be switched between the stimulus and control powers
%
% First argument is a 2-element vector specifying [startTime, endTime] in seconds.
%
% Second argument specifies what laser power % to use for the photostimulation (laser power at the 
%       beginning of the acquisition will be used as the imaging laser power).

hSI = src.hSI;

switch event.EventName
    
    case 'acqModeStart'
        
        disp('Starting loop')
        
        % Setup as acquisition starts
        hSI.extCustomProps.nFramesAcq = 0;
        hSI.extCustomProps.frameCounts = {[]};
        
        stimTimes = varargin{1};
        hSI.extCustomProps.stimROIPower = varargin{2};
        hSI.extCustomProps.scanPower = hSI.extCustomProps.stimPower;
        
        fps = hSI.hRoiManager.scanVolumeRate * hSI.hFastZ.numFramesPerVolume;
        disp(['Scanning at ', num2str(hSI.extCustomProps.scanPower), '% power'])
        if ~isempty(stimTimes)
            hSI.extCustomProps.stimStartFrame = ceil(stimTimes(1) * fps);
            hSI.extCustomProps.stimEndFrame = ceil(stimTimes(2) * fps);
            if hSI.extCustomProps.stimStartFrame < 2
                hSI.extCustomProps.stimStartFrame = 2;
            end
            if hSI.extCustomProps.stimEndFrame >= hSI.hFastZ.numVolumes
                hSI.extCustomProps.stimEndFrame = hSI.hFastZ.numVolumes - 1;
            end
            disp(['Photostim at ', num2str(hSI.extCustomProps.stimROIPower), '% power'])
            disp(['Stim on time: ', num2str(stimTimes(1)), ' sec'])
            disp(['Stim off time: ', num2str(stimTimes(2)), ' sec'])
            disp(['Stim on from frames ', ...
                num2str(hSI.extCustomProps.stimStartFrame), ' to ', ...
                num2str(hSI.extCustomProps.stimEndFrame)])
        else
            hSI.extCustomProps.stimStartFrame = [];
            hSI.extCustomProps.stimEndFrame = [];
            disp('No photostimulation in this block');
        end
        
        disp('--------------------------------------------')
        disp('Starting trial 0...')
        
    case 'acqModeDone'
        
        disp('Block finished')
        
        % Trim the empty array off the end of the log data variables
        hSI.extCustomProps.frameCounts(end) = [];
        hSI.extCustomProps.stimROIPower(end) = [];
        
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
        disp(['Starting trial ' num2str(hSI.hScan2D.logFileCounter), '...'])
        
    case 'frameAcquired'
        
        % Get stim start and end frames
        stimStartFrame = hSI.extCustomProps.stimStartFrame;
        stimEndFrame = hSI.extCustomProps.stimEndFrame;
        
        % Increment frame counter
        fileCount = hSI.hScan2D.logFileCounter;
        hSI.extCustomProps.nFramesAcq = hSI.extCustomProps.nFramesAcq + 1;
        hSI.extCustomProps.frameCounts{fileCount}(end + 1) ...
            = hSI.extCustomProps.nFramesAcq;
        
        % Modify laser power as necessary
        if hSI.extCustomProps.nFramesAcq >= stimStartFrame && hSI.extCustomProps.nFramesAcq <= stimEndFrame
            
            % Set laser to stim power on even numbered frames
            if mod(hSI.extCustomProps.nFramesAcq, 2)
                hSI.hBeams.powers = hSI.extCustomProps.stimROIPower;
            else
                hSI.hBeams.powers = hSI.extCustomProps.scanPower;
            end
            
        else
            % Just turn laser off on even frames if it's not within the stim period
            if mod(hSI.extCustomProps.nFramesAcq, 2)
               hSI.hBeams.powers = 0.3;                
            else
               hSI.hBeams.powers = hSI.extCustomProps.scanPower; 
            end
        end
        
        % Record laser powers
        hSI.extCustomProps.stimROIPower{hSI.hScan2D.logFileCounter}(end + 1) ...
            = hSI.hBeams.powers;
        
end%case


end