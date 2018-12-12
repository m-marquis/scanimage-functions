function [] = fly_tracker_server_MM_chrimson()
% fly_tracker_server 
% example for using the ScanImage API to set up a grab
hSI = evalin('base','hSI');             % get hSI from the base workspace

%%%%
% Start a server to listen on the socket to fly tracker
%%%%
clear t;
PORT = 30000;
t = tcpip('0.0.0.0', PORT, 'NetworkRole', 'server');
set(t, 'InputBufferSize', 30000); 
set(t, 'TransferDelay', 'off');
disp(['About to listen on port: ' num2str(PORT)]);
fopen(t);
pause(1.0);
disp('Client connected');

while 1
    
    % Read the message     
    while (t.BytesAvailable == 0)
        pause(0.1);
    end
        
    data = fscanf(t, '%s');
    data = strtrim(data);
    disp(data);

    if( strcmp(data,'END_OF_SESSION') == 1 )
        break;
    end
    
    % Extract trial info
    durCell = regexp(data, '(?<=dur_).*(?=_nTrials)', 'match');
    blockDur = str2num(durCell{:});
    nTrialCell = regexp(data, '(?<=nTrials_).*', 'match');
    nTrials = str2num(nTrialCell{:});
    hSI.hScan2D.logFileStem = data; % set the base file name for the Tiff file   
    
    % Set up stack acquisition
    hSI.hFastZ.enable = true;
    
    % Set volumes to scan
    volumes_per_second = hSI.hRoiManager.scanVolumeRate;
    disp(['VPS: ' num2str( volumes_per_second )]);
    hSI.hFastZ.numVolumes = int32(ceil((blockDur/nTrials) * volumes_per_second));
    disp(['Number of volumes per trial: ' num2str( hSI.hFastZ.numVolumes )]);
    
    % Configure acquitions
    hSI.hChannels.loggingEnable = true;    % Make sure logging is enabled  
    hSI.hScan2D.logFileCounter = 1;         % Set file counter to 1    
    hSI.acqsPerLoop = nTrials;             
    hSI.extTrigEnable = true;               % Enable external trigger 
    
    % Make sure user functions are enabled for laser power modulation
    hSI.extCustomProps.stimPower = hSI.hBeams.powers;
    for iFun = 1:numel(hSI.hUserFunctions.userFunctionsCfg)
        hSI.hUserFunctions.userFunctionsCfg(iFun).Enable = true;
    end
    
    % Start acquisitions
    pause(1.0);
    hSI.startLoop();
    
    % Signal back to the fly tracker client that it can start daq and image
    % acquisition.
    fprintf(t, 'SI51_Acq_1');
end

% Clean up
hSI.extTrigEnable = false;
hSI.hBeams.powers = hSI.extCustomProps.stimPower;

%close the socket
fclose(t);
delete(t);

end