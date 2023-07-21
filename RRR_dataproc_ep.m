function RRR_dataproc_ep(dirwork, dirsave, useParallel)
%
%Function to import data for ERN from the flanker, Go/NoGo, and Stroop and
% prepare those data for running through the ERP PCA Toolkit
%
%RRR_preICA(dirwork, dirsave, useParallel)
%
%Files from were recorded with EGI (.mff)
%
%
%Inputs
% dirwork - directory where raw files are located
% dirsave - directory where files should be saved
%
%Optional input
% useParallel - 0 (default), process files in serial; 1, process files in
%  parallel (requires MATLAB parallel processing toolbox)
%
%Output
% No variables are outputted to the Matlab workspace
% An eeglab dataset will be saved in the dirsave directory
%
%Required software and plugins (and tested versions)
% EEGLab (v2021.1), ERPLab (V8.20), ERP PCA Toolkit (v2.95),
%  and MFFMatlabIO (v3.7) are required for data processing and importing

%History
%by Peter Clayson (10/25/21)
%peter.clayson@gmail.com

%check whether EEGLab, ERPLab, ERP PCA Toolkit and MFFMatlabIO plugin
% are contained in the MATLAB path


fprintf('\nEnsuring dependents are found in the Matlab path\n');

%Check for EEGLab
if exist('eeglab.m','file') ~= 2
    
    dlg = {'Warning: EEGlab is not found. EEGlab may not be installed';...
        'or EEGLab may not be located in the MATLAB path.';...
        'This script requires EEGLab to run.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('EEGlab found\n');
end

%Check for ERPLab
if exist('eegplugin_erplab.m','file') ~= 2
    
    dlg = {'Warning: ERPLab is not found. ERPLab may not be installed';...
        'or ERPLab may not be located in the MATLAB path.';...
        'This script requires ERPLab to run.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('ERPlab found\n');
end

%Check for ERP PCA Toolkit plugin
if exist('ep.m','file') ~= 2
    
    dlg = {'Warning: ERP PCA Toolkit is not found.';...
        'This script requires the ERP PCA Toolkit to run.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('ERP PCA Toolkit found\n');
end


%Check for MFFMatlabIO plugin
if exist('eegplugin_mffmatlabio.m','file') ~= 2
    
    dlg = {'Warning: MFFMatlabIO plugin is not found.';...
        'This script requires the MFFMatlabIO plugin to load .mff files.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('MFFMatlabIO plugin found\n');
end

%If user does not indicate whether parallel processing should be used,
% assume that the user does NOT want to process data in parallel
if nargin < 3
    useParallel = 0;
end

%If user indicates that files should be processed in parallel, ensure that
% the parallel processing toolbox is intstalled. If the toolbox is not
% installed, run files serially.
if useParallel == 1
    if exist('parpool.m','file') ~= 2
        
        dlg = {'Warning: Parallel toolbox not installed';...
            'Data will be run serially rather than in parallel.'};
        
        for ii = 1:length(dlg)
            fprintf('%s\n',dlg{ii});
        end
        useParallel = 0;
        fprintf('\n\n');
    else
        fprintf('Parallel Toolbox installed\n');
    end
end

%Print the date and time that processing began
fprintf(...
    '\n******\nProcessing began on %s at %s \n******\n\n',...
    date, datestr(now, 'HH:MM:SS'));

%Set a variable to determine whether any files are actually found. It will
% remain 0 until appropriate files are located.
filesfound = 0;

%Pull .mff file names from dirwork.
mfflist=dir(fullfile(dirwork,'*.mff'));
if ~isempty(mfflist)
    if filesfound == 1
        filenames=[filenames,{mfflist.name}];
    elseif filesfound == 0
        filenames={mfflist.name};
        filesfound = 1;
    end
else
    dlg = {'Warning: There were no .mff files found in the specified directory.';...
        'Please specify the correct directory with .mff files to process.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
end

if filesfound == 0
    dlg = {'Warning: There were no .mff files found in the';...
        'specified directory (dirwork).';...
        'Please specify the correct directory with EEG files to process.'};
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
end

%start a timer
t1 = tic;

%create a log file where 'bad filenames' can be written for double checking
% Open file for writing
% Get current date and time
currentDateTime = datetime('now');

% Format date and time as a string
dateString = datestr(currentDateTime, 'yyyy-mm-dd_HH-MM-SS');

% Create filename using date string
fileName = ['badfiles_' dateString '.txt'];

% Open file for writing
fileID = fopen(fullfile(dirsave,fileName), 'w');

% Write some text to the file
fprintf(fileID, 'Log of bad files.\n');

% Close the file
fclose(fileID);

if useParallel == 1
    
    %use all physical cores for data processing
    nCores = feature('numCores');
    
    parpool(nCores);
    
    %Loop through all subjects in parallel
    parfor ii = 1:length(filenames)
        
        %indicate which participant is being processed and the progress
        fprintf('\n******\nProcessing participant %s\n******\n',filenames{ii});
        
        try
            
            %process the data file
            prepare_eeg(filenames{ii}, dirwork, dirsave);
            
        catch
            
            % Handle any exceptions that occurred
            fprintf('\n******\nAn error occurred for %s\n******\n',filenames{ii});
            
            % Open file for appending
            fileID = fopen(fullfile(dirsave,fileName), 'a');
            
            % Append bad file name
            fprintf(fileID, '%s\n', filenames{ii});
            
            % Close file
            fclose(fileID);
            
        end
        
    end
    
    %close the pool of workers
    delete(gcp);
    
elseif useParallel == 0
    
    %Loop through all subjects serially
    for ii = 1:length(filenames)
        
        %indicate which participant is being processed and the progress
        fprintf('\n******\nProcessing participant %s\n******\n',filenames{ii});
        
        try
            %process the data file
            prepare_eeg(filenames{ii}, dirwork, dirsave);
        catch
            % Handle any exceptions that occurred
            fprintf('\n******\nAn error occurred for %s\n******\n',filenames{ii});
            
            % Open file for appending
            fileID = fopen(fullfile(dirsave,fileName), 'a');
            
            % Append bad file name
            fprintf(fileID, '%s\n', filenames{ii});
            
            % Close file
            fclose(fileID);
        end
        
        %just keeping track of the participant being processed
        fprintf('Finished participant %d of %d\n', ii, length(filenames));
        
    end
    
end


%stop timer
t2 = toc(t1);

%incdicate how long it took to process the files
tMinutes = round(t2/60);
disp(['It took ' num2str(tMinutes) ' minutes to process these data']);

end

function prepare_eeg(subject, wrkdir, savefileshere)

%parse the name of the file, task, and order
[~,savename] = fileparts(subject);
name_split = strsplit(savename,'_');

%subject id
subjid = name_split{2};

%task for recording
task = name_split{3};

%order of presentation
torder = name_split{4};

%load .mff file
EEG = pop_mffimport(fullfile(wrkdir,subject), {'code'});

%put the name of file in the dataset
EEG.setname = savename;
EEG.subject = savename;
EEG.filename = savename;

%locate erp pca toolkit to make a path to the location of the channel
% location file
eplabdir = fileparts(which('ep.m'));

%attach channel locations
EEG = pop_chanedit(EEG, 'lookup',...
    fullfile(eplabdir,...
    'electrodes',...
    'old EGI Hydrocel',...
    'GSN-Hydrocel-129.ced'));


%Use ERPLab for filtering
% IIR Butterworth
% 4th order (24 dB/oct)
% .10 to 30 Hz half-amplitude cutoffs
% ignore online reference (129)
EEG = pop_basicfilter(EEG, 1:128,...
    'Cutoff', [.10 30],...
    'Design', 'butter',...
    'Filter', 'bandpass',...
    'Order', 4,...
    'Boundary','');

%get the nevents of events to cycle through
nevents = size(EEG.event,2);


if strcmpi(task,'flanker') %epoch flanker data
    
    task_out = 'flk';
    
    %define congruent/incongruent event keys
    cong_keys = {'>>>>>' '<<<<<'};
    inco_keys = {'>><>>' '<<><<'};
    
    
    for ii = 1:nevents
        %response-locked
        if strcmp(EEG.event(ii).code,'resp') && (ii+1 <= nevents) && ...
                (str2double(EEG.event(ii+1).mffkey_rtim) >= 100) && ...
                (str2double(EEG.event(ii+1).mffkey_rtim) < 700)
            if strcmp(EEG.event(ii+1).code,'TRSP')
                switch EEG.event(ii+1).mffkey_eval
                    case '0'
                        EEG.event(ii).type = 'err'; %error
                    case '1'
                        EEG.event(ii).type = 'cor'; %correct
                end
                EEG.event(ii).mffkeys = EEG.event(ii+1).mffkeysbackup;
            end
            
            %stimulus-locked
        elseif strcmp(EEG.event(ii).code,'Stm+') && (ii+3 <= nevents) && ...
                (str2double(EEG.event(ii+3).mffkey_rtim) >= 100) && ...
                (str2double(EEG.event(ii+3).mffkey_rtim) < 700)
            if strcmp(EEG.event(ii+3).code,'TRSP')  && ...
                    strcmp(EEG.event(ii+3).mffkey_eval,'1')
                switch EEG.event(ii+3).mffkey_Fla
                    case cong_keys
                        EEG.event(ii).type = 'con'; %congruent
                    case inco_keys
                        EEG.event(ii).type = 'inc'; %incongruent
                end
                EEG.event(ii).mffkeys = EEG.event(ii+3).mffkeysbackup;
            end
        end
    end
    
    
    %epoch response-locked
    EEG_resp = pop_epoch(EEG, {'cor' 'err'},...
        [-.4 .8], 'newname', 'epoched', 'epochinfo', 'yes');
    
    %remove unused epoch markers
    EEG_resp = pop_selectevent(EEG_resp, 'type', {'cor' 'err'},...
        'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    %baseline adjust
    EEG_resp = pop_rmbase(EEG_resp, [-400 -200]);
    
    
    %eopch stimulus-locked
    EEG_stim = pop_epoch(EEG, {'con' 'inc'},...
        [-.2 .6], 'newname', 'epoched', 'epochinfo', 'yes');
    
    %remove unused epoch markers
    EEG_stim = pop_selectevent(EEG_stim, 'type', {'con' 'inc'},...
        'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    %baseline adjust
    EEG_stim = pop_rmbase(EEG_stim, [-200 0]);
    
    
    
elseif strcmpi(task,'gng') %epoch go/nogo data
    
    task_out = 'gng';
    
    %define congruent/incongruent event keys
    go_keys = {'Target'};
    nogo_keys = {'Tilt_Left' 'Tilt_Right'};
    
    
    for ii = 1:nevents
        if strcmp(EEG.event(ii).code,'resp') && (ii+1 <= nevents) && ...
                (str2double(EEG.event(ii+1).mffkey_rtim) >= 100) && ...
                (str2double(EEG.event(ii+1).mffkey_rtim) < 700)
            if strcmp(EEG.event(ii+1).code,'TRSP')
                switch EEG.event(ii+1).mffkey_eval
                    case '0'
                        EEG.event(ii).type = 'err'; %error
                    case '1'
                        EEG.event(ii).type = 'cor'; %correct
                end
                
                EEG.event(ii).mffkeys = EEG.event(ii+1).mffkeysbackup;
                
            end
            
            
        elseif strcmp(EEG.event(ii).code,'Stm+') && (ii+3 <= nevents)
            if strcmp(EEG.event(ii+3).code,'TRSP') && ...
                    strcmp(EEG.event(ii+3).mffkey_eval,'1')
                if strcmp(EEG.event(ii+3).mffkey_Tilt,go_keys)
                    if (str2double(EEG.event(ii+3).mffkey_rtim) >= 100) && ...
                            (str2double(EEG.event(ii+3).mffkey_rtim) < 700)
                        EEG.event(ii).type = 'go'; %go trial
                        EEG.event(ii).mffkeys = EEG.event(ii+3).mffkeysbackup;
                    end
                end
                
              
            elseif strcmp(EEG.event(ii+2).code,'TRSP') && ...
                    strcmp(EEG.event(ii+2).mffkey_eval,'1') && ...
                    any(strcmp(EEG.event(ii+2).mffkey_Tilt,nogo_keys))
                
                EEG.event(ii).type = 'ng'; %go trial
                EEG.event(ii).mffkeys = EEG.event(ii+2).mffkeysbackup;
                
            end
        end
        
    end
    
    
    %epoch response-locked
    EEG_resp = pop_epoch(EEG, {'cor' 'err'},...
        [-.4 .8], 'newname', 'epoched', 'epochinfo', 'yes');
    
    %remove unused epoch markers
    EEG_resp = pop_selectevent(EEG_resp, 'type', {'cor' 'err'},...
        'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    %baseline adjust
    EEG_resp = pop_rmbase(EEG_resp, [-400 -200]);
    
    
    %epoch stimulus-locked
    EEG_stim = pop_epoch(EEG, {'go' 'ng'},...
        [-.2 .6], 'newname', 'epoched', 'epochinfo', 'yes');
    
    %remove unused epoch markers
    EEG_stim = pop_selectevent(EEG_stim, 'type', {'go' 'ng'},...
        'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    %baseline adjust
    EEG_stim = pop_rmbase(EEG_stim, [-200 0]);
    
    
elseif strcmpi(task,'stroop') %epoch stroop data
    
    task_out = 'str';
    
    for ii = 1:nevents
        if strcmp(EEG.event(ii).code,'resp') && (ii+1 <= nevents) && ...
                (str2double(EEG.event(ii+1).mffkey_rtim) >= 100) && ...
                (str2double(EEG.event(ii+1).mffkey_rtim) < 700)
            if strcmp(EEG.event(ii+1).code,'TRSP')
                switch EEG.event(ii+1).mffkey_eval
                    case '0'
                        EEG.event(ii).type = 'err'; %error
                    case '1'
                        EEG.event(ii).type = 'cor'; %correct
                end
                
                EEG.event(ii).mffkeys = EEG.event(ii+1).mffkeysbackup;
                
            end
            
        elseif strcmp(EEG.event(ii).code,'Stm+') && (ii+3 <= nevents) && ...
                (str2double(EEG.event(ii+3).mffkey_rtim) >= 100) && ...
                (str2double(EEG.event(ii+3).mffkey_rtim) < 700)
            if strcmp(EEG.event(ii+3).code,'TRSP')  && ...
                    strcmp(EEG.event(ii+3).mffkey_eval,'1')
                
                switch EEG.event(ii+3).mffkey_Wrd
                    case 'BLUE'
                        EEG.event(ii).type = 'neu'; %neutral trial
                    case 'GREEN'
                        if strcmpi(EEG.event(ii+3).mffkey_Hue,'green')
                            EEG.event(ii).type = 'con'; %congruent trial
                        else
                            EEG.event(ii).type = 'inc'; %incongruent trial
                        end
                    case 'RED'
                        if strcmpi(EEG.event(ii+3).mffkey_Hue,'red')
                            EEG.event(ii).type = 'con'; %congruent trial
                        else
                            EEG.event(ii).type = 'inc'; %incongruent trial
                        end
                end
                
                EEG.event(ii).mffkeys = EEG.event(ii+3).mffkeysbackup;
                
            end
        end
    end
    
    
    
    %epoch response-locked
    EEG_resp = pop_epoch(EEG, {'cor' 'err'},...
        [-.4 .8], 'newname', 'epoched', 'epochinfo', 'yes');
    
    %remove unused epoch markers
    EEG_resp = pop_selectevent(EEG_resp, 'type', {'cor' 'err'},...
        'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    %baseline adjust
    EEG_resp = pop_rmbase(EEG_resp, [-400 -200]);
    
    
    %eopch stimulus-locked
    EEG_stim = pop_epoch(EEG, {'con' 'inc' 'neu'},...
        [-.2 .6], 'newname', 'epoched', 'epochinfo', 'yes');
    
    %remove unused epoch markers
    EEG_stim = pop_selectevent(EEG_stim, 'type', {'con' 'inc' 'neu'},...
        'deleteevents','on','deleteepochs','on','invertepochs','off');
    
    %baseline adjust
    EEG_stim = pop_rmbase(EEG_stim, [-200 0]);
    
    
end

%Save as .set file
pop_saveset(EEG_resp, 'filename',[subjid '_' task_out '_' torder '_resp.set'],...
    'filepath', savefileshere);

%Save as .set file
pop_saveset(EEG_stim, 'filename',[subjid '_' task_out '_' torder '_stim.set'],...
    'filepath', savefileshere);


end
