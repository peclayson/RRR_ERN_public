function RRR_ERN_scoring(wrkdir,savefile)

if ~iscell(wrkdir)
    rawfilesloc = dir(fullfile(wrkdir,'*_rb.ept'));
    files = struct2cell(rawfilesloc)';
    files = files(:,1:2);
elseif iscell(wrkdir)
    for ii = 1:size(wrkdir,2)
        if ii == 1
            rawfilesloc = dir(fullfile(wrkdir{ii},'*_rb.ept'));
            files = struct2cell(rawfilesloc)';
            files = files(:,1:2);
        else
            rawfilesloc = dir(fullfile(wrkdir{ii},'*_rb.ept'));
            newfiles = struct2cell(rawfilesloc)';
            newfiles = newfiles(:,1:2);
            
            files = [files; newfiles];
        end
    end
end

ERPinfo.resp.ern_ROI = {'E6'};
ERPinfo.resp.ern_meanwind = [0 100];
ERPinfo.resp.pe_ROI = {'E62'};
ERPinfo.resp.pe_meanwind = [200 400];
ERPinfo.resp.events = {'cor' 'err'};

ERPinfo.stim.flk.n1_ROI = {'E6'};
ERPinfo.stim.flk.n1_meanwind = [90 130];
ERPinfo.stim.flk.n2_ROI = {'E6'};
ERPinfo.stim.flk.n2_meanwind = [240 300];

ERPinfo.stim.str.n1_ROI = {'E6'};
ERPinfo.stim.str.n1_meanwind = [90 120];
ERPinfo.stim.str.n2_ROI = {'E6'};
ERPinfo.stim.str.n2_meanwind = [240 300];


ERPinfo.stim.gng.n1_ROI = {'E6'};
ERPinfo.stim.gng.n1_meanwind = [100 130];
ERPinfo.stim.gng.n2_ROI = {'E6'};
ERPinfo.stim.gng.n2_meanwind = [240 300];


ERPinfo.stim.flk.events = {'con' 'inc'};
ERPinfo.stim.str.events = {'con' 'inc' 'neu'};
ERPinfo.stim.gng.events = {'go' 'ng'};

global EPmain EPtictoc

EPtictoc.start=[];
EPtictoc.step=1;
EPtictoc.stop=0;

ep_tictoc('begin');

% Make N by 2 matrix of fieldname + value type
variable_names_types = [["subjid", "cell"]; ...
    ["task", "cell"]; ...
    ["torder", "cell"]; ...
    ["event", "cell"]; ...
    ["ern", "double"]; ...
    ["pe", "double"]];
%     ["urtrl", "double"]...
 


% Make table using fieldnames & value types from above
resp_master = table(...
    'Size',[(size(files,1)*900),size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));

variable_names_types = [["subjid", "cell"]; ...
    ["task", "cell"]; ...
    ["torder", "cell"]; ...
    ["event", "cell"]; ...
    ["n1", "double"];...
    ["n2", "double"]];
%     ["urtrl", "double"]...
    


% Make table using fieldnames & value types from above
stim_master = table(...
    'Size',[(size(files,1)*900),size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));

resp_tbl_ind = 1;
stim_tbl_ind = 1;
resp_individ = [];
stim_individ = [];

%Print the date and time that processing began
fprintf(...
    '\n******\nProcessing began on %s at %s \n******\n\n',...
    date, datestr(now, 'HH:MM:SS'));


%start a timer
t1 = tic;

for ii = 1:size(files,1)
    %indicate which participant is being processed and the progress
    fprintf('\n******\nProcessing participant %s\n******\n',files{ii});
    
    
    fileloc = fullfile(char(files(ii,2)),char(files(ii,1)));
    
    %load ep dataset
    EPdata = ep_readData('file',fileloc,'format','ep_mat');
    
    if contains(EPdata.dataName,'resp')
        resp_individ = resp_scoretrial(EPdata, ERPinfo, EPdata.dataName);
        
        if ~isempty(resp_individ)
            tbl_nrow = height(resp_individ);
            
            resp_master(resp_tbl_ind:(resp_tbl_ind+tbl_nrow-1),:) = resp_individ;
            resp_tbl_ind = resp_tbl_ind + tbl_nrow;
            
        end
        
    elseif contains(EPdata.dataName,'stim')
        stim_individ = stim_scoretrial(EPdata, ERPinfo, EPdata.dataName);
        
        if ~isempty(stim_individ)
            tbl_nrow = height(stim_individ);
            
            stim_master(stim_tbl_ind:(stim_tbl_ind+tbl_nrow-1),:) = stim_individ;
            stim_tbl_ind = stim_tbl_ind + tbl_nrow;
            
        end 
    end
    
end



resp_master = resp_master(1:(resp_tbl_ind-1),:);
writetable(resp_master,fullfile(savefile,...
    strcat('RRR_ERN_',datestr(now, 'mmddyy'),'_resp.csv')));

stim_master = stim_master(1:(stim_tbl_ind-1),:);
writetable(stim_master,fullfile(savefile,...
    strcat('RRR_ERN_',datestr(now, 'mmddyy'),'_stim.csv')));


%stop timer
t2 = toc(t1);

%incdicate how long it took to process the files
tMinutes = round(t2/60);
disp(['It took ' num2str(tMinutes) ' minutes to process these data.']);


end

function individ = resp_scoretrial(EPdata, ERPinfo, subjid)

str = strsplit(subjid,'_');
subjid = str{1};
task = str{2};
torder = str{3};

start_ind = 1;

for ii = start_ind:length(EPdata.cellNames)
    
    event = EPdata.cellNames(ii);
    
    if strcmp(event,'trial')
        event = EPdata.events{1,ii}.value;
    end
    
    if (EPdata.analysis.badTrials(ii) == 0)
        
        ern_mean = zeros(0,length(ERPinfo.resp.ern_ROI));
        
        for jj = 1:length(ERPinfo.resp.ern_ROI)
            
            chnind = find(strcmp(EPdata.chanNames,...
                ERPinfo.resp.ern_ROI{jj})); 
            
            ern_mean(jj) = mean(EPdata.data(chnind,...
                knnsearch(EPdata.timeNames,ERPinfo.resp.ern_meanwind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.resp.ern_meanwind(2)),...
                ii)); %#ok<*FNDSB>
            
        end
        
        if length(ern_mean) > 1
            ern_mean = mean(ern_mean);
        end
        
        
        pe_mean = zeros(0,length(ERPinfo.resp.pe_ROI));
        
        for jj = 1:length(ERPinfo.resp.pe_ROI)
            
            chnind = find(strcmp(EPdata.chanNames,...
                ERPinfo.resp.pe_ROI{jj})); 
            
            pe_mean(jj) = mean(EPdata.data(chnind,...
                knnsearch(EPdata.timeNames,ERPinfo.resp.pe_meanwind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.resp.pe_meanwind(2)),...
                ii)); %#ok<*FNDSB>
            
        end
        
        if length(pe_mean) > 1
            pe_mean = mean(pe_mean);
        end
        
        
        if all(size(EPdata.events{1,ii}) == 1)
            sub = EPdata.events{1,ii};
        else
            sub = EPdata.events{1,ii}(find(strcmp({EPdata.events{1,ii}.type}, 'resp')== 1));
        end
        
        keys = sub.keys;
        mffkeys = keys(find(strcmp({keys.code}, 'mffkeys')==1)).data;
        mffkeys_parsed = regexp(mffkeys,'\d*','Match');
                
        if ~exist('individ','var')
            
            individ = table;
            individ.subjid = cellstr(subjid);
            individ.task = cellstr(task);
            individ.torder = cellstr(torder);
            individ.event = cellstr(event);
            individ.ern = ern_mean;
            individ.pe = pe_mean;
%             individ.urtrl = str2double(mffkeys_parsed{6});
            
        elseif exist('individ','var')
            
            row = table;
            row.subjid = cellstr(subjid);
            row.task = cellstr(task);
            row.torder = cellstr(torder);
            row.event = cellstr(event);
            row.ern = ern_mean;
            row.pe = pe_mean;
%             row.urtrl = str2double(mffkeys_parsed{6});
            
            individ = vertcat(individ,row); %#ok<AGROW>
            
        end
        
    end
    
end

if ~exist('individ','var')
    individ = [];
end

end


function individ = stim_scoretrial(EPdata, ERPinfo, subjid)


str = strsplit(subjid,'_');
subjid = str{1};
task = str{2};
torder = str{3};

start_ind = 1;

for ii = start_ind:length(EPdata.cellNames)
    
    event = EPdata.cellNames(ii);
    
    if strcmp(event,'trial')
        event = EPdata.events{1,ii}.value;
    end
    
    if (EPdata.analysis.badTrials(ii) == 0) 
        
        n1_mean = zeros(0,length(ERPinfo.stim.(task).n1_ROI));
        
        for jj = 1:length(ERPinfo.stim.(task).n1_ROI)
            
            chnind = find(strcmp(EPdata.chanNames,...
                ERPinfo.stim.(task).n1_ROI{jj})); 
            
            n1_mean(jj) = mean(EPdata.data(chnind,...
                knnsearch(EPdata.timeNames,ERPinfo.stim.(task).n1_meanwind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.stim.(task).n1_meanwind(2)),...
                ii)); %#ok<*FNDSB>
            
        end
        
        if length(n1_mean) > 1
            n1_mean = mean(n1_mean);
        end
        
        n2_mean = zeros(0,length(ERPinfo.stim.(task).n2_ROI));
        
        for jj = 1:length(ERPinfo.stim.(task).n2_ROI)
            
            chnind = find(strcmp(EPdata.chanNames,...
                ERPinfo.stim.(task).n2_ROI{jj})); 
            
            n2_mean(jj) = mean(EPdata.data(chnind,...
                knnsearch(EPdata.timeNames,ERPinfo.stim.(task).n2_meanwind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.stim.(task).n2_meanwind(2)),...
                ii)); %#ok<*FNDSB>
            
        end
        
        if length(n2_mean) > 1
            n2_mean = mean(n2_mean);
        end
        
        if all(size(EPdata.events{1,ii}) == 1)
            sub = EPdata.events{1,ii};
        else
            sub = EPdata.events{1,ii}(find(strcmp({EPdata.events{1,ii}.type}, 'Stm+')== 1));
        end
        
        keys = sub.keys;
        mffkeys = keys(find(strcmp({keys.code}, 'mffkeys')==1)).data;
        mffkeys_parsed = regexp(mffkeys,'\d*','Match');
        
        
        if ~exist('individ','var')
            
            individ = table;
            individ.subjid = cellstr(subjid);
            individ.task = cellstr(task);
            individ.torder = cellstr(torder);
            individ.event = cellstr(event);
            individ.n1 = n1_mean;
            individ.n2 = n2_mean;
%             individ.urtrl = str2double(mffkeys_parsed{6});
            
        elseif exist('individ','var')
            
            row = table;
            row.subjid = cellstr(subjid);
            row.task = cellstr(task);
            row.torder = cellstr(torder);
            row.event = cellstr(event);
            row.n1 = n1_mean;
            row.n2 = n2_mean;
%             row.urtrl = str2double(mffkeys_parsed{6});
            
            individ = vertcat(individ,row); %#ok<AGROW>
            
        end
        
    end
    
end

if ~exist('individ','var')
    individ = [];
end

end


