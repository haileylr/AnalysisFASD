%% Format data files
% get LFP and position data -0.5 - +1s relative to choice point entry
% also save out trial accuracy, session choice accuraccy, trial delay
% length, session/rat number and time spent in the choice point
 
clear;
addpath(getCurrentPath);
load('data_directories')
 
% get data directory
codeDir = getCurrentPath;
codePiece = strsplit(codeDir,'\','CollapseDelimiters',false);
codePiece(end-1:end)=[];
dataDir = join(codePiece,'\');
dataDir = dataDir{1};
dataDir = horzcat(dataDir,'\data');
 
% conditions
conditions{1} = 'AE';
conditions{2} = 'SI';
 
getNewData = 1;
if getNewData == 1
     % this script is to prep your paths for mvgc
    for condi = 1:length(conditions)
 
        % get datafolders (session names) into a cell array
        Datafolders = [dataDir,'\',conditions{condi}];
        cd(Datafolders)
        
        % some stuff to make sure we're extracting what we think and also to
        % get our rat names out
        dir_content = [];
        dir_content = dir(Datafolders);
        if isempty(dir_content)
            continue
        end
        dir_content = extractfield(dir_content,'name');
        remIdx = contains(dir_content,'.mat') | contains(dir_content,'.');
        dir_content(remIdx)=[];
        ratIDs = dir_content;
 
        for rati = 1:length(ratIDs)
 
            % create datafolder per rat
            ratDir = [];
            ratDir = [Datafolders,'\',ratIDs{rati},'\DA'];
            cd(ratDir);
 
            % get session datafolder
            dir_content = [];
            dir_content = dir(ratDir);
            if isempty(dir_content)
                continue
            end
            dir_content = extractfield(dir_content,'name');
            remIdx = contains(dir_content,'.mat') | contains(dir_content,'.');
            dir_content(remIdx)=[];
            % get sessions
            sessions = dir_content;
            sessionName{condi}{rati} = sessions; %update 9/20 -- need to order sessions
 
            for sessi = 1:length(sessions)
 
                try
                    % create datafolder
                    datafolder = [Datafolders,'\',ratIDs{rati},'\DA\',sessions{sessi}];
                    cd(datafolder);
 
                    % load sequenceArray
                    load('sequenceArray');
 
                    % first calculate timespent at choice - remove trials with failed choice exit/ greater than 10% tracking error on stem
                    tempfailedTS = []; failedTS = []; choiceExitEntry = [];
                    idxExit  = find(contains(sequenceCell(1,:),'gaEntry'));
                    idxEntry = find(contains(sequenceCell(1,:),'cpEntry'));
                    idxFailedTS = find(contains(sequenceCell(1,:),'trackingErrorStem'));
                    choiceExitEntry = [sequenceCell(3:end,idxExit), sequenceCell(3:end,idxEntry)];
                    
                    % for double trial sessions, will need to remove empty
                    % cell or will get error message
                    [tempfailedTS, idxRem] = emptyCellErase(sequenceCell(3:end,idxFailedTS)); 
                    choiceExitEntry(idxRem,:) = [];
                    
                    % index trials to remove from mean
                    failedTS = cell2mat(tempfailedTS)==1;
                    choiceExitEntry(failedTS,:) = [];
                    tsChoice{condi}(rati,sessi) = nanmean((cell2mat(choiceExitEntry(:,1))-cell2mat(choiceExitEntry(:,2)))./1e6);
 
                    % get lfp and position from choice point 
                    % format variables for looping
                    lfp = cell([size(sequenceCell,1)-2 1]);
                    pos = cell([size(sequenceCell,1)-2 1]);
                    numTrials = size(sequenceCell,1)-2;
                    for triali = 3:numTrials+2
                        %find cp entry time 
                        cpEntryIdx = []; cpTime = [];
                        cpEntryIdx = find(contains(sequenceCell(1,:),'cpEntry'));
                        cpTime     = [sequenceCell{triali,cpEntryIdx}]; 
                        
                        % get position data-- sequence cell has whole trial
                        idxPOS  = []; temppos = [];
                        idxPOS  = find(contains(sequenceCell(1,:),'trialPosition'));
                        temppos = sequenceCell{triali,idxPOS};
                        
                        %skip previously excluded trials
                        if isempty(temppos)
                            continue
                        end
                        
                        %position -0.5 to +1s relative to cp entry
                        timepos = []; aroundCP2 = []; posAroundCP = [];
                        timepos     = temppos(3,:);
                        aroundCP2   = dsearchn(timepos', cpTime');
                        posAroundCP = temppos(:, aroundCP2-(30*0.5):aroundCP2+(30*1));
                        
                        %getting position around cp from all trials
                        pos{triali-2} = posAroundCP;
                        
                        %get lfp from rats used for lfp analysis (behavior
                        %only rats will not have LFP column in sequence cell)
                        idxLFP = [];
                        idxLFP  = find(contains(sequenceCell(1,:),'trialLFP'));
                        if ~isempty(idxLFP)  %for rats with LFP
                            %get trial lfp-- sequence cell has whole trial
                            templfp = [];
                            templfp = sequenceCell{triali,idxLFP};
                            
                            %lfp -0.5 to +1s relative to cp entry
                            timelfp = []; aroundCP1 = []; lfpAroundCP = [];
                            timelfp     = templfp(3,:);
                            aroundCP1   = dsearchn(timelfp', cpTime');
                            lfpAroundCP = templfp(:, aroundCP1-(2000*0.5):aroundCP1+(2000*1));
                            
                            %getting lfp (if exitsts) around cp from all trials
                            lfp{triali-2} = lfpAroundCP;
                        end
                    end
                    
                    if ~isempty(idxLFP)
                        % exclude lfp and corresponding position trial data if lfp trial
                        % was marked for removal 
                        tempexclude = []; idxEX = [];
                        idxEX       = find(contains(sequenceCell(1,:),'excludeLFP_stem2choice'));
                        tempexclude = cell2mat(sequenceCell(3:end,idxEX));
                        
                        extrial = [];
                        extrial = tempexclude(:,1)==1 | tempexclude(:,2)==1;
                        for triali = 1:length(extrial)
                            if extrial(triali)==1
                                lfp{triali} = [];
                                pos{triali} = [];
                            end
                        end
                        lfpData{condi}{rati,sessi} = lfp; % store lfp (only from rats with lfp)
                     end
 
                    posData{condi}{rati,sessi} = pos; % store position
                    
                    % get accuracy data
                    tempAcc = []; tempAcc = cell2mat(sequenceCell(3:end,contains(sequenceCell(1,:),'accBinary')));
                    accData{condi}{rati,sessi} = (1-mean(tempAcc))*100; % store session accuracy
                    
                   % get and store delay data
                    tempDelay = [];
                    tempDelay = sequenceCell(3:end,contains(sequenceCell(1,:),'delayTime'));
                    % empty cell for trials w/ failed sb entry or exit
                    idxFailedExit  = find(contains(sequenceCell(1,:),'trackingErrorStem'));
                    idxFailedEntry = find(contains(sequenceCell(1,:),'trackingErrorSBentry'));
                    
                    %remove any empty cells due to double trials, then
                    %remove trials with delay errors
                    tempexlude = []; extrial = []; idxRem = [];
                    [failedExit,idxRem] = emptyCellErase(sequenceCell(3:end,idxFailedExit));
                    failedEntry = sequenceCell(3:end,idxFailedEntry);
                    failedEntry(idxRem)= [];
                    tempDelay(idxRem) = [];
 
                    tempexclude = cell2mat([failedEntry, failedExit]);
                    extrial = tempexclude(:,1)==1 | tempexclude(:,2)==1;
                    for triali = 1:length(extrial)
                        if extrial(triali)==1
                            tempDelay{triali} = [];
                        end
                    end
                    delayData{condi}{rati, sessi} = tempDelay; % store delay
                    
                    % get and store accuracy (boolean) data to separate
                    % correct/errors later if so desired
                    tempBoolean = [];
                    tempBoolean = cell2mat(sequenceCell(3:end,contains(sequenceCell(1,:),'accBinary')));
                    accBoolean{condi}{rati, sessi} = tempBoolean;
                    
                catch
                    disp('Could not extract lfp and position data');
                    pause;
                end
            end
            disp(['Completed ',ratIDs{rati}]);
        end   
        disp(['Completed with ',conditions{condi},' rat ',ratIDs{rati}])
    end
    cd(paths.dataOUTpath);
    date = char(datetime('today'));
    save(['data_analysis_aroundCP_',date],'lfpData','posData','tsChoice','accData', 'delayData', 'accBoolean', 'sessionName');
else
    load('data_analysis_15-Jun-2023'); %file no longer exists 
end
