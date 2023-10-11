%% SCRIPT_analysisLFP
load('data_analysis_aroundCP_22-Sep-2023.mat')

% prep cells - maintain session structure
posAESI{1} = cell(size(posData{1}));
posAESI{2} = cell(size(posData{2}));

lfpAESI{1} = cell(size(lfpData{1}));
lfpAESI{2} = cell(size(lfpData{2}));

delayAESI{1} = cell(size(posData{1}));
delayAESI{2} = cell(size(posData{2}));

% use accBoolean for acc, index trials later

% Clean
% failed choice exit
delayData{2}{6,7}(35)= [];
accBoolean{2}{6,7}(35)= [];
posData{2}{6,7}(35)= [];
lfpData{2}{6,7}(35)= [];
% remove tagged double trials
lfpData{1}{1,20}(17)=[];
posData{1}{1,20}(17)=[];
posData{1}{8,14}(10)=[];

% empty -> nan cells for delay and lfp
for condi= 1:length(posData)
    for rati=1:size(posData{condi},1)
        for sessi=1:size(posData{condi},2)
            posAESI{condi}{rati,sessi}= empty2nan(posData{condi}{rati,sessi});
            delayAESI{condi}{rati,sessi} = empty2nan(delayData{condi}{rati,sessi});
            try  lfpAESI{condi}{rati,sessi} = empty2nan(lfpData{condi}{rati,sessi});
            catch
                continue
            end
        end
    end
end

% Coherence 
% get chronux parameters
params = getCustomParams;
srate  = 2000;
params.Fs = srate;
params.tapers = [2 3];
params.fpass= [0 50];

f = 1:0.5:50;
srate = 2000;
cohAESI = []; powerPFC=[]; powerHPC = [];
for condi = 1:length(lfpAESI)
    for rati = 1:size(lfpAESI{condi},1)
        for sessi= 1:size(lfpAESI{condi},2)
            if isempty(lfpAESI{condi}{rati,sessi})
                continue
            end

            for triali = 1:length(lfpAESI{condi}{rati,sessi})
                if isnan(lfpAESI{condi}{rati,sessi}{triali})
                    cohAESI{condi}{rati,sessi}{triali}= NaN;
                    powerPFC{condi}{rati,sessi}{triali}= NaN;
                    powerHPC{condi}{rati,sessi}{triali}= NaN;
                else
                    pfc = []; hpc = []; lfp = [];
                    lfp = lfpAESI{condi}{rati,sessi}{triali};
                    pfc = detrend(lfp(1,:),3);
                    hpc = detrend(lfp(2,:),3);

                    %welch
                    %get coherence from all trials
                    cohAESI{condi}{rati,sessi}{triali} = mscohere(pfc, hpc, [], [], f, srate);

                    %get power from all trials
                    S1 = []; S2 = [];
                    S1 = pwelch(pfc,[],[],f,srate);
                    S2 = pwelch(hpc,[],[],f,srate);
                    powerPFC{condi}{rati,sessi}{triali} = log10(S1);
                    powerHPC{condi}{rati,sessi}{triali} = log10(S2);

                    %multitaper
                    % [cohAESI{condi}{rati,sessi}{triali},[],powerPFC{condi}{rati,sessi}{triali}, ...
                    %    powerHPC{condi}{rati,sessi}{triali}] =coherencyc(pfc,hpc,params);
                end
            end
        end
        disp(['Completed with ','rat ',num2str(rati),'/',num2str(size(lfpAESI{condi},1))])
    end
end

%get actual session number
sessionNumber{1}= NaN(size(posData{1}));
sessionNumber{2} = NaN(size(posData{2}));
for condi=1:length(sessionName)
    for rati=1:length(sessionName{condi})
        for sessi = 1:length(sessionName{condi}{rati})
            numbers = str2double(extract(sessionName{condi}{rati}{sessi}, digitsPattern)); % get session number
            sessionNumber{condi}(rati, sessi)= numbers(1); 
        end
    end
end

lfpAnalyze{1}=cell(size(lfpData{1},1),max(sessionNumber{1},[],'all'));
lfpAnalyze{2}=cell(size(lfpData{2},1),max(sessionNumber{2},[],'all'));
varNames = [{'Delay'},{'Accuracy'},{'Coherence'},{'PFC Power'}, {'HPC Power'},{'Session'},{'Trial'},{'Position'}];
for condi = 1:length(lfpAESI)
    for rati = 1:size(lfpAESI{condi},1)
        array={};
        templfpAnalyze=[];
        for sessi = 1:size(lfpAESI{condi},2)
            if isempty(lfpAESI{condi}{rati,sessi})
                continue
            end
  
            temp = {};
            temp(:,1) = delayAESI{condi}{rati,sessi}; % get delay
            temp(:,2) = num2cell(accBoolean{condi}{rati,sessi}); %get trial choice accuracy
            temp(:,3) = cohAESI{condi}{rati,sessi}; %get coherence
            temp(:,4) = powerPFC{condi}{rati,sessi};
            temp(:,5) = powerHPC{condi}{rati,sessi};
            temp(1:size(temp,1),6) = num2cell(sessionNumber{condi}(rati,sessi)); % get session/rat indicator 
            temp(:,7) = num2cell(1:size(temp,1)); % trial information
            temp(:,8) = posAESI{condi}{rati,sessi}; %get positon data

            % erase rows w/ NaN
            idxNan= isnan(cell2mat(temp(:,1))) | cell2mat(cellfun(@nnz,cellfun(@isnan,temp(:,8), 'Unif',0), 'Unif',0));
            temp(idxNan,:)=[];

            % keep session structure
            % across rats
            templfpAnalyze{sessi} = vertcat(varNames,temp);
        end

        orderS = sessionNumber{condi}(rati,:); %get session order
        idxnan = isnan(orderS); %remove nonexistent sessions
        orderS(idxnan)=[];

        lfpAnalyze{condi}(rati,orderS)= templfpAnalyze; %get coherence for each session in order
    end
end


%% Position Coherence Analysis 
%get x and y position data from each trial for AE and SI 
xMeans=[]; yMeans=[]; xData=[]; yData=[];
for condi= 1:length(lfpAnalyze)
    for rati= 1:size(lfpAnalyze{condi},1)
        for sessi = 1:length(lfpAnalyze{condi}(rati,:))
            if isempty(lfpAnalyze{condi}{rati,sessi})
                continue
            end

            idxPos = find(contains(lfpAnalyze{condi}{rati,sessi}(1,:),'Position'));
            tempX=[]; tempY=[];
            for triali = 2:length(lfpAnalyze{condi}{rati,sessi})
                try % sometimes srate may have pulled an extra sample
                    % get x/y data
                    tempX(triali,:)= lfpAnalyze{condi}{rati,sessi}{triali,idxPos}(1,:);
                    tempY(triali,:) = lfpAnalyze{condi}{rati,sessi}{triali,idxPos}(2,:);
                end
            end
            % get average x and y within each session
            xMeans{rati}{sessi,1} = mean(tempX,1);
            yMeans{rati}{sessi,1} = mean(tempY,1);
        end
    end
    xData{condi}=emptyCellErase(vertcat(xMeans{:}));
    yData{condi}=emptyCellErase(vertcat(yMeans{:}));
end

% load example position data
datafolder = '/Users/haileyrosenblum/Desktop/Matlab/getVT/DA3_20-53';
cd(datafolder)
[x,y,t] = getVTdata(datafolder,'interp','VT1.mat');

xAE = vertcat(xData{1}{:});
yAE = vertcat(yData{1}{:});
xSI = vertcat(xData{2}{:});
ySI = vertcat(yData{2}{:});

%Plot median position for LFP analysis
figure('color','w')
    plot(x,y,'Color',[.6 .6 .6]); hold on;
    plot(median(xAE,1),median(yAE,1),'r','LineWidth',2)
    plot(median(xSI,1),median(ySI,1),'b','LineWidth',2)
    xlim([150 350])
    ylim([150 350])
    title('Median position for LFP analysis')

%Plot mean position for LFP analysis
figure('Color','w');
    plot(x,y,'Color',[.6 .6 .6]); hold on;
    plot(mean(xAE,1),mean(yAE,1),'r','LineWidth',2)
    plot(mean(xSI,1),mean(ySI,1),'b','LineWidth',2)
    xlim([150 350])
    ylim([150 350])
    title('Mean position for LFP analysis')
    axis off
    
%% Prep Variables for Coherence and Power Analysis AE vs SI 
%get mean coherence, power from each session (1:0.5:50Hz)
cohCat=[]; pfcCat=[]; hpcCat=[]; ratCoh=[];ratPowPFC=[];ratPowHPC=[];sessCoh=[];sessPowPFC=[];sessPowHPC=[];

cohCatShort=[];cohCatMed=[]; cohCatLong=[]; pfcCatShort=[]; pfcCatMed=[];pfcCatLong=[];hpcCatShort=[];hpcCatMed=[]; hpcCatLong=[];

ratCohShort=[]; ratCohMed=[]; ratCohLong=[]; ratPowPfcShort=[]; ratPowPfcMed=[]; ratPowPfcLong=[];ratPowHpcShort=[];ratPowHpcMed=[];ratPowHpcLong=[];
sessCohShort=[]; sessCohMed=[]; sessCohLong=[]; sessPowPfcShort=[]; sessPowPfcMed=[]; sessPowPfcLong=[];sessPowHPCShort=[];sessPowMed=[];sessPowLong=[];

method=1; %decide which delay method to use
f = 1:0.5:50;
for condi = 1:length(lfpAnalyze)
    for rati = 1:size(lfpAnalyze{condi},1)
        for sessi= 1:size(lfpAnalyze{condi},2)
            if isempty(lfpAnalyze{condi}{rati,sessi})
                continue
            end
            idxCoh = find(contains(lfpAnalyze{condi}{rati,sessi}(1,:),'Coherence'));
            idxPFC = find(contains(lfpAnalyze{condi}{rati,sessi}(1,:),'PFC Power'));
            idxHPC = find(contains(lfpAnalyze{condi}{rati,sessi}(1,:),'HPC Power'));
            idxDelay = find(contains(lfpAnalyze{condi}{rati,sessi}(1,:),'Delay'));

            % ACROSS DELAYS
            %get coherence and power for the session
            cohSession=[]; pfcSession=[]; hpcSession=[];
            cohSession = lfpAnalyze{condi}{rati,sessi}(2:end,idxCoh);
            pfcSession = lfpAnalyze{condi}{rati,sessi}(2:end,idxPFC);
            hpcSession = lfpAnalyze{condi}{rati,sessi}(2:end,idxHPC);

            %find session mean
            cohCat{condi}{rati,sessi} = mean(vertcat(cohSession{:}),1);
            pfcCat{condi}{rati,sessi} = mean(vertcat(pfcSession{:}),1);
            hpcCat{condi}{rati,sessi} = mean(vertcat(hpcSession{:}),1);

            % DELAYS
            %get index of all trials corresponding to each delay
            idxShort=[]; idxMed=[]; idxLong=[]; delayTrials=[];
            delayTrials = cell2mat(lfpAnalyze{condi}{rati,sessi}(2:end,idxDelay));
            idxShort = find(delayTrials >= delayMinMax{method}(1,1) & delayTrials <= delayMinMax{method}(1,2))+1;
            idxMed   = find(delayTrials >= delayMinMax{method}(2,1) & delayTrials <= delayMinMax{method}(2,2))+1;
            idxLong  = find(delayTrials >= delayMinMax{method}(3,1) & delayTrials <= delayMinMax{method}(3,2))+1;

            %if there are less than 3 trials at any delay, skip
            if numel(idxShort)<3 || numel(idxMed)<3 || numel(idxLong)<3
                cohCatShort{condi}{rati,sessi} = [];
                cohCatMed{condi}{rati,sessi}   = [];
                cohCatLong{condi}{rati,sessi}  = [];
                pfcCatShort{condi}{rati,sessi} = [];
                pfcCatMed{condi}{rati,sessi}   = [];
                pfcCatLong{condi}{rati,sessi}  = [];
                hpcCatShort{condi}{rati,sessi} = [];
                hpcCatMed{condi}{rati,sessi}   = [];
                hpcCatLong{condi}{rati,sessi}  = [];
                continue;
            end

            %group coherence and power sorted by delay
            shortCoh=[]; medCoh=[]; longCoh=[];
            shortPowPFC=[]; shortPowHPC=[];medPowPFC=[];medPowHPC=[];longPowPFC=[];longPowHPC=[];
            shortCoh = lfpAnalyze{condi}{rati,sessi}(idxShort,idxCoh);
            medCoh   = lfpAnalyze{condi}{rati,sessi}(idxMed, idxCoh);
            longCoh  = lfpAnalyze{condi}{rati,sessi}(idxLong, idxCoh);

            shortPowPFC = lfpAnalyze{condi}{rati,sessi}(idxShort,idxPFC);
            medPowPFC   = lfpAnalyze{condi}{rati,sessi}(idxMed,idxPFC);
            longPowPFC  = lfpAnalyze{condi}{rati,sessi}(idxLong,idxPFC);

            shortPowHPC = lfpAnalyze{condi}{rati,sessi}(idxShort,idxHPC);
            medPowHPC   = lfpAnalyze{condi}{rati,sessi}(idxMed,idxHPC);
            longPowHPC  = lfpAnalyze{condi}{rati,sessi}(idxLong,idxHPC);

            %find session mean
            cohCatShort{condi}{rati, sessi} = mean(vertcat(shortCoh{:}),1);
            cohCatMed{condi}{rati, sessi} = mean(vertcat(medCoh{:}),1);
            cohCatLong{condi}{rati, sessi} = mean(vertcat(longCoh{:}),1);

            pfcCatShort{condi}{rati, sessi} = mean(vertcat(shortPowPFC{:}),1);
            pfcCatMed{condi}{rati, sessi} = mean(vertcat(medPowPFC{:}),1);
            pfcCatLong{condi}{rati, sessi} = mean(vertcat(longPowPFC{:}),1);

            hpcCatShort{condi}{rati, sessi} = mean(vertcat(shortPowHPC{:}),1);
            hpcCatMed{condi}{rati, sessi} = mean(vertcat(medPowHPC{:}),1);
            hpcCatLong{condi}{rati, sessi} = mean(vertcat(longPowHPC{:}),1);
        end

        % RATS
        % OVER DELAYS
        %organize data for rat analysis, get average coherence for each rat
        tempRatCoh=[];tempRatPowPFC=[];tempRatPowHPC=[];
        tempRatCoh    = mean(vertcat(cohCat{condi}{rati,:}),1);
        tempRatPowPFC = mean(vertcat(pfcCat{condi}{rati,:}),1);
        tempRatPowHPC = mean(vertcat(hpcCat{condi}{rati,:}),1);

        %normalize
        ratCoh{condi}(rati,:)    = normalize(tempRatCoh,2,'range');
        ratPowPFC{condi}(rati,:) = normalize(tempRatPowPFC,2,'range');
        ratPowHPC{condi}(rati,:) = normalize(tempRatPowHPC,2,'range');

        % DELAYS
        %organize data for rat analysis, get average coherence for each rat
        tempRatCohShort=[];tempRatCohMed=[];tempRatCohLong=[];
        tempRatPowPfcShort=[]; tempRatPowPfcMed=[]; tempRatPowPfcLong=[];
        tempRatPowHpcShort=[]; tempRatPowHpcMed=[]; tempRatPowHpcLong=[];

        tempRatCohShort  = mean(vertcat(cohCatShort{condi}{rati,:}),1);
        tempRatCohMed  = mean(vertcat(cohCatMed{condi}{rati,:}),1);
        tempRatCohLong  = mean(vertcat(cohCatLong{condi}{rati,:}),1);

        tempRatPowPfcShort = mean(vertcat(pfcCatShort{condi}{rati,:}),1);
        tempRatPowPfcMed = mean(vertcat(pfcCatMed{condi}{rati,:}),1);
        tempRatPowPfcLong = mean(vertcat(pfcCatLong{condi}{rati,:}),1);

        tempRatPowHpcShort = mean(vertcat(hpcCatShort{condi}{rati,:}),1);
        tempRatPowHpcMed = mean(vertcat(hpcCatMed{condi}{rati,:}),1);
        tempRatPowHpcLong = mean(vertcat(hpcCatLong{condi}{rati,:}),1);

        %normalize
        try
            ratCohShort{condi}(rati,:) = normalize(tempRatCohShort,2,'range');
            ratCohMed{condi}(rati,:)   = normalize(tempRatCohMed,2,'range');
            ratCohLong{condi}(rati,:)  = normalize(tempRatCohLong,2,'range');

            ratPowPfcShort{condi}(rati,:) = normalize(tempRatPowPfcShort,2,'range');
            ratPowPfcMed{condi}(rati,:) = normalize(tempRatPowPfcMed,2,'range');
            ratPowPfcLong{condi}(rati,:) = normalize(tempRatPowPfcLong,2,'range');

            ratPowHpcShort{condi}(rati,:) = normalize(tempRatPowHpcShort,2,'range');
            ratPowHpcMed{condi}(rati,:) = normalize(tempRatPowHpcMed,2,'range');
            ratPowHpcLong{condi}(rati,:) = normalize(tempRatPowHpcLong,2,'range');
        catch
            ratCohShort{condi}(rati,1:length(f))=NaN;
            ratCohMed{condi}(rati,1:length(f))=NaN;
            ratCohLong{condi}(rati,1:length(f))=NaN;
            ratPowPfcShort{condi}(rati,1:length(f))=NaN;
            ratPowPfcMed{condi}(rati,1:length(f))=NaN;
            ratPowPfcLong{condi}(rati,1:length(f))=NaN;
            ratPowHpcShort{condi}(rati,1:length(f))=NaN;
            ratPowHpcMed{condi}(rati,1:length(f))=NaN;
            ratPowHpcLong{condi}(rati,1:length(f))=NaN;
        end

    end
    % SESSIONS
    % ACROSS DELAYS
    % organize data for session analysis, noramlize
    sessCoh{condi}= normalize(vertcat(cohCat{condi}{:}),2,'range');
    sessPowPFC{condi}= normalize(vertcat(pfcCat{condi}{:}),2,'range');
    sessPowHPC{condi}= normalize(vertcat(hpcCat{condi}{:}),2,'range');

    %DELAYS
    % organize data for session analysis, noramlize
    sessCohShort{condi}= normalize(vertcat(cohCatShort{condi}{:}),2,'range');
    sessCohMed{condi}= normalize(vertcat(cohCatMed{condi}{:}),2,'range');
    sessCohLong{condi}= normalize(vertcat(cohCatLong{condi}{:}),2,'range');

    sessPowPfcShort{condi}= normalize(vertcat(pfcCatShort{condi}{:}),2,'range');
    sessPowPfcMed{condi}= normalize(vertcat(pfcCatMed{condi}{:}),2,'range');
    sessPowPfcLong{condi}= normalize(vertcat(pfcCatLong{condi}{:}),2,'range');

    sessPowHpcShort{condi}= normalize(vertcat(hpcCatShort{condi}{:}),2,'range');
    sessPowHpcMed{condi}= normalize(vertcat(hpcCatMed{condi}{:}),2,'range');
    sessPowHpcLong{condi}= normalize(vertcat(hpcCatLong{condi}{:}),2,'range');
end

% ACROSS DELAY
% prep for figures
% N= sessions
cohPowAEsess.coh = sessCoh{1};cohPowAEsess.pPFC = sessPowPFC{1}; cohPowAEsess.pHPC = sessPowHPC{1};
cohPowSIsess.coh = sessCoh{2};cohPowSIsess.pPFC = sessPowPFC{2};cohPowSIsess.pHPC = sessPowHPC{2};

% N= rats
cohPowAErat.coh = ratCoh{1};cohPowAErat.pPFC = ratPowPFC{1}; cohPowAErat.pHPC = ratPowHPC{1};
cohPowSIrat.coh = ratCoh{2}; cohPowSIrat.pPFC = ratPowPFC{2}; cohPowSIrat.pHPC = ratPowHPC{2};

% DELAYS
% prep for figures
% N= sessions
    % short
    cohPowAEsessShort.coh = sessCohShort{1}; cohPowAEsessShort.pPFC = sessPowPfcShort{1}; cohPowAEsessShort.pHPC = sessPowHpcShort{1};
    cohPowSIsessShort.coh = sessCohShort{2}; cohPowSIsessShort.pPFC = sessPowPfcShort{2};cohPowSIsessShort.pHPC = sessPowHpcShort{2};

    % medium
    cohPowAEsessMed.coh = sessCohMed{1}; cohPowAEsessMed.pPFC = sessPowPfcMed{1}; cohPowAEsessMed.pHPC = sessPowHpcMed{1};
    cohPowSIsessMed.coh = sessCohMed{2}; cohPowSIsessMed.pPFC = sessPowPfcMed{2};cohPowSIsessMed.pHPC = sessPowHpcMed{2};
    
    % long
    cohPowAEsessLong.coh = sessCohLong{1}; cohPowAEsessLong.pPFC = sessPowPfcLong{1}; cohPowAEsessLong.pHPC = sessPowHpcLong{1};
    cohPowSIsessLong.coh = sessCohLong{2}; cohPowSIsessLong.pPFC = sessPowPfcLong{2};cohPowSIsessLong.pHPC = sessPowHpcLong{2};

% N= rats
    %short
    cohPowAEratShort.coh = ratCohShort{1}; cohPowAEratShort.pPFC = ratPowPfcShort{1}; cohPowAEratShort.pHPC = ratPowHpcShort{1};
    cohPowSIratShort.coh = ratCohShort{2}; cohPowSIratShort.pPFC = ratPowPfcShort{2};cohPowSIratShort.pHPC = ratPowHpcShort{2};

    % medium
    cohPowAEratMed.coh = ratCohMed{1}; cohPowAEratMed.pPFC = ratPowPfcMed{1}; cohPowAEratMed.pHPC = ratPowHpcMed{1};
    cohPowSIratMed.coh = ratCohMed{2}; cohPowSIratMed.pPFC = ratPowPfcMed{2};cohPowSIratMed.pHPC = ratPowHpcMed{2};
    
    % long
    cohPowAEratLong.coh = ratCohLong{1}; cohPowAEratLong.pPFC = ratPowPfcLong{1}; cohPowAEratLong.pHPC = ratPowHpcLong{1};
    cohPowSIratLong.coh = ratCohLong{2}; cohPowSIratLong.pPFC = ratPowPfcLong{2};cohPowSIratLong.pHPC = ratPowHpcLong{2};


%% Coherence and Power Across Delays (N= Rats and Sessions)
%sessions
f = 1:0.5:50;
titleCoh = 'N=sessions'; titlePow = 'N=sessions';
cohPowFigures(cohPowAEsess, cohPowSIsess, f, titleCoh, titlePow)

%rats
titleCoh = 'N=rats'; titlePow = 'N=rats';
cohPowFigures(cohPowAErat, cohPowSIrat, f, titleCoh, titlePow)

%% Coherence AE vs SI (by delay)
%Short
    % sessions
    titleCoh = 'Short Delay N=sessions'; titlePow = 'Short Delay N=sessions';
    cohPowFigures(cohPowAEsessShort, cohPowSIsessShort, f, titleCoh, titlePow)

    %rats
    titleCoh = 'Short Delay N=rats'; titlePow = 'Short Delay N=rats';
    cohPowFigures(cohPowAEratShort, cohPowSIratShort, f, titleCoh, titlePow)
%Med
    % sesssions
    titleCoh = 'Medium Delay N=sessions'; titlePow = 'Medium Delay N=sessions';
    cohPowFigures(cohPowAEsessMed, cohPowSIsessMed, f, titleCoh, titlePow)

    %rats
    titleCoh = 'Medium Delay N=rats'; titlePow = 'Medium Delay N=rats';
    cohPowFigures(cohPowAEratMed, cohPowSIratMed, f, titleCoh, titlePow)
%Long
   % sessions
    titleCoh = 'Long Delay N=sessions'; titlePow = 'Long Delay N=sessions';
    cohPowFigures(cohPowAEsessLong, cohPowSIsessLong, f, titleCoh, titlePow)

    %rats
    titleCoh = 'Long Delay N=rats'; titlePow = 'Long Delay N=rats';
    cohPowFigures(cohPowAEratLong, cohPowSIratLong, f, titleCoh, titlePow)


