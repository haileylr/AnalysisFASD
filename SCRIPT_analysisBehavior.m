%% This script includes analyses on time spent in choice point, delay length, and choice accuracy 
dataDirectory = '/Users/haileyrosenblum/Desktop/Matlab/Data output';
cd(dataDirectory)
load("data_analysis_aroundCP_22-Sep-2023.mat")

% remove SI rat 6 session 7 trial 35 - failed choice exit 
delayData{2}{6,7}(35)= [];
accBoolean{2}{6,7}(35)= [];

%clean data
accAE   = accData{1}(:);
accSI   = accData{2}(:);
tsChoiceSI = tsChoice{2}(:);
tsChoiceAE = tsChoice{1}(:);

[~,idxAE] = emptyCellErase(accAE); 
[~,idxSI] = emptyCellErase(accSI); 

accAE(idxAE)=[];  tsChoiceAE(idxAE)=[];
accSI(idxSI)=[];  tsChoiceSI(idxSI)=[];


%remove empty trials
accB = accBoolean;

delayAeSi=[]; 
delayAeSi{1} = cell(size(delayData{1}));
delayAeSi{2} = cell(size(delayData{2}));
for condi= 1:length(delayData)
    for rati=1:size(delayData{condi},1)
        for sessi=1:size(delayData{condi},2)
            if isempty(delayData{condi}{rati,sessi})
                continue
            end
            idxRem = [];
            [delayAeSi{condi}{rati,sessi}, idxRem] = emptyCellErase(delayData{condi}{rati,sessi});
            accB{condi}{rati,sessi}(idxRem) = [];
        end
    end
end


%% Time Spent in Choice Point
tsChoiceAE_rat = []; tsChoiceSI_rat = [];
tsChoiceAE_rat = tsChoice{1};
idx0 = find(tsChoiceAE_rat ==0);
tsChoiceAE_rat(idx0)=NaN;
tsChoiceAE_rat(8,:)=[]; %remove the outlier

idx0= [];
tsChoiceSI_rat = tsChoice{2};
idx0 = find(tsChoiceSI_rat ==0);
tsChoiceSI_rat(idx0)=NaN;

% BY RAT
tsAE_rat = []; tsSI_rat = [];
tsAE_rat = mean(tsChoiceAE_rat, 2, 'omitnan');
tsSI_rat = mean(tsChoiceSI_rat, 2, 'omitnan');

isoutlier(tsSI_rat);
isoutlier(tsAE_rat); %if included rat 8 is an outlier

data2plot = []; data2plot{1} = tsSI_rat; data2plot{2} = tsAE_rat;
figure('color','w')
multiBarPlot(data2plot,[{'SI'} {'AE'}],'Time spent in choice point (s)','y')
title('N=rats')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'TS at choice',numCorrections);

% BY SESSION
%testing sessions with outlier rat removed
adtest(vertcat(tsChoiceSI_rat(:))) %not normally distributed --> nonparametric test 
adtest(vertcat(tsChoiceAE_rat(:)))

data2plot = []; data2plot{1} = tsChoiceSI_rat(:); data2plot{2} = tsChoiceAE_rat(:);
figure('color','w')
multiBarPlot(data2plot,[{'SI'} {'AE'}],'Time spent in choice point (s)','vertical','n')
title('N=sessions')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'TS at choice',numCorrections);
[p,h,stats] = ranksum(tsChoiceAE,tsChoiceSI); 


%% Overall Choice Accuracy
%RAT 
%get index to return to rat
accAE = cell2mat(empty2nan(accData{1}));
accAE(8,:)=[]; %remove outlier rat 8
accSI = cell2mat(empty2nan(accData{2}));

accAE_rat = []; accSI_rat = [];
accAE_rat = mean(accAE, 2, 'omitnan');
accSI_rat = mean(accSI, 2, 'omitnan');

data2plot = []; data2plot{1} = accSI_rat; data2plot{2} = accAE_rat;
figure('color','w')
multiBarPlot(data2plot,[{'SI'} {'AE'}],'Overall Choice Accuracy','y')
title('N=rats')
ylim([50,90])
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Overall Choice Accuracy',numCorrections);

%SESSIONS
data2plot = []; data2plot{1} = accSI(:); data2plot{2} = accAE(:);
figure('color','w')
multiBarPlot(data2plot,[{'SI'} {'AE'}],'Overall Choice Accuracy','n')
title('N=rats')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
title('N=sessions')
ylim([50,90])
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Overall Choice Accuracy',numCorrections);


%OR based on all sessions (including sessions exluded due to poor
%tracking/signal)
%RATS
load('allAccData')
allAccData{1}(8,:)=[]; %remove outlier rat
data2plot = []; data2plot{1} = mean(allAccData{1},2,'omitnan'); data2plot{2} = mean(allAccData{2},2,'omitnan');
figure('color','w')
multiBarPlot(data2plot,[{'AE'} {'SI'}],'Overall Choice Accuracy (%)','y')
title('N=rats')
ylim([50,90])
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Overall Choice Accuracy',numCorrections);

%SESSIONS
data2plot = []; data2plot{1} = allAccData{1}(:); data2plot{2} = allAccData{2}(:);
figure('color','w')
multiBarPlot(data2plot,[{'AE'} {'SI'}],'Overall Choice Accuracy (%)','n')
title('N=rats')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
title('N=sessions')
ylim([50,90])
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Overall Choice Accuracy',numCorrections);


%% Proportion of "poor" sessions and Days to 80%
%order session choice accuracy, fill in blanks with behavior sheets
accOrdered{1}=cell(11,22);
accOrdered{2}=cell(7,18);
for condi=1:length(accData)
    for rati=1:size(accData{condi},1)
        orderS=[]; 
        orderS = sessionNumber{condi}(rati,:); %get session order
        idxnan = isnan(orderS); %remove nonexistent sessions
        orderS(idxnan)=[];

        accSess =emptyCellErase(accData{condi}(rati,:)); %remove nonexistent sessions
        accOrdered{condi}(rati,orderS)= accSess;
    end
end

%move to excel and fill in blanks from behavior sheets
AE_allSessions=readmatrix("AE_allSessions.xlsx");
SI_allSessions=readmatrix("SI_allSessions.xlsx");
allAccData= [];
allAccData{1}=AE_allSessions;
allAccData{2}=SI_allSessions;

%all rats have data up to DA9
poorSess=[]; days2_80=[];
for condi=1:length(allAccData)
    for rati= 1:size(allAccData{condi},1)
        poorSess{condi}(rati) = nnz(allAccData{condi}(rati,:)<80)/nnz(~isnan(allAccData{condi}(rati,:))); %prop below 80
        days2_80{condi}(rati)= find(allAccData{condi}(rati,:)>=80,1); %days to reach 80
    end
end

%exclude rat 8
poorSess{1}(8)=[];
days2_80{1}(8)=[];

%plot poor sessions (if want equal number of sessions for each rat change
%to (rati,1:9)
data2plot = []; data2plot{1} = poorSess{2}; data2plot{2} = poorSess{1};
figure('color','w')
multiBarPlot(data2plot,[{'SI'} {'AE'}],'Proportion of sessions below 80% choice accuracy','y')
title('N=rats')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Proportion of sessions below 80% choice accuracy',numCorrections);

%plot days to 80
data2plot = []; data2plot{1} = days2_80{2}; data2plot{2} = days2_80{1};
figure('color','w')
multiBarPlot(data2plot,[{'SI'} {'AE'}],'Days to reach 80%','y')
title('N=rats')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Days to reach 80%',numCorrections);
ylim([1 13])


%% Method 1 (more restrictive delay definition) Choice Accuracy by Delay
%KMEANS -- all delays for AE and SI rats 
delayTimes=[];  cluster=[]; idxDelay=[]; delays=[];ogOrder=[];
for condi= 1:length(delayData)
    for rati= 1:size(delayData{condi},1)
    tempDelayTimes=[]; delayTimes=[];
    tempDelayTimes= cell2mat(vertcat(delayData{condi}{rati,:}));

    %remove delays above 80s
    delayTimes = tempDelayTimes(tempDelayTimes <85);

    %kmeans- cluster 3 delay groups
    idx= []; C = [];
    rng('default');
    [idx,C] = kmeans(delayTimes,3);
    
    %save out order before sorting...clusters will not be assigned in
    %ascending order automatically
    [cluster{condi}(rati,:), ogOrder{condi}(rati,:)]= sort(C);
    idxDelay{condi}{rati}=idx;

    %also store delay
    delays{condi}{rati}=delayTimes;
    end
end


%return cluster data --- but clusters may not be in correct order --- assign? S,M,L
clusterDelay=[]; minTrials=[];
for condi=1:length(cluster)
    for rati=1:length(cluster{condi})
        %get cluster index of short, medium, long for each rat
        short=[]; med=[]; long=[];
        short = ogOrder{condi}(rati,1);
        med   = ogOrder{condi}(rati,2);
        long  = ogOrder{condi}(rati,3);
        
        % idx of trials corresponding to short, medium, long clusters
        returnShort=[]; returnMed=[]; returnLong=[];
        newArray = zeros(length(idxDelay{condi}{rati}), 1);
        returnShort = find(idxDelay{condi}{rati}==short);
        returnMed   = find(idxDelay{condi}{rati}==med);
        returnLong  = find(idxDelay{condi}{rati}==long);
        
        %assign trial clusters (1=short, 2=medium, 3=long)
        newArray(returnShort) = 1;
        newArray(returnMed)   = 2;
        newArray(returnLong)  = 3;

        %get number of trials in the smallest cluster for each rat
        minTrials{condi}{rati} = min([length(find((newArray==1))), length(find((newArray==2))), length(find((newArray==3)))]);

        clusterDelay{condi}{rati} = newArray;
    end
end

%find rat with the lowest number of any trial type-- then half that to get
%the minimum sample size for random sampling iterations
sampleSz = min([min(cell2mat(minTrials{1})), min(cell2mat(minTrials{1}))])/2;


%separate delay clusters for each rat
c1= [];c2 =[];c3= [];cluster1 =[];cluster2 =[];cluster3 =[];
for condi=1:length(clusterDelay)
    for rati=1:length(clusterDelay{condi})
        idxShort =[]; idxMed =[]; idxLong =[];
        idxShort = find(clusterDelay{condi}{rati}==1);
        idxMed = find(clusterDelay{condi}{rati}==2);
        idxLong = find(clusterDelay{condi}{rati}==3);
        
        c1{condi}{rati}= delays{condi}{rati}(idxShort);
        c2{condi}{rati} = delays{condi}{rati}(idxMed);
        c3{condi}{rati} = delays{condi}{rati}(idxLong);
    end
end
cluster1 = horzcat(c1{:});
cluster2 = horzcat(c2{:});
cluster3 = horzcat(c3{:});

% random sampling
thMinS = []; thMaxS = []; thMinM = []; thMaxM = []; thMinL = []; thMaxL = []; 
shortMin =[]; shortMax = []; medMin = []; medMax = []; longMin = []; longMax = [];
rng('default');
for n=1:1000
    randS = []; randM =[]; randL = [];
    for rati = 1:length(cluster1)
        % random index
        idxRandS = randperm(length(cluster1{rati}),sampleSz);
        idxRandM = randperm(length(cluster2{rati}),sampleSz);
        idxRandL = randperm(length(cluster3{rati}),sampleSz);

        %use random index to get delays
        randS{rati} = cluster1{rati}(idxRandS);
        randM{rati} = cluster2{rati}(idxRandM);
        randL{rati} = cluster3{rati}(idxRandL);
    end
    
    % made vector of all delays in each randomly sampled cluster (across all rats)
    randS_cat = []; randM_cat = []; randL_cat = [];
    randS_cat = vertcat(randS{:});
    randM_cat = vertcat(randM{:});
    randL_cat = vertcat(randL{:});

    %zscore -- test different values
    z = 0.5:0.25:4;
    for i = 1:length(z)
        randS_z = []; randM_z = []; randL_z = [];
        randS_z = zscore(log(randS_cat));
        randM_z = zscore(log(randM_cat));
        randL_z = zscore(log(randL_cat));

        % find z score +/- n sd
        % note to self: 0.5 is too restrictive--way too much data being exluded
        thMinS(n, i) = randS_cat(dsearchn(randS_z, -z(i)));
        thMaxS(n, i) = randS_cat(dsearchn(randS_z, z(i)));

        thMinM(n, i) = randM_cat(dsearchn(randM_z, -z(i)));
        thMaxM(n, i) = randM_cat(dsearchn(randM_z, z(i)));

        thMinL(n, i) = randL_cat(dsearchn(randL_z, -z(i)));
        thMaxL(n, i) = randL_cat(dsearchn(randL_z, z(i)));
    end
end

%get mean thresholds at each z 
shortMin = mean(thMinS,1);
shortMax = mean(thMaxS,1);

medMin = mean(thMinM,1);
medMax = mean(thMaxM,1);

longMin = mean(thMinL,1);
longMax = mean(thMaxL,1);

% idx sd of interest -- here we are lookinng at +/-0.75 sd (reference delayThresholds SCRIPT)
idxMinMax = find(z== 0.75);

%prep cell to save delay info
delayMinMax=[];
delayMinMax{1}= [shortMin(idxMinMax), shortMax(idxMinMax); medMin(idxMinMax), medMax(idxMinMax); longMin(idxMinMax), longMax(idxMinMax)];


%prep cells
delayShort{1}  = nan(size(delayData{1}));
delayShort{2}  = nan(size(delayData{2}));
delayMedium{1} = nan(size(delayData{1}));
delayMedium{2} = nan(size(delayData{2}));
delayLong{1}   = nan(size(delayData{1}));
delayLong{2}   = nan(size(delayData{2}));

accShort{1} = nan(size(delayData{1}));
accShort{2} = nan(size(delayData{2}));
accMed{1}   = nan(size(delayData{1}));
accMed{2}   = nan(size(delayData{2}));
accLong{1}  = nan(size(delayData{1}));
accLong{2}  = nan(size(delayData{2}));

for condi = 1:length(delayAeSi)
    for rati= 1:size(delayAeSi{condi},1)
        for sessi= 1:length(delayAeSi{condi}(rati,:))
              if isempty(delayAeSi{condi}{rati,sessi})
                  continue
              end
               
              trialDelays=[]; tempAcc=[];
              trialDelays = cell2mat(delayAeSi{condi}{rati,sessi});
              tempAcc = accB{condi}{rati,sessi};

              %index short medium long delays
              tempShort= []; tempMed= []; tempLong= [];
              tempShort = find(shortMin(idxMinMax) < trialDelays & trialDelays < shortMax(idxMinMax));
              tempMed   = find(medMin(idxMinMax) < trialDelays & trialDelays < medMax(idxMinMax));
              tempLong  = find(longMin(idxMinMax) < trialDelays & trialDelays < longMax(idxMinMax));
              
              %skip session if there aren't at least 3 trials of each delay
              %length
              if length(tempShort)< 3 || length(tempMed)<3 || length(tempLong)<3
                  continue
              end

              %Find mean accuracy for short medium long delays for each
              %session and store
              accShort{condi}(rati,sessi) = (1-nanmean(tempAcc(tempShort)))*100;
              accMed{condi}(rati,sessi)   = (1-nanmean(tempAcc(tempMed)))*100;
              accLong{condi}(rati,sessi)  = (1-nanmean(tempAcc(tempLong)))*100;
        end
    end
end


%SESSIONS
%remove outlier AE rat 8
accShort_or = accShort{1};
accShort_or(8,:)=[];
accMed_or = accMed{1};
accMed_or(8,:)=[];
accLong_or = accLong{1};
accLong_or(8,:)=[];

%remove nan
aShortSI = rmmissing(accShort{2}(:));
aShortAE = rmmissing(accShort_or(:));
aMedSI = rmmissing(accMed{2}(:));
aMedAE = rmmissing(accMed_or(:));
aLongSI = rmmissing(accLong{2}(:));
aLongAE = rmmissing(accLong_or(:));

% Accuracy Error Bar plot
meanAccAE = [mean(aShortAE), mean(aMedAE), mean(aLongAE)];
meanAccSI = [mean(aShortSI), mean(aMedSI), mean(aLongSI)];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(aShortAE,1), stderr(aMedAE,1), stderr(aLongAE,1)];
errSI = [stderr(aShortSI,1), stderr(aMedSI,1), stderr(aLongSI,1)];
figure('Color', 'w'); hold on;
errorbar(group, meanAccAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanAccSI,errSI, 'b', 'LineWidth', 2);
ylabel('Choice Accuracy')
title('Method 1, N=sessions')

%ttest
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(aShortAE,aShortSI,parametric,stat_test,'Choice Accuracy Short Delay',numCorrections);
readStats(aMedAE,aMedSI,parametric,stat_test,'Choice Accuracy Medium Delay',numCorrections);
readStats(aLongAE,aLongSI,parametric,stat_test,'Choice Accuracy Long Delay',numCorrections);

%RATS
%Acc
accS_AErat = mean(accShort_or,2, 'omitnan');
accS_SIrat = mean(accShort{2},2, 'omitnan');

accM_AErat = mean(accMed_or,2, 'omitnan');
accM_SIrat = mean(accMed{2},2, 'omitnan');

accL_AErat = mean(accLong_or,2, 'omitnan');
accL_SIrat = mean(accLong{2},2, 'omitnan');

% Accuracy Error Bar plot
meanAccAE = [mean(accS_AErat,'omitnan'), mean(accM_AErat,'omitnan'), mean(accL_AErat,'omitnan')];
meanAccSI = [mean(accS_SIrat), mean(accM_SIrat), mean(accL_SIrat)];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(accS_AErat,1), stderr(accM_AErat,1), stderr(accL_AErat,1)];
errSI = [stderr(accS_SIrat,1), stderr(accM_SIrat,1), stderr(accL_SIrat,1)];
figure('Color', 'w'); hold on;
errorbar(group, meanAccAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanAccSI,errSI, 'b', 'LineWidth', 2);
ylabel('Choice Accuracy')
title('Method 1, N=rats')

%ttest
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(accS_SIrat, accS_AErat,parametric,stat_test,'Choice Accuracy Short Delay',numCorrections);
readStats(accM_SIrat, accM_AErat,parametric,stat_test,'Choice Accuracy Medium Delay',numCorrections);
readStats(accL_SIrat, accL_AErat,parametric,stat_test,'Choice Accuracy Long Delay',numCorrections);


%% Method 2 (less restrictive delay definition) Choice Accuracy by Delay
% idx sd of interest -- here we are lookinng at +/-2.25 sd (reference delayThresholds SCRIPT)
idxMinMax = find(z== 2.25);

%add second method to delay data cell and save
delayMinMax{2}= [10, shortMax(idxMinMax); 30, medMax(idxMinMax); 60, longMax(idxMinMax)];
save('delayMinMax', 'delayMinMax')

%prep cells
delayShort{1}  = nan(size(delayData{1}));
delayShort{2}  = nan(size(delayData{2}));
delayMedium{1} = nan(size(delayData{1}));
delayMedium{2} = nan(size(delayData{2}));
delayLong{1}   = nan(size(delayData{1}));
delayLong{2}   = nan(size(delayData{2}));

accShort{1} = nan(size(delayData{1}));
accShort{2} = nan(size(delayData{2}));
accMed{1}   = nan(size(delayData{1}));
accMed{2}   = nan(size(delayData{2}));
accLong{1}  = nan(size(delayData{1}));
accLong{2}  = nan(size(delayData{2}));

for condi = 1:length(delayAeSi)
    for rati= 1:size(delayAeSi{condi},1)
        for sessi= 1:length(delayAeSi{condi}(rati,:))
              if isempty(delayAeSi{condi}{rati,sessi})
                  continue
              end
               
              trialDelays=[]; tempAcc=[];
              trialDelays = cell2mat(delayAeSi{condi}{rati,sessi});
              tempAcc = accB{condi}{rati,sessi};

              %index short medium long delays
              tempShort= []; tempMed= []; tempLong= [];
              tempShort = find(10 < trialDelays & trialDelays < shortMax(idxMinMax));
              tempMed   = find(30 < trialDelays & trialDelays < medMax(idxMinMax));
              tempLong  = find(60 < trialDelays & trialDelays < longMax(idxMinMax));
              
              %skip session if there aren't at least 3 trials of each delay
              %length
              if length(tempShort)< 3 || length(tempMed)<3 || length(tempLong)<3
                  continue
              end

              %Find mean accuracy for short medium long delays for each
              %session and store
              accShort{condi}(rati,sessi) = (1-nanmean(tempAcc(tempShort)))*100;
              accMed{condi}(rati,sessi)   = (1-nanmean(tempAcc(tempMed)))*100;
              accLong{condi}(rati,sessi)  = (1-nanmean(tempAcc(tempLong)))*100;
        end
    end
end


%SESSIONS
%remove outlier AE rat 8
accShort_or = accShort{1};
accShort_or(8,:)=[];
accMed_or = accMed{1};
accMed_or(8,:)=[];
accLong_or = accLong{1};
accLong_or(8,:)=[];

%remove nan
aShortSI = rmmissing(accShort{2}(:));
aShortAE = rmmissing(accShort_or(:));
aMedSI = rmmissing(accMed{2}(:));
aMedAE = rmmissing(accMed_or(:));
aLongSI = rmmissing(accLong{2}(:));
aLongAE = rmmissing(accLong_or(:));

% Accuracy Error Bar plot
meanAccAE = [mean(aShortAE), mean(aMedAE), mean(aLongAE)];
meanAccSI = [mean(aShortSI), mean(aMedSI), mean(aLongSI)];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(aShortAE,1), stderr(aMedAE,1), stderr(aLongAE,1)];
errSI = [stderr(aShortSI,1), stderr(aMedSI,1), stderr(aLongSI,1)];
figure('Color', 'w'); hold on;
errorbar(group, meanAccAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanAccSI,errSI, 'b', 'LineWidth', 2);
ylabel('Choice Accuracy')
title('Method 2, N=sessions')

%ttest
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(aShortAE,aShortSI,parametric,stat_test,'Choice Accuracy Short Delay',numCorrections);
readStats(aMedAE,aMedSI,parametric,stat_test,'Choice Accuracy Medium Delay',numCorrections);
readStats(aLongAE,aLongSI,parametric,stat_test,'Choice Accuracy Long Delay',numCorrections);

%RATS
%Acc
accS_AErat = mean(accShort_or,2, 'omitnan');
accS_SIrat = mean(accShort{2},2, 'omitnan');

accM_AErat = mean(accMed_or,2, 'omitnan');
accM_SIrat = mean(accMed{2},2, 'omitnan');

accL_AErat = mean(accLong_or,2, 'omitnan');
accL_SIrat = mean(accLong{2},2, 'omitnan');

% Accuracy Error Bar plot
meanAccAE = [mean(accS_AErat,'omitnan'), mean(accM_AErat,'omitnan'), mean(accL_AErat,'omitnan')];
meanAccSI = [mean(accS_SIrat), mean(accM_SIrat), mean(accL_SIrat)];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(accS_AErat,1), stderr(accM_AErat,1), stderr(accL_AErat,1)];
errSI = [stderr(accS_SIrat,1), stderr(accM_SIrat,1), stderr(accL_SIrat,1)];
figure('Color', 'w'); hold on;
errorbar(group, meanAccAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanAccSI,errSI, 'b', 'LineWidth', 2);
ylabel('Choice Accuracy')
title('Method 2, N=rats')

%ttest
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(accS_SIrat, accS_AErat,parametric,stat_test,'Choice Accuracy Short Delay',numCorrections);
readStats(accM_SIrat, accM_AErat,parametric,stat_test,'Choice Accuracy Medium Delay',numCorrections);
readStats(accL_SIrat, accL_AErat,parametric,stat_test,'Choice Accuracy Long Delay',numCorrections);



