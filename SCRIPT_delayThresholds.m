%% Find delay thresholds for short, medium, and long delays using kmeans clustering and bootstrapping procedure
dataDirectory = '/Users/haileyrosenblum/Desktop/Matlab/Data output';
cd(dataDirectory)
load("data_analysis_aroundCP_22-Sep-2023.mat")

% remove SI rat 6 session 7 trial 35 - failed choice exit 
delayData{2}{6,7}(35)= [];
accBoolean{2}{6,7}(35)= [];

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


%return cluster data --- but clusters may not be in correct order --- assign S,M,L
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
        % note to self: 0.5 is too restrictive-- too much data being exluded
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

figure('Color', 'w'); %check for normal distribution
hist(randS_cat)
xlabel('delay length(s)')
ylabel('frequency')

    
%% Clean delay and accuracy data
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

%% Method 1: using random sampling+kmeans to determine delay min/max at short, medium, and long
%prep cells
delayShort{1}  = cell(size(delayData{1}));
delayShort{2}  = cell(size(delayData{2}));
delayMedium{1} = cell(size(delayData{1}));
delayMedium{2} = cell(size(delayData{2}));
delayLong{1}   = cell(size(delayData{1}));
delayLong{2}   = cell(size(delayData{2}));

accShort{1} = cell(size(delayData{1}));
accShort{2} = cell(size(delayData{2}));
accMed{1}   = cell(size(delayData{1}));
accMed{2}   = cell(size(delayData{2}));
accLong{1}  = cell(size(delayData{1}));
accLong{2}  = cell(size(delayData{2}));

numTrialShort{1} = cell(size(delayData{1}));
numTrialShort{2} = cell(size(delayData{2}));
numTrialMed{1}   = cell(size(delayData{1}));
numTrialMed{2}   = cell(size(delayData{2}));
numTrialLong{1}  = cell(size(delayData{1}));
numTrialLong{2}  = cell(size(delayData{2}));

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
              for minMax = 1:length(shortMin)
                  tempShort= []; tempMed= []; tempLong= [];
                  tempShort = find(shortMin(minMax) < trialDelays & trialDelays < shortMax(minMax));
                  tempMed   = find(medMin(minMax) < trialDelays & trialDelays < medMax(minMax));
                  tempLong  = find(longMin(minMax) < trialDelays & trialDelays < longMax(minMax));

                  %skip session if there aren't at least 3 trials of each delay
                  %length
                  if length(tempShort)< 3 || length(tempMed)<3 || length(tempLong)<3
                      delayShort{condi}{rati,sessi}(minMax) = NaN;
                      delayMedium{condi}{rati,sessi}(minMax) = NaN;
                      delayLong{condi}{rati,sessi}(minMax) = NaN;

                      accShort{condi}{rati,sessi}(minMax) = NaN;
                      accMed{condi}{rati,sessi}(minMax) = NaN;
                      accLong{condi}{rati,sessi}(minMax) = NaN;

                      numTrialShort{condi}{rati,sessi}(minMax) = NaN;
                      numTrialMed{condi}{rati,sessi}(minMax) = NaN;
                      numTrialLong{condi}{rati,sessi}(minMax) = NaN;
                      continue
                  end

                  %Get short medium long delays from each session
                  delayShort{condi}{rati,sessi}(minMax)  = mean(trialDelays(tempShort));
                  delayMedium{condi}{rati,sessi}(minMax) = mean(trialDelays(tempMed));
                  delayLong{condi}{rati,sessi}(minMax)   = mean(trialDelays(tempLong));

                  %Find mean accuracy for short medium long delays for each
                  %session and store
                  accShort{condi}{rati,sessi}(minMax) = (1-nanmean(tempAcc(tempShort)))*100;
                  accMed{condi}{rati,sessi}(minMax)   = (1-nanmean(tempAcc(tempMed)))*100;
                  accLong{condi}{rati,sessi}(minMax)  = (1-nanmean(tempAcc(tempLong)))*100;

                  %Get number of trials contributing to mean
                  numTrialShort{condi}{rati,sessi}(minMax) = length(tempShort);
                  numTrialMed{condi}{rati,sessi}(minMax)   = length(tempMed);
                  numTrialLong{condi}{rati,sessi}(minMax)  = length(tempLong);
              end
        end
    end
end

%SESSIONS
numTrial{1} = nan(size(delayData{1}));
numTrial{2}  = nan(size(delayData{2}));

accShort_cat = []; accMed_cat = []; accLong_cat = []; numTrialShort_cat = []; numTrialMed_cat = []; numTrialLong_cat = [];
meanTrialShort = []; meanTrialMed = []; meanTrialLong = []; numSessShort = []; numSessMed = []; numSessLong = [];
for condi = 1:length(accShort)
    accShort_cat{condi} = vertcat(accShort{condi}{:});
    accMed_cat{condi}   = vertcat(accMed{condi}{:});
    accLong_cat{condi}  = vertcat(accLong{condi}{:});


    numTrialShort_cat{condi} = vertcat(numTrialShort{condi}{:});
    numTrialMed_cat{condi}   = vertcat(numTrialMed{condi}{:});
    numTrialLong_cat{condi}  = vertcat(numTrialLong{condi}{:});

    % num trials
    % by rat over all conditions, check at +-/ 0.75 sd / idx position 2 
    for rati = 1:size(accShort{condi},1)
        for sessi=1:size(accShort{condi},2)
            if isempty(accShort{condi}{rati,sessi})
                continue
            end

            for sd=find(z==0.75)
                trialTotal = [];
                trialTotal = sum([numTrialShort{condi}{rati, sessi}(sd), numTrialMed{condi}{rati, sessi}(sd), numTrialLong{condi}{rati, sessi}(sd)]);
                if isnan(trialTotal)
                    trialTotal = 0;
                end
                numTrial{condi}(rati,sessi) = trialTotal;
            end
        end
    end

    % meanTrialShort{condi} = mean(numTrialShort_cat{condi},1,'omitnan');
    % meanTrialMed{condi}   = mean(numTrialMed_cat{condi},1,'omitnan');
    % meanTrialLong{condi}  = mean(numTrialLong_cat{condi},1,'omitnan');

    %num sessions (should be same for each condition across delays) --
    %check
    for i = 1:length(z)
        numSessShort{condi}(i) = nnz(~isnan(numTrialShort_cat{condi}(:,i)));
        numSessMed{condi}(i)   = nnz(~isnan(numTrialMed_cat{condi}(:,i)));
        numSessLong{condi}(i)  = nnz(~isnan(numTrialLong_cat{condi}(:,i)));
    end
end

%check how rat data is being represented in our calculations...
%visualize trials per session with at least 3 trials of short, medium, and
%long delays
%AE
data2plot = []; data2plot{1}=numTrial{1}(1,:); data2plot{2}=numTrial{1}(2,:); data2plot{3}=numTrial{1}(3,:); ...
    data2plot{4}=numTrial{1}(4,:); data2plot{5}=numTrial{1}(5,:); data2plot{6}=numTrial{1}(6,:);  data2plot{7}=numTrial{1}(7,:);...
    data2plot{8}=numTrial{1}(8,:); data2plot{9}=numTrial{1}(9,:); data2plot{10}=numTrial{1}(10,:);data2plot{11}=numTrial{1}(11,:); 
scat_color{1} = 'b'; scat_color{2} = 'r'; scat_color{3} = 'g'; scat_color{4} = 'm'; scat_color{5} = 'c'; scat_color{6} = 'k'; scat_color{7} = 'y'; scat_color{8} = 'w'; scat_color{9} = 'b'; scat_color{10} = 'r'; scat_color{11} = 'g'; 
verticalScatter(data2plot,scat_color,[{'rat 1'},{'rat 2'},{'rat 3'},{'rat 4'},{'rat 5'},{'rat 6'},{'rat 7'},{'rat 8'},{'rat 9'},{'rat 10'},{'rat 11'}],'Number of usable trials in session', 'mean');
title('AE +/- 0.75sd: Number of trials per session (only including sessions with at least 3 trials of short, medium. and long delays)')

data2plot = []; data2plot{1}=numTrial{2}(1,:); data2plot{2}=numTrial{2}(2,:); data2plot{3}=numTrial{2}(3,:); ...
    data2plot{4}=numTrial{2}(4,:); data2plot{5}=numTrial{2}(5,:); data2plot{6}=numTrial{2}(6,:);  data2plot{7}=numTrial{2}(7,:);
scat_color{1} = 'b'; scat_color{2} = 'r'; scat_color{3} = 'g'; scat_color{4} = 'm'; scat_color{5} = 'c'; scat_color{6} = 'k'; scat_color{7} = 'y'; 
verticalScatter(data2plot,scat_color,[{'rat 1'},{'rat 2'},{'rat 3'},{'rat 4'},{'rat 5'},{'rat 6'},{'rat 7'}],'Number of usable trials in session', 'mean');
title('SI +/- 0.75sd: Number of trials per session (only including sessions with at least 3 trials of short, medium. and long delays)')


%total number of sessions
numSessTotal = sum(vertcat(numSessShort{:}),1);

%avg number of trials per session
meanTrialShort = mean(vertcat(numTrialShort_cat{:}), 1, 'omitnan');
meanTrialMed   = mean(vertcat(numTrialMed_cat{:}),1,'omitnan');
meanTrialLong  = mean(vertcat(numTrialLong_cat{:}),1,'omitnan');


%ttest and avg num trials per session
tAccShort = []; tAccMed = []; tAccLong = [];
for i = 1:length(z)
    %ttest2
    [~, tAccShort(i), ~, tstatShort(i)] = ttest2(accShort_cat{1}(:,i), accShort_cat{2}(:,i));
    [~, tAccMed(i), ~, tstatMed(i)]   = ttest2(accMed_cat{1}(:,i), accMed_cat{2}(:,i));
    [~, tAccLong(i), ~, tstatLong(i)]  = ttest2(accLong_cat{1}(:,i), accLong_cat{2}(:,i));
end

%make x values delay ranges
labelS = {}; labelM = {}; labelL = {};
for i = 1: length(shortMin)
    tempLabel = [];
    tempLabel = {num2str(shortMin(i)) num2str('-') num2str(shortMax(i))};
    labelS(1,i) = join(tempLabel);
end
for i = 1: length(medMin)
    tempLabe  = [];
    tempLabel = {num2str(medMin(i)) num2str('-') num2str(medMax(i))};
    labelM(i) = join(tempLabel);
end
for i = 1: length(longMin)
    tempLabel = [];
    tempLabel = {num2str(longMin(i)) num2str('-') num2str(longMax(i))};
    labelL(i) = join(tempLabel);
end

%medium delay 
figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccMed);
xtickangle(45)
xlabel('Delay range (seconds)');
ylabel('choice accuracy p value')
title('Medium delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,meanTrialMed)
ylim([0, 12]);
ylabel('avg number of trials per session');
xticks(z)
xticklabels(labelM)

%could also use tstat instead of p

figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccMed);
ylabel('choice accuracy p value')
xlabel('Delay range (seconds)');
title('Medium delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,numSessTotal)
ylabel('total sessions');
xticks(z)
xticklabels(labelM)

%Long delay 
figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccLong);
xtickangle(45)
ylabel('choice accuracy p value')
xlabel('Delay range (seconds)');
title('Long delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,meanTrialLong)
ylim([0, 12]);
ylabel('avg number of trials per session');
xticks(z)
xticklabels(labelL)


figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccLong);
ylabel('choice accuracy p value')
xlabel('Delay range (seconds)');
title('Long delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,numSessTotal)
ylabel('total sessions');
xticks(z)
xticklabels(labelL)

%Short delay 
figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccShort);
xtickangle(45)
ylabel('choice accuracy p value')
xlabel('Delay range (seconds)');
title('Short delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,meanTrialShort)
ylim([0, 12]);
ylabel('avg number of trials per session');
xticks(z)
xticklabels(labelS)

figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccShort);
ylabel('choice accuracy p value')
xlabel('Delay range (seconds)');
title('Short delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,numSessTotal)
ylabel('total sessions');
xticks(z)
xticklabels(labelS)


% see if what we are defining as short, medium, and long is the same for AE
% and SI sessions
shortMeans = []; medMeans = []; longMeans = []; delayS_th = []; delayM_th = []; delayL_th = [];
for condi = 1:length(delayShort)
    for rati = 1:size(delayShort{condi},1)
        for sessi= 1:size(delayShort{condi},2)
            %skip nonexistsent sessions
            if isempty(delayShort{condi}{rati,sessi})
                continue
            end

            %get mean short, medium, and long delay for each sd threshold
            for threshold = 1:length(delayShort{condi}{rati,sessi})
                shortMeans{condi}{rati}(sessi,threshold) = delayShort{condi}{rati,sessi}(threshold);
                medMeans{condi}{rati}(sessi,threshold) = delayMedium{condi}{rati,sessi}(threshold);
                longMeans{condi}{rati}(sessi,threshold) = delayLong{condi}{rati, sessi}(threshold);
            end
        end
    end
    delayS_th{condi} = vertcat(shortMeans{condi}{:});
    delayM_th{condi} = vertcat(medMeans{condi}{:});
    delayL_th{condi} = vertcat(longMeans{condi}{:});
end

shortSig = []; medSig = []; longSig = []; idx_sd = [];
%lots of ttest2
for threshold = 1:size(delayS_th{condi},2)
    shortSig(threshold) = ttest2(delayS_th{1}(:,threshold), delayS_th{2}(:,threshold));
    medSig(threshold) = ttest2(delayM_th{1}(:,threshold), delayM_th{2}(:,threshold));
    longSig(threshold) = ttest2(delayS_th{1}(:,threshold), delayS_th{2}(:,threshold));

    if shortSig(1,threshold)==0 && medSig(1,threshold) ==0 && longSig(1,threshold)==0
        idx_sd(end+1) = threshold;
        disp(['Not significant: +/- ' num2str(z(threshold)), ' standard deviations'])
    end
end

% let's check
% looping
rows = 2; cols = 5;
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
figure('Color', 'w');
for i = 1:length(idx_sd)
    subplot(rows, cols, i); hold on;
    set(gcf,'Position', (get(0, 'Screensize')));

    meanDelay_AE = []; meanDelay_SI = [];
    meanDelay_AE = [mean(delayS_th{1}(:,i),'omitnan'), mean(delayM_th{1}(:,i),'omitnan'), mean(delayL_th{1}(:,i),'omitnan')];
    meanDelay_SI = [mean(delayS_th{2}(:,i),'omitnan'), mean(delayM_th{2}(:,i),'omitnan'), mean(delayL_th{2}(:,i),'omitnan')];

    errDelayAE = [stderr(delayS_th{1}(:,i),1), stderr(delayM_th{1}(:,i),1), stderr(delayL_th{1}(:,i),1)];
    errDelaySI = [stderr(delayS_th{2}(:,i),1), stderr(delayM_th{2}(:,i),1), stderr(delayL_th{2}(:,i),1)];

    errorbar(group, meanDelay_AE, errDelayAE, 'r', 'LineWidth', 2);
    errorbar(group, meanDelay_SI, errDelaySI, 'b', 'LineWidth', 2);
    ylabel('Delay length (s)')
    title(['+/-' num2str(z(idx_sd(i))) 'sd'])
end



%% Method 1: Using predetermined threshold to choice accuracy

sd = find(z==0.75); %get index of chosen threshold
accShort_cat{1}(:,sd);
accMed_cat{1}(:,sd);
accLong_cat{1}(:,sd);

% Accuracy Error Bar plot
meanAccAE = [mean(accShort_cat{1}(:,sd),'omitnan'), mean(accMed_cat{1}(:,sd),'omitnan'), mean(accLong_cat{1}(:,sd),'omitnan')];
meanAccSI = [mean(accShort_cat{2}(:,sd),'omitnan'), mean(accMed_cat{2}(:,sd),'omitnan'), mean(accLong_cat{2}(:,sd),'omitnan')];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(accShort_cat{1}(:,sd),1), stderr(accMed_cat{1}(:,sd),1), stderr(accLong_cat{1}(:,sd),1)];
errSI = [stderr(accShort_cat{2}(:,sd),1), stderr(accMed_cat{2}(:,sd),1), stderr(accLong_cat{2}(:,sd),1)];
figure('Color', 'w'); hold on;
errorbar(group, meanAccAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanAccSI,errSI, 'b', 'LineWidth', 2);
title('N=sessions')
ylabel('Choice Accuracy')
xlabel('Delays set by data')

stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(accShort_cat{1}(:,sd),accShort_cat{2}(:,sd),parametric,stat_test,'Choice Accuracy Short Delay',numCorrections);
readStats(accMed_cat{1}(:,sd),accMed_cat{2}(:,sd),parametric,stat_test,'Choice Accuracy Medium Delay',numCorrections);
readStats(accLong_cat{1}(:,sd),accLong_cat{2}(:,sd),parametric,stat_test,'Choice Accuracy Long Delay',numCorrections);


%% Method 2: What if we extend delay range a bit? min= defined range, max=value defined by kmeans+random sampling
delayShort{1}  = cell(size(delayData{1}));
delayShort{2}  = cell(size(delayData{2}));
delayMedium{1} = cell(size(delayData{1}));
delayMedium{2} = cell(size(delayData{2}));
delayLong{1}   = cell(size(delayData{1}));
delayLong{2}   = cell(size(delayData{2}));

accShort{1} = cell(size(delayData{1}));
accShort{2} = cell(size(delayData{2}));
accMed{1}   = cell(size(delayData{1}));
accMed{2}   = cell(size(delayData{2}));
accLong{1}  = cell(size(delayData{1}));
accLong{2}  = cell(size(delayData{2}));

numTrialShort{1} = cell(size(delayData{1}));
numTrialShort{2} = cell(size(delayData{2}));
numTrialMed{1}   = cell(size(delayData{1}));
numTrialMed{2}   = cell(size(delayData{2}));
numTrialLong{1}  = cell(size(delayData{1}));
numTrialLong{2}  = cell(size(delayData{2}));

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
              for minMax = 1:length(shortMin)
                  tempShort= []; tempMed= []; tempLong= [];
                  tempShort = find(10 < trialDelays & trialDelays < shortMax(minMax));
                  tempMed   = find(30 < trialDelays & trialDelays < medMax(minMax));
                  tempLong  = find(60 < trialDelays & trialDelays < longMax(minMax));

                  %skip session if there aren't at least 3 trials of each delay
                  %length
                  if length(tempShort)< 3 || length(tempMed)<3 || length(tempLong)<3
                      delayShort{condi}{rati,sessi}(minMax) = NaN;
                      delayMedium{condi}{rati,sessi}(minMax) = NaN;
                      delayLong{condi}{rati,sessi}(minMax) = NaN;

                      accShort{condi}{rati,sessi}(minMax) = NaN;
                      accMed{condi}{rati,sessi}(minMax) = NaN;
                      accLong{condi}{rati,sessi}(minMax) = NaN;

                      numTrialShort{condi}{rati,sessi}(minMax) = NaN;
                      numTrialMed{condi}{rati,sessi}(minMax) = NaN;
                      numTrialLong{condi}{rati,sessi}(minMax) = NaN;
                      continue
                  end

                  %Get short medium long delays from each session
                  delayShort{condi}{rati,sessi}(minMax)  = mean(trialDelays(tempShort));
                  delayMedium{condi}{rati,sessi}(minMax) = mean(trialDelays(tempMed));
                  delayLong{condi}{rati,sessi}(minMax)   = mean(trialDelays(tempLong));

                  %Find mean accuracy for short medium long delays for each
                  %session and store
                  accShort{condi}{rati,sessi}(minMax) = (1-nanmean(tempAcc(tempShort)))*100;
                  accMed{condi}{rati,sessi}(minMax)   = (1-nanmean(tempAcc(tempMed)))*100;
                  accLong{condi}{rati,sessi}(minMax)  = (1-nanmean(tempAcc(tempLong)))*100;

                  %Get number of trials contributing to mean
                  numTrialShort{condi}{rati,sessi}(minMax) = length(tempShort);
                  numTrialMed{condi}{rati,sessi}(minMax)   = length(tempMed);
                  numTrialLong{condi}{rati,sessi}(minMax)  = length(tempLong);
              end
        end
    end
end

% see if what we are defining as short, medium, and long is the same for AE
% and SI sessions
numTrial{1} = nan(size(delayData{1}));
numTrial{2}  = nan(size(delayData{2}));

shortMeans = []; medMeans = []; longMeans = []; delayS_th = []; delayM_th = []; delayL_th = [];
for condi = 1:length(delayShort)
    for rati = 1:size(delayShort{condi},1)
        for sessi= 1:size(delayShort{condi},2)
            %skip nonexistsent sessions
            if isempty(delayShort{condi}{rati,sessi})
                continue
            end

            %get mean short, medium, and long delay for each sd threshold
            for threshold = 1:length(delayShort{condi}{rati,sessi})
                shortMeans{condi}{rati}(sessi,threshold) = delayShort{condi}{rati,sessi}(threshold);
                medMeans{condi}{rati}(sessi,threshold) = delayMedium{condi}{rati,sessi}(threshold);
                longMeans{condi}{rati}(sessi,threshold) = delayLong{condi}{rati, sessi}(threshold);
            end
        end
    end
    delayS_th{condi} = vertcat(shortMeans{condi}{:});
    delayM_th{condi} = vertcat(medMeans{condi}{:});
    delayL_th{condi} = vertcat(longMeans{condi}{:});
end

shortSig = []; medSig = []; longSig = []; idx_sd = [];
%lots of ttest2
for threshold = 1:size(delayS_th{condi},2)
    shortSig(threshold) = ttest2(delayS_th{1}(:,threshold), delayS_th{2}(:,threshold));
    medSig(threshold) = ttest2(delayM_th{1}(:,threshold), delayM_th{2}(:,threshold));
    longSig(threshold) = ttest2(delayS_th{1}(:,threshold), delayS_th{2}(:,threshold));

    if shortSig(1,threshold)==0 && medSig(1,threshold) ==0 && longSig(1,threshold)==0
        idx_sd(end+1) = threshold;
        disp(['Not significant: +/- ' num2str(z(threshold)), ' standard deviations'])
    end
end

% let's check
%looping
rows = 2; cols = 5;
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
figure('Color', 'w');
for i = 1:length(idx_sd)
    subplot(rows, cols, i); hold on;
    set(gcf,'Position', (get(0, 'Screensize')));

    meanDelay_AE = []; meanDelay_SI = [];
    meanDelay_AE = [mean(delayS_th{1}(:,i),'omitnan'), mean(delayM_th{1}(:,i),'omitnan'), mean(delayL_th{1}(:,i),'omitnan')];
    meanDelay_SI = [mean(delayS_th{2}(:,i),'omitnan'), mean(delayM_th{2}(:,i),'omitnan'), mean(delayL_th{2}(:,i),'omitnan')];

    errDelayAE = [stderr(delayS_th{1}(:,i),1), stderr(delayM_th{1}(:,i),1), stderr(delayL_th{1}(:,i),1)];
    errDelaySI = [stderr(delayS_th{2}(:,i),1), stderr(delayM_th{2}(:,i),1), stderr(delayL_th{2}(:,i),1)];

    errorbar(group, meanDelay_AE, errDelayAE, 'r', 'LineWidth', 2);
    errorbar(group, meanDelay_SI, errDelaySI, 'b', 'LineWidth', 2);
    ylabel('Delay length (s)')
    title(['+/-' num2str(z(idx_sd(i))) 'sd'])
end

%SESSIONS
accShort_cat = []; accMed_cat = []; accLong_cat = []; numTrialShort_cat = []; numTrialMed_cat = []; numTrialLong_cat = [];
meanTrialShort = []; meanTrialMed = []; meanTrialLong = []; numSessShort = []; numSessMed = []; numSessLong = [];
for condi = 1:length(accShort)
    accShort_cat{condi} = vertcat(accShort{condi}{:});
    accMed_cat{condi}   = vertcat(accMed{condi}{:});
    accLong_cat{condi}  = vertcat(accLong{condi}{:});


    numTrialShort_cat{condi} = vertcat(numTrialShort{condi}{:});
    numTrialMed_cat{condi}   = vertcat(numTrialMed{condi}{:});
    numTrialLong_cat{condi}  = vertcat(numTrialLong{condi}{:});

    %num sessions (should be same for each condition across delays) --
    %check
    for i = 1:length(z)
        numSessShort{condi}(i) = nnz(~isnan(numTrialShort_cat{condi}(:,i)));
        numSessMed{condi}(i)   = nnz(~isnan(numTrialMed_cat{condi}(:,i)));
        numSessLong{condi}(i)  = nnz(~isnan(numTrialLong_cat{condi}(:,i)));
    end

    % %num trials
    %by rat over all conditions, check at +-/ 0.75 sd / idx position 2
    for rati = 1:size(accShort{condi},1)
        for sessi=1:size(accShort{condi},2)
            if isempty(accShort{condi}{rati,sessi})
                continue
            end

            for sd=find(z==2.25)
                trialTotal = [];
                trialTotal = sum([numTrialShort{condi}{rati, sessi}(sd), numTrialMed{condi}{rati, sessi}(sd), numTrialLong{condi}{rati, sessi}(sd)]);
                if isnan(trialTotal)
                    trialTotal = 0;
                end
                numTrial{condi}(rati,sessi) = trialTotal;
            end
        end
    end

end

%check how rat data is being represented in our calculations...
%visualize trials per session with at least 3 trials of short, medium, and
%long delays
%AE
data2plot = []; data2plot{1}=numTrial{1}(1,:); data2plot{2}=numTrial{1}(2,:); data2plot{3}=numTrial{1}(3,:); ...
    data2plot{4}=numTrial{1}(4,:); data2plot{5}=numTrial{1}(5,:); data2plot{6}=numTrial{1}(6,:);  data2plot{7}=numTrial{1}(7,:);...
    data2plot{8}=numTrial{1}(8,:); data2plot{9}=numTrial{1}(9,:); data2plot{10}=numTrial{1}(10,:);data2plot{11}=numTrial{1}(11,:); 
figure;
scat_color{1} = 'b'; scat_color{2} = 'r'; scat_color{3} = 'g'; scat_color{4} = 'm'; scat_color{5} = 'c'; scat_color{6} = 'k'; scat_color{7} = 'y'; scat_color{8} = 'w'; scat_color{9} = 'b'; scat_color{10} = 'r'; scat_color{11} = 'g'; 
verticalScatter(data2plot,scat_color,[{'rat 1'},{'rat 2'},{'rat 3'},{'rat 4'},{'rat 5'},{'rat 6'},{'rat 7'},{'rat 8'},{'rat 9'},{'rat 10'},{'rat 11'}],'Number of usable trials in session', 'mean');
title('AE +/- 2.25sd: Number of trials per session (only including sessions with at least 3 trials of short, medium. and long delays)')

data2plot = []; data2plot{1}=numTrial{2}(1,:); data2plot{2}=numTrial{2}(2,:); data2plot{3}=numTrial{2}(3,:); ...
    data2plot{4}=numTrial{2}(4,:); data2plot{5}=numTrial{2}(5,:); data2plot{6}=numTrial{2}(6,:);  data2plot{7}=numTrial{2}(7,:);
figure;
scat_color{1} = 'b'; scat_color{2} = 'r'; scat_color{3} = 'g'; scat_color{4} = 'm'; scat_color{5} = 'c'; scat_color{6} = 'k'; scat_color{7} = 'y'; 
verticalScatter(data2plot,scat_color,[{'rat 1'},{'rat 2'},{'rat 3'},{'rat 4'},{'rat 5'},{'rat 6'},{'rat 7'}],'Number of usable trials in session', 'mean');
title('SI +/- 2.25sd: Number of trials per session (only including sessions with at least 3 trials of short, medium. and long delays)')


%total number of sessions
numSessTotal = sum(vertcat(numSessShort{:}),1);

%avg number of trials per session
meanTrialShort = mean(vertcat(numTrialShort_cat{:}), 1, 'omitnan');
meanTrialMed   = mean(vertcat(numTrialMed_cat{:}),1,'omitnan');
meanTrialLong  = mean(vertcat(numTrialLong_cat{:}),1,'omitnan');


%ttest and avg num trials per session
tAccShort = []; tAccMed = []; tAccLong = [];
for i = 1:length(z)
    %ttest2
    [~, tAccShort(i), ~, tstatShort(i)] = ttest2(accShort_cat{1}(:,i), accShort_cat{2}(:,i));
    [~, tAccMed(i), ~, tstatMed(i)]   = ttest2(accMed_cat{1}(:,i), accMed_cat{2}(:,i));
    [~, tAccLong(i), ~, tstatLong(i)]  = ttest2(accLong_cat{1}(:,i), accLong_cat{2}(:,i));
end

%make x values delay ranges
labelS = {}; labelM = {}; labelL = {};
for i = 1: length(shortMin)
    tempLabel = [];
    tempLabel = {'10' num2str('-') num2str(shortMax(i))};
    labelS(1,i) = join(tempLabel);
end
for i = 1: length(medMin)
    tempLabe  = [];
    tempLabel = {'30' num2str('-') num2str(medMax(i))};
    labelM(i) = join(tempLabel);
end
for i = 1: length(longMin)
    tempLabel = [];
    tempLabel = {'60' num2str('-') num2str(longMax(i))};
    labelL(i) = join(tempLabel);
end

%medium delay 
figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccMed);
xtickangle(45)
xlabel('Delay range (seconds)');
ylabel('p value')
title('Medium delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,meanTrialMed)
ylim([0, 12]);
ylabel('avg number of trials per session');
xticks(z)
xticklabels(labelM)

figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccMed);
ylabel('choice accuracy p value')
xlabel('Delay range (seconds)');
title('Medium delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,numSessTotal)
ylabel('total sessions');
xticks(z)
xticklabels(labelM)

%Long delay 
figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccLong);
xtickangle(45)
ylabel('p value')
xlabel('Delay range (seconds)');
title('Long delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,meanTrialLong)
ylim([0, 12]);
ylabel('avg number of trials per session');
xticks(z)
xticklabels(labelL)

%Short delay 
figure('Color', 'w'); hold on;
yyaxis left
plot(z, tAccShort);
xtickangle(45)
ylabel('p value')
xlabel('Delay range (seconds)');
title('Short delay: significance, delay thresholds, and amount of data included')
yyaxis right
plot(z,meanTrialShort)
ylim([0, 12]);
ylabel('avg number of trials per session');
xticks(z)
xticklabels(labelS)

%% Method 2: Using predetermined threshold to examine choice accuracy

sd = find(z==2.25); %get index of chosen threshold
accShort_cat{1}(:,sd);
accMed_cat{1}(:,sd);
accLong_cat{1}(:,sd);

% Accuracy Error Bar plot
meanAccAE = [mean(accShort_cat{1}(:,sd),'omitnan'), mean(accMed_cat{1}(:,sd),'omitnan'), mean(accLong_cat{1}(:,sd),'omitnan')];
meanAccSI = [mean(accShort_cat{2}(:,sd),'omitnan'), mean(accMed_cat{2}(:,sd),'omitnan'), mean(accLong_cat{2}(:,sd),'omitnan')];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(accShort_cat{1}(:,sd),1), stderr(accMed_cat{1}(:,sd),1), stderr(accLong_cat{1}(:,sd),1)];
errSI = [stderr(accShort_cat{2}(:,sd),1), stderr(accMed_cat{2}(:,sd),1), stderr(accLong_cat{2}(:,sd),1)];
figure('Color', 'w'); hold on;
errorbar(group, meanAccAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanAccSI,errSI, 'b', 'LineWidth', 2);
title('N=sessions')
ylabel('Choice Accuracy')
xlabel('Delay min set by experiment, max set by data')

stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(accShort_cat{1}(:,sd),accShort_cat{2}(:,sd),parametric,stat_test,'Choice Accuracy Short Delay',numCorrections);
readStats(accMed_cat{1}(:,sd),accMed_cat{2}(:,sd),parametric,stat_test,'Choice Accuracy Medium Delay',numCorrections);
readStats(accLong_cat{1}(:,sd),accLong_cat{2}(:,sd),parametric,stat_test,'Choice Accuracy Long Delay',numCorrections);


%% for testing using manually set delays 
allShort{1} = nan(size(delayData{1}));
allShort{2} = nan(size(delayData{2}));
allMed{1}   = nan(size(delayData{1}));
allMed{2}   = nan(size(delayData{2}));
allLong{1}  = nan(size(delayData{1}));
allLong{2}  = nan(size(delayData{2}));

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
            tempShort = find(10 < trialDelays & trialDelays < 24);
            tempMed   = find(30 < trialDelays & trialDelays < 45);
            tempLong  = find(60 < trialDelays & trialDelays < 79);

            %skip session if there aren't at least 3 trials of eacf delay
            %length
            if length(tempShort)< 3 || length(tempMed)<3 || length(tempLong)<3
                continue
            end

            %Get short medium long delays from each session
            allShort{condi}(rati,sessi)  = mean(trialDelays(tempShort));
            allMed{condi}(rati,sessi) = mean(trialDelays(tempMed));
            allLong{condi}(rati,sessi)   = mean(trialDelays(tempLong));

            %Find mean accuracy for short medium long delays for each
            %session
            accShort{condi}(rati,sessi) = (1-nanmean(tempAcc(tempShort)))*100;
            accMed{condi}(rati,sessi)   = (1-nanmean(tempAcc(tempMed)))*100;
            accLong{condi}(rati,sessi)  = (1-nanmean(tempAcc(tempLong)))*100;
        end
    end
end


% Delay Error Bar plot
meanDelayAE = [mean(allShort{1},'all', 'omitnan'), mean(allMed{1},'all','omitnan'), mean(allLong{1},'all','omitnan')];
meanDelaySI = [mean(allShort{2},'all','omitnan'), mean(allMed{2},'all','omitnan'), mean(allLong{2},'all','omitnan')];
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = [stderr(vertcat(allShort{1}(:)),1), stderr(vertcat(allMed{1}(:)),1), stderr(vertcat(allMed{1}(:)),1)];
errSI = [stderr(vertcat(allShort{2}(:)),1), stderr(vertcat(allMed{1}(:)),1), stderr(vertcat(allMed{1}(:)),1)];
figure('Color', 'w'); hold on;
errorbar(group, meanDelayAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, meanDelaySI,errSI, 'b', 'LineWidth', 2);

[hld,pld] = ttest2(vertcat(allShort{1}(:)), vertcat(allShort{2}(:)))
[hmd,pmd] = ttest2(vertcat(allMed{1}(:)), vertcat(allMed{2}(:)))
[hsd,psd] = ttest2(vertcat(allLong{1}(:)), vertcat(allLong{2}(:)))

