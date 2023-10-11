%% VTE analysis
clear;
datafolder = '/Users/haileyrosenblum/Desktop/Matlab/Data output';
cd(datafolder);
load('data_analysis_stem2cpExit_22-Sep-2023.mat');

%first remove the double trial empty cell from AE R1 S20 -- row 17 (trial 18 in
%sequence cell) and AE R8 S5 (trial 11 in sequence cell) or else index won't match w/ accuracy and delay cells
posData{1}{1,20}(17)=[];
%posData{1}{8,14}(10)=[]; outlier rat, note to self: forgot to
%remove this trial before going through VTE, so there will be an extra NaN
%in the cell when everything is removed -- update fixed when saved sortVTE

% remove SI rat 6 session 7 trial 35 - failed choice exit 
posData{2}{6,7}(35)= [];

% empty trials --> NaN
% NOTE TO SELF: More delay trials are removed when those cells are cleaned
% (failed SB entry and failed stem entry). Address by keepping VTE data structure same as data
% files before cleaning empty cells, index later 
posAESI{1} = cell(size(posData{1}));
posAESI{2} = cell(size(posData{2}));

for condi= 1:length(posData)
    for rati=1:size(posData{condi},1)
        for sessi=1:size(posData{condi},2)
            if isempty(posData{condi}{rati,sessi})
                continue
            end
            posAESI{condi}{rati,sessi} = empty2nan(posData{condi}{rati,sessi});
        end
    end
end

%% Get zlnIdPhi and shuffle AE/SI data
window_sec = 1; %window of time to consider for adaptive estimation - 1 sec works
postSmoothing = 0.5; %smoothing window time (0.5 sec was from Redish)
vt_srate = 30; %video tracking sampling 
display = 0; %set to 0 for no display

lnIdPhi = []; posSmooth=[];
for condi = 1:length(posAESI)
    for rati = 1:size(posAESI{condi},1)
        for sessi= 1:length(posAESI{condi}(rati,:))
            if isempty(posAESI{condi}{rati,sessi})
                continue
            end

            for triali= 1:length(posAESI{condi}{rati,sessi})
                if isnan(posAESI{condi}{rati,sessi}{triali})
                    lnIdPhi{condi}{rati,sessi}(triali)= NaN;
                    posSmooth{condi}{rati,sessi}{triali} = NaN;
                    continue
                end

                posX = []; posY = [];
                posX = posAESI{condi}{rati,sessi}{triali}(1,:);
                posY = posAESI{condi}{rati,sessi}{triali}(2,:);

                %get x and y from stem to cp exit
                TempX = []; TempY = []; idxMid = [];
                idxMid = posX <= 450;
                TempX = posX(idxMid); % x position data as a vector
                TempY = posY(idxMid); % y position data as a vector

                %time spent
                tsTemp = []; ts = [];
                tsTemp = posAESI{condi}{rati,sessi}{triali}(3,:);
                ts = tsTemp(idxMid);

                % smooth position data
                X = []; Y = [];
                X = smoothdata(TempX, 'gaussian', 30);
                Y = smoothdata(TempY, 'gaussian', 30);
                posSmooth{condi}{rati,sessi}{triali} = [X;Y;ts];


                %get IdPhi for each trial, then ln of IdPhi
                IdPhi = [];
                IdPhi  = IdPhi_RedishFun(X,Y,window_sec,postSmoothing,vt_srate,display);
                lnIdPhi{condi}{rati,sessi}(triali) = log(IdPhi);
            end
        end
       
    end
end

%organize idphi, position, condition, rat, session, trial information --
%preparing to look at VTE trials
zlnIdPhi=[];
for condi = 1:length(lnIdPhi)
    temp_zlnIdPhi = [];
    for rati = 1:size(lnIdPhi{condi},1)
        array={};
        for sessi = 1:length(lnIdPhi{condi}(rati,:))
            if isempty(lnIdPhi{condi}{rati,sessi})
                continue
            end

            temp = {};
            temp(1,:) = num2cell(lnIdPhi{condi}{rati,sessi}); % get idphi
            temp(2,:) = posSmooth{condi}{rati,sessi}; % get position
            temp(3,:) = num2cell(repmat(condi,[1 size(temp,2)]));  % get condition
            temp(4,:) = num2cell(repmat(rati,[1 size(temp,2)]));  % rat information
            temp(5,:) = num2cell(repmat(sessi,[1 size(temp,2)])); % get session indicator (rowi)
            temp(6,:) = num2cell(1:size(temp,2)); % trial information

            % concatenate across loops, get each rat's trial info
            array = horzcat(array,temp);
        end

        %zscore lnIdPhi by rat
        znan = [];
        znan = zscore_NaN(cell2mat(array(1,:)));
        
        %new cell with updated z scored lnIdPhi and condi,rat,sess,trial
        %info
        temp_zlnIdPhi{condi}{rati} = vertcat(num2cell(znan), array(2:end,:));
    end
    % across rats
    zlnIdPhi{condi} = horzcat(temp_zlnIdPhi{condi}{:});
end

%Blind to condition
%shuffle data
zlnIdPhi_all = horzcat(zlnIdPhi{:});
rng('default')
shuffleIdPhi = zlnIdPhi_all(:,randperm(length(zlnIdPhi_all)));

%% What is VTE threshold?
%THRESHOLD METHOD 1: Stout 2021 first clear deflection of the normal distribution
% try fitting a complex 
fig = figure('color','w');
h1 = histogram(cell2mat(shuffleIdPhi(1,:))); hold on;
h1.FaceColor = [0.8500 0.3250 0.0980];
ylimits = ylim;
title('Distribution of z[ln(IdPhi)] estimates')
box off
ylabel('Number of Trials')
xlabel('z[ln(IdPhi)]')
h1.BinWidth = 0.1;
zIdPhi_bins = h1.BinEdges;
%axis tight;
%xlim([-2 4])
pd = fitdist(cell2mat(shuffleIdPhi(1,:))','Kernel','Kernel','normal');
stepVal = 1/length(shuffleIdPhi(1,:));
x_curve = h1.BinEdges(1):stepVal:h1.BinEdges(end);
y_curve = pdf(pd,x_curve);
yyaxis right; hold on;
plot(x_curve,y_curve,'k','LineWidth',2); 
xline(0.3, '--b', "LineWidth", 2)
%VTE = zlnIdPhi = ~0.3

%THRESHOLD METHOD 2 2: Papale 2016 
%take mean zIdphi of 5 most populous bins, then take all the samples below that mean, 
%reflect them around the mean, and used that pseudo-sample to define an expected Gaussian 
%distribution of non-VTE laps
edges = -4:0.05:10;
% ACROSS GROUPS
zscoredIdPhi = cell2mat(shuffleIdPhi(1,:));
bincounts = []; idxBin = [];
[bincounts, ~, idxBin] = histcounts(zscoredIdPhi, edges);

%get 5 most populous bins
top5 = []; idxTop5 = [];
top5 = find(ismember(bincounts, maxk(bincounts, 5))==1);
idxTop5 = ismember(idxBin,top5);

%get mean zlnIdPhi of most populous bins
meanbins = mean(zscoredIdPhi(idxTop5));

%find samples below the mean, and then reflect them around the mean
belowMean = []; belowMeanreflect = [];
belowMean = sort(zscoredIdPhi(zscoredIdPhi < meanbins));
belowMeanreflect = sort(-belowMean,'ascend')+2*meanbins;

%expected gaussian and plot
expDist = []; pd = []; y = [];
expDist = [belowMean, meanbins, belowMeanreflect];
pd = fitdist(expDist','Normal');
y = pdf(pd,expDist);
figure('Color', 'w');
plot(expDist, y)

%plot actual distribution + expected distribution
figure('color','w'); hold on
histogram(zscoredIdPhi, 'BinWidth', 0.05,'FaceColor',[0.8500 0.3250 0.0980], 'Normalization','pdf'); 
plot(expDist, y, 'LineWidth', 2, 'Color', 'k')
xlim([-1.5 4])
xline(0.025, '--b', "LineWidth", 4)
title('z(ln(IdPhi)) Distribution Collapsed Across Groups')



%% Check VTE trials--Remove trials that are incorrectly classified as VTE
% METHOD 1
%plot multiple trials at a time
%sample session to plot as base
datafolder = '/Users/haileyrosenblum/Desktop/Matlab/getVT/DA10_20-78';
cd(datafolder); % change directory so  matlab can see data
missing_data = 'interp';
[x_sample,y_sample] = getVTdata(datafolder,missing_data,'VT1.mat');

% set threshold, index z(ln(IdPhi) and position/time data above threshold
vte = 0.3;
idxVTE = find(cell2mat(shuffleIdPhi(1,:)) >= vte);

% plot 
rowVal = 2;
colVal = 3;
figure('color','w'); counter = 0; idx2rem_step1 = [];
set(gcf,'Position', (get(0, 'ScreenSize')));
for i = 1:length(idxVTE)

    counter = counter+1;
    subplot(rowVal,colVal,counter); hold on;

    % this data is plotted against a general session
    plot(x_sample,y_sample,'Color',[.8 .8 .8]); box off;
    hold on;    

    % plot position with color heat indicating angular velocity
    head_velocity = [];
    [~,head_velocity] = kinematics2D(shuffleIdPhi{2,idxVTE(i)}(1,:),shuffleIdPhi{2,idxVTE(i)}(2,:),shuffleIdPhi{2,idxVTE(i)}(3,:),'y'); % 'y' to convert to seconds
    head_velocity = normalize(head_velocity,'range');
    s = scatter(shuffleIdPhi{2,idxVTE(i)}(1,1:end-1),shuffleIdPhi{2,idxVTE(i)}(2,1:end-1),[],head_velocity(:));
    s.Marker = '.';
    s.SizeData = 75;   
    c= colorbar;
    c.Label.String = 'normalized velocity';

    title(['Trial ', num2str(idxVTE(i)), ' zlnIdPhi = ',num2str(shuffleIdPhi{1,idxVTE(i)})])

    if counter == 6 || i == length(idxVTE)
        idx2rem_step1 = horzcat(idx2rem_step1, str2num(input('Enter trials falsely identified as VTE: ','s')));
        %pause;
        close;
        set(gcf,'Position', (get(0,'ScreenSize')));
        counter = 0;
    end
end

%save step 1 remove data
save('idx2rem_step1', 'idx2rem_step1') 

% STEP 2
%check rejected VTEs, if they fail again, then index to remove
idx2rem_step2 = [];
for i = 1:length(idx2rem_step1)
    figure('Color', 'w')
    set(gcf,'Position', (get(0, 'Screensize')));
    plot(x_sample,y_sample,'Color',[.8 .8 .8]); hold on;
    title(['Trial ', num2str(idx2rem_step1(i)), ' zlnIdPhi = ',num2str(shuffleIdPhi{1,idx2rem_step1(i)})])

    [~,head_velocity] = kinematics2D(shuffleIdPhi{2,idx2rem_step1(i)}(1,:),shuffleIdPhi{2,idx2rem_step1(i)}(2,:),shuffleIdPhi{2,idx2rem_step1(i)}(3,:),'y'); % 'y' to convert to seconds

    for k = 1:length(shuffleIdPhi{2,idx2rem_step1(i)}(1,:))-1
        head_velocity = normalize(head_velocity,'range');
        scatter(shuffleIdPhi{2,idx2rem_step1(i)}(1,k),shuffleIdPhi{2,idx2rem_step1(i)}(2,k),[],head_velocity(k));
        drawnow
        if k == length(shuffleIdPhi{2,idx2rem_step1(i)}(1,:))-1
            idx2rem_step2 = horzcat(idx2rem_step2, str2num(input('Enter trials falsely identified as VTE: ','s')));
            close;
        end
    end
end


%make variable to save step 2 remove data
save('idx2rem_step2', 'idx2rem_step2')

% METHOD 2 
% Checking trials w/ < threshold IdPhi, but where the rat enters both L and R goal arms
% hold variables
%look under threshold
vte = 0.3;

%identify parameters for left choice point
minY = 195; addY = abs(minY-242);
minX = 165; addX = abs(minX-240);
lCP_fld = [minX minY addX addY]; %[180 215 395 50]; % x,y (first corner) , x,y (second right top corner)

% same for right choice point
minY = 249; addY = abs(minY-295);
minX = 165; addX = abs(minX-240);
rCP_fld = [minX minY addX addY]; %[180 215 395 50]; % x,y (first corner) , x,y (second right top corner)

% now for each trial identify if the rat was in one bin before entering the
% other
cd('/Users/haileyrosenblum/Desktop/Matlab/Functions')
trialnum = [];
for triali = 1:length(shuffleIdPhi)
    if shuffleIdPhi{1,triali} <= vte %if zlnIdPhi is less than threshold value
        %is rat in left and right goal arm?
        [oopsTrial] = [];
        [oopsTrial] = hget_oopsTrials(shuffleIdPhi{2,triali}(1,:),shuffleIdPhi{2,triali}(2,:),lCP_fld,rCP_fld); 
        if oopsTrial == 1
            trialnum(end+1) = triali; %get trial index for when rat enters both goal arms
        end
    end
end

%check new VTE
rowVal = 2;
colVal = 3;
figure('Color', 'w'); hold on;
set(gcf,'Position', (get(0, 'Screensize')));
addVTE = []; counter = 0;
for triali = 90:length(trialnum)
    counter = counter+1;
    subplot(rowVal,colVal,counter); hold on;

    plot(x_sample,y_sample,'Color',[.8 .8 .8]); box off; hold on;  
    title(['Trial ', num2str(trialnum(triali)), ' zlnIdPhi = ', num2str(shuffleIdPhi{1, trialnum(triali)})])

    set(gcf,'Position', (get(0, 'Screensize')));
    plot(shuffleIdPhi{2,trialnum(triali)}(1,:),shuffleIdPhi{2,trialnum(triali)}(2,:),'r', 'LineWidth', 2)
    rectangle ('position', rCP_fld);  % right reward field
    rectangle ('position', lCP_fld);  % left reward field

    if counter == 6 || triali == length(trialnum)
        addVTE = horzcat(addVTE, str2num(input('Enter trials that should be counted as VTE: ','s')));        %pause;
        close;
        set(gcf,'Position', (get(0, 'Screensize')));
        counter = 0;
    end
end

save('addVTE', 'addVTE')

%save all remove/include variables together
VTE_clean.step1 = idx2rem_step1;
VTE_clean.step2 = idx2rem_step2;
VTE_clean.addVTE = addVTE;
VTE_clean.info = ['data variable contains index of VTE trials to remove/include based on set threshold. ' ...
    'data(1,:) is a vector containg index of trials above zlnIdPhi threshold marked for removal after first glance through trials' ...
    'data(2,:) is a vector containig index of VTE trials to remove after watching position data (based on index from step 1).' ...
    'data(3,:) is a vector containing index of trials below zlnIdPhi threshold to include as VTE trials'];
save('VTE_clean', 'VTE_clean')
    

%% Return data to original structure
% add remove/inlude info to row in shuffleIdPhi variable
vteTrials = zeros(1,length(shuffleIdPhi));
vteTrials(idxVTE)=1;           %label trials above threshold
vteTrials(VTE_clean.step2)=0;  %remove falsely identified VTEs
vteTrials(VTE_clean.addVTE)=1; %add VTE trials with zlnIdPhi below threshold
find(vteTrials==1);      

%add row to shuffleVTE with VTE trial info: 1 = VTE trial, 0 = non-VTE trial
shuffleVTE = vertcat(shuffleIdPhi, num2cell(vteTrials));

%try to return data to original order
sortVTE = [];
for condi = 1:numel(unique(cell2mat(shuffleVTE(3,:))))
    %get all of AE (1) or SI (2) data
    idxCondi = []; condiIdPhi = [];
    idxCondi = find(cell2mat(shuffleVTE(3,:))==condi);
    condiIdPhi = shuffleVTE(:,idxCondi);

    for rati = 1:numel(unique(cell2mat(condiIdPhi(4,:))))
        %get all of rat n's data
        idxRat = []; ratIdPhi = [];
        idxRat   = find(cell2mat(condiIdPhi(4,:))==rati);
        ratIdPhi = condiIdPhi(:,idxRat);

        for sessi=1:numel(unique(cell2mat(ratIdPhi(5,:))))
            %get rat n's data for session i
            sessIdPhi = []; idxSess= [];
            idxSess   = find(cell2mat(ratIdPhi(5,:))==sessi);
            sessIdPhi    = ratIdPhi(:,idxSess);

            %order the trials
            idxTrial = [];
            [~,idxTrial] = sort(cell2mat(sessIdPhi(6,:)));
            sortVTE{condi}{rati,sessi} = sessIdPhi(:,idxTrial);
        end
    end
end %yay!

%forgot to exlclude this trial earlier (double trial NaN, need index to
%match with delay and accuracy cells for later...)
sortVTE{1}{8,14}(:,10)=[];
save('sortVTE', 'sortVTE')

%save VTE info in this structure
zlnIdPhi_struct = [];
zlnIdPhi_struct.AE.data = sortVTE{1};
zlnIdPhi_struct.AE.info = ['data variable contains zlnIdphi scores. Within each cell: ' ...
    'data(1,:) is zlnIdPhi data. data(2,:) is position data. data(3,:) is condition index (1=AE) data(4,:) is ' ...
    'rat index. data(5,:) is session index. data(6,:) is trial index. data(7,:) is VTE (1) vs non-VTE (0).'];
zlnIdPhi_struct.SI.data = sortVTE{2};
zlnIdPhi_struct.SI.info = ['data variable contains zlnIdphi scores. Within each cell:' ...
    'data(1,:) is zlnIdPhi data. data(2,:) is position data. data(3,:) is condition index (2=SI) data(4,:) is ' ...
    'rat index. data(5,:) is session index. data(6,:) is trial index. data(7,:) is VTE (1) vs non-VTE (0).'];


%% Cell array with vte, accuracy, delay, coherence, session, trial data
datafile={'data_analysis_aroundCP_22-Sep-2023.mat', 'sortVTE'};
for i = 1:numel(datafile)
   load(datafile{i}) 
end

% prep cells - maintain session structure
posAESI{1} = cell(size(posData{1}));
posAESI{2} = cell(size(posData{2}));

lfpAESI{1} = cell(size(posData{1}));
lfpAESI{2} = cell(size(posData{2}));

delayAESI{1} = cell(size(posData{1}));
delayAESI{2} = cell(size(posData{2}));

% use accBoolean for acc, index trials later

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

% Clean
% failed choice exit
delayAESI{2}{6,7}(35)= [];
accBoolean{2}{6,7}(35)= [];
posAESI{2}{6,7}(35)= [];
lfpAESI{2}{6,7}(35)= [];
% remove tagged double trials
lfpAESI{1}{1,20}(17)=[];
posAESI{1}{1,20}(17)=[];
posAESI{1}{8,14}(10)=[];

% Coherence 
cohAESI = [];
f = 1:0.5:50;
srate = 2000;
for condi = 1:length(lfpAESI)
    for rati = 1:size(lfpAESI{condi},1)
        for sessi= 1:size(lfpAESI{condi},2)
            if isempty(lfpAESI{condi}{rati,sessi})
                continue
            end

            for triali = 1:length(lfpAESI{condi}{rati,sessi})
                if isnan(lfpAESI{condi}{rati,sessi}{triali})
                    cohAESI{condi}{rati,sessi}{triali}= NaN;
                else
                    pfc = []; hpc = []; lfp = [];
                    lfp = lfpAESI{condi}{rati,sessi}{triali};
                    pfc = detrend(lfp(1,:),3);
                    hpc = detrend(lfp(2,:),3);

                    %get coherence from all trials
                    cohAESI{condi}{rati,sessi}{triali} = mscohere(pfc, hpc, [], [], f, srate);
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

% cell array sorted  by rat with zlnIdPhi, VTE, delay, choice accuracy,
% coherence, session, and trial data
vteAnalyze=[];
varNames = [{'zlnIdPhi'},{'VTE'},{'Delay'},{'Accuracy'},{'Coherence'},{'Session'},{'Trial'},{'Position'}];
for condi = 1:length(sortVTE)
    for rati = 1:size(sortVTE{condi},1)
        array={};
        for sessi = 1:size(sortVTE{condi},2)
            if isempty(sortVTE{condi}{rati,sessi})
                continue
            end

            temp = {};
            temp(:,1) = sortVTE{condi}{rati,sessi}(1,:); %get zlnIdPhi
            temp(:,2) = sortVTE{condi}{rati,sessi}(7,:); %get VTE 
            temp(:,3) = delayAESI{condi}{rati,sessi}; % get delay
            temp(:,4) = num2cell(accBoolean{condi}{rati,sessi}); %get trial choice accuracy
            try
                temp(:,5) = cohAESI{condi}{rati,sessi}; %get coherence
            catch
                temp(:,5) = num2cell(nan(size(temp,1),1)); %nan for behavior-only rats
            end
            temp(1:size(temp,1),6) = num2cell(sessionNumber{condi}(rati,sessi)); % get session/rat indicator 
            temp(:,7) = num2cell(1:size(temp,1)); % trial information
            temp(:,8) = posAESI{condi}{rati,sessi}; %get positon data

            %erase rows w/ NaN
            idxNan= isnan(cell2mat(temp(:,1))) | isnan(cell2mat(temp(:,3)));
            temp(idxNan,:)=[];

            % concatenate across loops, get each rat's trial info
            array = vertcat(array,temp);
        end
        % across rats
        vteAnalyze{condi}{rati} = vertcat(varNames, array);
    end
end

%save('vteAnalyze', 'vteAnalyze')

%% What rats have what sessions
sessionExist=[];
for condi= 1:length(sessionNumber)
    maxSession= max(sessionNumber{condi},[],'all');
    for rati= 1:size(sessionNumber{condi},1)
        for session= 1:maxSession
            if nnz(sessionNumber{condi}(rati,:)==session)==1
                sessionExist{condi}(session,rati)=1;
            else
                sessionExist{condi}(session,rati)=0;
            end
        end
    end
end

sessionExist{1}(:,8)=[]; %exclude rat 8
figure('Color','w')
sessionBar = bar(sessionExist{1},'stacked', FaceColor='flat')
sessionBar = bar(sessionExist{1},'stacked', FaceColor='flat')
sessionBar(1).CData = [0.75, 0, 0.75];
sessionBar(2).CData = [0, 0.75, 0.75];
sessionBar(3).CData = [0.75, 0.75, 0];
sessionBar(4).CData = [0.25, 0.25, 0.25];
sessionBar(5).CData = [0.6350, 0.0780, 0.1840];
sessionBar(6).CData = [0.3010, 0.7450, 0.9330];
sessionBar(7).CData = [0.4660, 0.6740, 0.1880];
sessionBar(8).CData = [0.4940, 0.1840, 0.5560];
sessionBar(9).CData = [0.9290, 0.6940, 0.1250];
sessionBar(10).CData = [0.8500, 0.3250, 0.0980];
ylim([0,10])
xlabel('Session Number')
ylabel('Rats w/ data')
legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5', 'Rat6', 'Rat7', 'Rat9', 'Rat10', 'Rat11')
title('sessions: behavior AE')

figure('Color','w')
sessionBar = bar(sessionExist{1}(:,1:7),'stacked', FaceColor='flat')
sessionBar = bar(sessionExist{1}(:,1:7),'stacked', FaceColor='flat')
sessionBar(1).CData = [0.75, 0, 0.75];
sessionBar(2).CData = [0, 0.75, 0.75];
sessionBar(3).CData = [0.75, 0.75, 0];
sessionBar(4).CData = [0.25, 0.25, 0.25];
sessionBar(5).CData = [0.6350, 0.0780, 0.1840];
sessionBar(6).CData = [0.3010, 0.7450, 0.9330];
sessionBar(7).CData = [0.4660, 0.6740, 0.1880];
ylim([0,10])
xlabel('Session Number')
ylabel('Rats w/ data')
legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5', 'Rat6', 'Rat7')
title('sessions: lfp AE')

figure('Color','w')
sessionBar = bar(sessionExist{2},'stacked', FaceColor='flat')
sessionBar(1).CData = [0.75, 0, 0.75];
sessionBar(2).CData = [0, 0.75, 0.75];
sessionBar(3).CData = [0.75, 0.75, 0];
sessionBar(4).CData = [0.25, 0.25, 0.25];
sessionBar(5).CData = [0.6350, 0.0780, 0.1840];
sessionBar(6).CData = [0.3010, 0.7450, 0.9330];
sessionBar(7).CData = [0.4660, 0.6740, 0.1880];
xlabel('Session Number')
ylabel('Rats w/ data')
ylim([0,10])
legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5', 'Rat6', 'Rat7')
title('sessions: behavior SI')


%% Set up VTE analysis variables: lnIdPhi, proportion VTE trials, VTE choice accuracy, VTE proportion by delay, VTE accuracy by delay
datafolder = '/Users/haileyrosenblum/Desktop/Matlab/Data output';
cd(datafolder);
datafile={'data_analysis_stem2cpExit_22-Sep-2023.mat' 'allAccData' 'vteAnalyze' 'delayMinMax'};
for i = 1:numel(datafile)
   load(datafile{i}) 
end

lnIdPhi=[]; propVTE=[];accVTE=[]; vteDelay1=[]; vteDelayAcc1 =[]; vteDelay2=[]; vteDelayAcc2=[]; numVTE=[]; numDelayVTE1=[]; numDelayVTE2=[];
for condi=1:length(vteAnalyze)
    for rati = 1:length(vteAnalyze{condi})

        idxIdPhi = find(contains(vteAnalyze{condi}{rati}(1,:),'zlnIdPhi'));
        lnIdPhi{condi}(rati) = mean(cell2mat(vteAnalyze{condi}{rati}(2:end,idxIdPhi))); %check lnIdPhi by rat

        idxVTE = find(contains(vteAnalyze{condi}{rati}(1,:),'VTE'));
        propVTE{condi}(rati) = mean(cell2mat(vteAnalyze{condi}{rati}(2:end,idxVTE))); %proportion VTE trials by rat

        idxAcc = find(contains(vteAnalyze{condi}{rati}(1,:), 'Accuracy'));
        idxVTE = find(contains(vteAnalyze{condi}{rati}(1,:),'VTE'));

         %find VTE trials, then find corresponding accuracy information -->
         %%get percent correct on VTE trials
        idxVTEtrials = find(cell2mat(vteAnalyze{condi}{rati}(2:end,idxVTE))==1)+1;
        accVTE{condi}(rati) = 1-mean(cell2mat(vteAnalyze{condi}{rati}(idxVTEtrials,idxAcc))); 

        numVTE{condi}(rati,1)=numel(idxVTEtrials); %curious to see actual number of vte
        numVTE{condi}(rati,2)=numel(cell2mat(vteAnalyze{condi}{rati}(2:end,idxVTE)));

        %Check VTEs by delay 
        delayTrials=[];
        idxDelay = find(contains(vteAnalyze{condi}{rati}(1,:), 'Delay'));

        %delay range method 1
        %get index of all trials corresponding to each delay
        idxShort=[]; idxMed=[]; idxLong=[];
        delayTrials = cell2mat(vteAnalyze{condi}{rati}(2:end,idxDelay));
        idxShort = find(delayTrials >= delayMinMax{1}(1,1) & delayTrials <= delayMinMax{1}(1,2))+1;
        idxMed   = find(delayTrials >= delayMinMax{1}(2,1) & delayTrials <= delayMinMax{1}(2,2))+1;
        idxLong  = find(delayTrials >= delayMinMax{1}(3,1) & delayTrials <= delayMinMax{1}(3,2))+1;

        %get vte trial information for each delay range
        shortVTE=[]; medVTE=[]; longVTE=[];
        shortVTE = cell2mat(vteAnalyze{condi}{rati}(idxShort,idxVTE))==1;
        medVTE   = cell2mat(vteAnalyze{condi}{rati}(idxMed, idxVTE))==1;
        longVTE  = cell2mat(vteAnalyze{condi}{rati}(idxLong, idxVTE))==1;
        
        %get proportion of vte trials at each delay
        vteDelay1{condi}(rati,1) = mean(shortVTE);
        vteDelay1{condi}(rati,2) = mean(medVTE);
        vteDelay1{condi}(rati,3) = mean(longVTE);

        %how many trials do we have for accuracy analysis?
        numDelayVTE1{condi}(rati,1)=nnz(shortVTE);
        numDelayVTE1{condi}(rati,2)=nnz(medVTE);
        numDelayVTE1{condi}(rati,3)=nnz(longVTE);
    
        %factor in accuracy to vte delay analysis
        accuracy = [];
        accuracy = cell2mat(vteAnalyze{condi}{rati}(2:end,idxAcc))==1;

        %vte percent correct at each delay
        shortAcc=[]; medAcc=[]; longAcc=[];
        shortAcc= accuracy(idxShort-1); medAcc= accuracy(idxMed-1); longAcc= accuracy(idxLong-1);
        vteDelayAcc1{condi}(rati,1)= mean(shortAcc(shortVTE));
        vteDelayAcc1{condi}(rati,2)= mean(medAcc(medVTE));
        vteDelayAcc1{condi}(rati,3)= mean(longAcc(longVTE));

        %exlude if rat has under 3 trials to contribute
        idxUnder3=[];
        idxUnder3 = find(numDelayVTE1{condi}(rati,:) <3);
        vteDelayAcc1{condi}(rati,idxUnder3)= NaN;

        %delay range method 2 (repeat as above)
        %get index of all trials corresponding to each delay
        idxShort=[]; idxMed=[]; idxLong=[];
        idxShort = find(delayTrials >= delayMinMax{2}(1,1) & delayTrials <= delayMinMax{2}(1,2))+1;
        idxMed = find(delayTrials >= delayMinMax{2}(2,1) & delayTrials <= delayMinMax{2}(2,2))+1;
        idxLong = find(delayTrials >= delayMinMax{2}(3,1) & delayTrials <= delayMinMax{2}(3,2))+1;

        %get vte trial information for each delay range
        shortVTE=[]; medVTE=[]; longVTE=[];
        shortVTE = cell2mat(vteAnalyze{condi}{rati}(idxShort, idxVTE))==1;
        medVTE   = cell2mat(vteAnalyze{condi}{rati}(idxMed, idxVTE))==1;
        longVTE  = cell2mat(vteAnalyze{condi}{rati}(idxLong, idxVTE))==1;

         %get proportion of vte trials at each delay
        vteDelay2{condi}(rati,1)= mean(shortVTE);
        vteDelay2{condi}(rati,2)= mean(medVTE);
        vteDelay2{condi}(rati,3)= mean(longVTE);

        %how many trials do we have for accuracy analysis?
        numDelayVTE2{condi}(rati,1)=nnz(shortVTE);
        numDelayVTE2{condi}(rati,2)=nnz(medVTE);
        numDelayVTE2{condi}(rati,3)=nnz(longVTE);
        
        %vte percent correct at each delay
        shortAcc=[]; medAcc=[]; longAcc=[];
        shortAcc= accuracy(idxShort-1); medAcc= accuracy(idxMed-1); longAcc= accuracy(idxLong-1);
        vteDelayAcc2{condi}(rati,1)= mean(shortAcc(shortVTE));
        vteDelayAcc2{condi}(rati,2)= mean(medAcc(medVTE));
        vteDelayAcc2{condi}(rati,3)= mean(longAcc(longVTE));   

        %exlude if rat has under 3 trials to contribute
        idxUnder3=[];
        idxUnder3 = find(numDelayVTE2{condi}(rati,:) <3);
        vteDelayAcc2{condi}(rati,idxUnder3)= NaN;
    end
end

%remove outlier rat 8
lnIdPhi{1}(8)=[];
propVTE{1}(8)=[];
accVTE{1}(8)=[];
vteDelay1{1}(8,:)=[];
vteDelayAcc1{1}(8,:)=[];
vteDelay2{1}(8,:)=[];
vteDelayAcc2{1}(8,:)=[];

%% IdPhi
data2plot = []; data2plot{1} = lnIdPhi{1}; data2plot{2} = lnIdPhi{2};
figure('color','w')
multiBarPlot(data2plot,[{'AE'} {'SI'}],'Mean zlnIdPhi','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Mean lnIdPhi',numCorrections);

%% Prop VTE
data2plot = []; data2plot{1} = propVTE{1}; data2plot{2} = propVTE{2};
figure('color','w')
multiBarPlot(data2plot,[{'AE'} {'SI'}],'Proportion VTE Per Rat','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Proportion VTE Per Rat',numCorrections);

%% Acc VTE 
data2plot = []; data2plot{1} = accVTE{1}; data2plot{2} = accVTE{2};
figure('color','w')
multiBarPlot(data2plot,[{'AE'} {'SI'}],'Percent correct on VTE trials Per Rat (%)','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Percent correct on VTE trials Per Rat',numCorrections);

%% VTE by Delay Method 1
%error bar
vteByDelayAE = mean(vteDelay1{1},1);
vteByDelaySI = mean(vteDelay1{2},1);
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = stderr(vteDelay1{1},1);
errSI = stderr(vteDelay1{2},1);
figure('Color', 'w'); hold on;
errorbar(group, vteByDelayAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, vteByDelaySI,errSI, 'b', 'LineWidth', 2);
ylim([0.04 0.22])
ylabel('Proportion VTE Trials')
title('Method 1, N=rats')

%bar 
data2plot = []; data2plot{1} = vteDelay1{1}(:,1); data2plot{2} = vteDelay1{2}(:,1); data2plot{3} = vteDelay1{1}(:,2); data2plot{4} = vteDelay1{2}(:,2);...
     data2plot{5} = vteDelay1{1}(:,3); data2plot{6} = vteDelay1{2}(:,3);
figure('color','w')
multiBarPlot(data2plot,[{'AE Short'} {'SI Short'} {'AE Med'} {'SI Med'} {'AE Long'} {'SI Long'}],'Proportion VTE Trials by Delay Method 1','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Proportion VTE Trials by Delay Method 1',numCorrections);
readStats(data2plot{3},data2plot{4},parametric,stat_test,'Proportion VTE Trials by Delay Method 1',numCorrections);
readStats(data2plot{5},data2plot{6},parametric,stat_test,'Proportion VTE Trials by Delay Method 1',numCorrections);

%Stacked bar-- how many trials do we have for accuracy analysis for each rat?
%plan: remove rat from accuracy analysis for the don't have at least 3 vte
%trials at each delay length
numDelayVTE1{1}(8,:)=[]; %exclude rat 8
figure('Color','w')
subplot 211; bar(numDelayVTE1{1},'stacked', FaceColor='flat');
ylabel('Number of VTE trials at each delay')
xlabel('Rat')
legend('short', 'mediuim', 'long')
subplot 212; bar(numDelayVTE1{2},'stacked', FaceColor='flat');
xlabel('Rat')
ylabel('Number of VTE trials at each delay')
legend('short', 'mediuim', 'long')

numDelayVTE2{1}(8,:)=[]; %exclude rat 8
figure('Color','w')
subplot 211; bar(numDelayVTE2{1},'stacked', FaceColor='flat');
ylabel('Number of VTE trials at each delay')
xlabel('Rat')
legend('short', 'mediuim', 'long')
subplot 212; bar(numDelayVTE2{2},'stacked', FaceColor='flat');
xlabel('Rat')
ylabel('Number of VTE trials at each delay')
legend('short', 'mediuim', 'long')


% vte choice accuracy by delay
data2plot = []; data2plot{1} = vteDelayAcc1{1}(:,1); data2plot{2} = vteDelayAcc1{2}(:,1); data2plot{3} = vteDelayAcc1{1}(:,2); data2plot{4} = vteDelayAcc1{2}(:,2);...
     data2plot{5} = vteDelayAcc1{1}(:,3); data2plot{6} = vteDelayAcc1{2}(:,3);
figure('color','w')
multiBarPlot(data2plot,[{'AE Short'} {'SI Short'} {'AE Med'} {'SI Med'} {'AE Long'} {'SI Long'}],'Percent correct on VTE trials Per Rat Method 1','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Percent correct on VTE trials Per Rat Method 1',numCorrections);
readStats(data2plot{3},data2plot{4},parametric,stat_test,'Percent correct on VTE trials Per Rat Method 1',numCorrections);
readStats(data2plot{5},data2plot{6},parametric,stat_test,'Percent correct on VTE trials Per Rat Method 1',numCorrections);


%% VTE by Delay Method 2
%error bar
vteByDelayAE = mean(vteDelay2{1},1);
vteByDelaySI = mean(vteDelay2{2},1);
group = categorical({'Short', 'Med', 'Long'}, {'Short', 'Med', 'Long'});
errAE = stderr(vteDelay2{1},1);
errSI = stderr(vteDelay2{2},1);
figure('Color', 'w'); hold on;
errorbar(group, vteByDelayAE, errAE, 'r', 'LineWidth', 2);
errorbar(group, vteByDelaySI,errSI, 'b', 'LineWidth', 2);
ylim([0.06 0.22])
ylabel('Proportion VTE Trials')
title('Method 2, N=rats')

%bar
data2plot = []; data2plot{1} = vteDelay2{1}(:,1); data2plot{2} = vteDelay2{2}(:,1); data2plot{3} = vteDelay2{1}(:,2); data2plot{4} = vteDelay2{2}(:,2);...
     data2plot{5} = vteDelay2{1}(:,3); data2plot{6} = vteDelay2{2}(:,3);
figure('color','w')
multiBarPlot(data2plot,[{'AE Short'} {'SI Short'} {'AE Med'} {'SI Med'} {'AE Long'} {'SI Long'}],'Proportion VTE Trials by Delay Method 2','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Proportion VTE Trials by Delay Method 2',numCorrections);
readStats(data2plot{3},data2plot{4},parametric,stat_test,'Proportion VTE Trials by Delay Method 2',numCorrections);
readStats(data2plot{5},data2plot{6},parametric,stat_test,'Proportion VTE Trials by Delay Method 2',numCorrections);


% vte error trials by delay
data2plot = []; data2plot{1} = vteDelayAcc2{1}(:,1); data2plot{2} = vteDelayAcc2{2}(:,1); data2plot{3} = vteDelayAcc2{1}(:,2); data2plot{4} = vteDelayAcc2{2}(:,2);...
     data2plot{5} = vteDelayAcc2{1}(:,3); data2plot{6} = vteDelayAcc2{2}(:,3);
figure('color','w')
multiBarPlot(data2plot,[{'AE Short'} {'SI Short'} {'AE Med'} {'SI Med'} {'AE Long'} {'SI Long'}],'Percent correct on VTE trials Per Rat Method 2','y')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Percent correct on VTE trials Per Rat Method 2',numCorrections);
readStats(data2plot{3},data2plot{4},parametric,stat_test,'Percent correct on VTE trials Per Rat Method 2',numCorrections);
readStats(data2plot{5},data2plot{6},parametric,stat_test,'Percent correct on VTE trials Per Rat Method 2',numCorrections);


%% Session choice accuracy and proportion VTE
%prep cells
propVTE{1} = nan(size(posData{1}));
propVTE{2} = nan(size(posData{2}));

sessAcc{1} = nan(size(posData{1}));
sessAcc{2} = nan(size(posData{2}));

sessionNumVTE{1} = nan(size(posData{1}));
sessionNumVTE{2} = nan(size(posData{2}));
for condi=1:length(vteAnalyze)
    for rati = 1:length(vteAnalyze{condi})
        idxVTE=[]; idxSess=[]; idxAcc=[]; numSess=[];
        %find column containing VTE, Session data
        idxVTE = find(contains(vteAnalyze{condi}{rati}(1,:),'VTE'));
        idxSess = find(contains(vteAnalyze{condi}{rati}(1,:),'Session'));
        idxAcc = find(contains(vteAnalyze{condi}{rati}(1,:),'Accuracy'));

        %which sessions does this rat have?
        numSess = unique(cell2mat(vteAnalyze{condi}{rati}(2:end, idxSess)));

        for sessi=1:length(numSess)
            %find rows corresponing to data from sessi
            idxSessN= find(cell2mat(vteAnalyze{condi}{rati}(2:end,idxSess))==numSess(sessi))+1;

            %get proportion VTE and choice accuracy from that session
            propVTE{condi}(rati,sessi)= mean(cell2mat(vteAnalyze{condi}{rati}(idxSessN,idxVTE)));
            sessAcc{condi}(rati,sessi)= 1-mean(cell2mat(vteAnalyze{condi}{rati}(idxSessN,idxAcc)));

            %how many VTE/session? try choice acc?
            sessionNumVTE{condi}(rati,sessi)=numel(find(cell2mat(vteAnalyze{condi}{rati}(idxSessN,idxVTE))==1));
        end
    end
end

%exclude rat 8
propVTE{1}(8,:)=[];
sessAcc{1}(8,:)=[];
accData{1}(8,:)=[];

sessVTE_AE=[]; sessVTE_SI=[]; sessAcc_AE=[]; sessAcc_SI=[]; 
sessVTE_AE = vertcat(propVTE{1}(:));
sessVTE_SI = vertcat(propVTE{2}(:));
sessAcc_AE = vertcat(sessAcc{1}(:));
sessAcc_SI = vertcat(sessAcc{2}(:));

sessVTE_AE(isnan(sessVTE_AE))=[];
sessVTE_SI(isnan(sessVTE_SI))=[];
sessAcc_AE(isnan(sessAcc_AE))=[];
sessAcc_SI(isnan(sessAcc_SI))=[];

%scatter plot
figure('Color', 'w'); hold on 
scatter(sessVTE_AE, sessAcc_AE, 'r')
scatter(sessVTE_SI, sessAcc_SI, 'b')
ls = lsline;
set(ls(1), 'color', 'b')
set(ls(2), 'color', 'r')
xlabel('Proportion VTE per session')
ylabel('Choice Accuracy (%)')

[r, p] = corrcoef(sessVTE_AE, sessAcc_AE)
[r, p] = corrcoef(sessVTE_SI, sessAcc_SI)

%% Proportion VTE (sessions)
data2plot = []; data2plot{1} = sessVTE_AE; data2plot{2} = sessVTE_SI;
figure('color','w')
multiBarPlot(data2plot,[{'AE'} {'SI'}],'Proportion VTE (N=sessions)','n')
stat_test = 'ttest2'; parametric = 'y'; numCorrections = [];
readStats(data2plot{1},data2plot{2},parametric,stat_test,'Proportion VTE (N=sessions)',numCorrections);

%% VTE and learning
% get each rat's VTEs by session
for condi=1:length(vteAnalyze)
    for rati=1:length(vteAnalyze{condi})




%data length for later
dist_ae = cellfun(@length,dataAE_lfp{1}(2:end,1));


