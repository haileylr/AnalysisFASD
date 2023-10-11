%% Modulation index -- Phase amplitued coupling
% investigating teta-gamma coupling in sessions from alcohol exposed vs
% sham intubated rats

%based on random sampling, we found 4s of data was enough for reliable
%estimates of MI
load('data_analysis_inCP_07-Jul-2023.mat');  

%clean LFP
lfpAE = lfpData{1}(:);
lfpSI = lfpData{2}(:);

[~,idxAE] = emptyCellErase(lfpAE); 
[~,idxSI] = emptyCellErase(lfpSI);

lfpAE(idxAE)=[]; 
lfpSI(idxSI)=[]; 

for i = 1:length(lfpAE)
    lfpAE{i} = emptyCellErase(lfpAE{i});
end
for i = 1:length(lfpSI)
    lfpSI{i} = emptyCellErase(lfpSI{i});
end

%get MI
%AE rats
MI_alcohol = [];
signal_data.phase_bandpass     = ([6 12]);  %Frequency band for filtering the phase signal
signal_data.amplitude_bandpass = ([30 55]); %Frequency band for filtering the amplitude signal
signal_data.srate              = 2000;      %Sampling rate (Hz)

for sessioni= 1:length(lfpAE)
    pfc =[]; hpc = [];
    for triali= 1:length(lfpAE{sessioni})
        pfc{triali} = detrend(lfpAE{sessioni}{triali}(1,:),3);
        hpc{triali} = detrend(lfpAE{sessioni}{triali}(2,:),3);
    end

    %concatenate lfp across trials
    pfc_session=[]; hpc_session= [];
    pfc_session = horzcat(pfc{:});
    hpc_session = horzcat(hpc{:});

    %only include sessions w/ at least x seconds of data
    if length(pfc_session) < 60000
        continue
    end

    % Filter phase and amplitude signals
    fP_filt = []; fA_filt = [];
    fP_filt = skaggs_filter_var(hpc_session,signal_data.phase_bandpass(:,1),signal_data.phase_bandpass(:,2),signal_data.srate);
    fA_filt = skaggs_filter_var(pfc_session,signal_data.amplitude_bandpass(:,1),signal_data.amplitude_bandpass(:,2),signal_data.srate);

    %Hilbert transform -- extract phase and amplitude filtered signals
    fP_hilbert =[]; fP_phase =[]; fP_amp = []; fA_hilbert = []; fA_phase = []; fA_amp = [];
    fP_phase= hilbertPhase(fP_filt);
    fP_hilbert = hilbert(fP_filt);
    fP_amp = abs(fP_hilbert);

    fA_phase= hilbertPhase(fA_filt);
    fA_hilbert = hilbert(fA_filt);
    fA_amp = abs(fA_hilbert);

    %Bin by phase and find mean amplitude in each bin
    N = 18;  %number of bins (Tort 2010)
    PhaseBins = 0 : 360/N : 360-(360/N);
    Bin= []; index = [];
    for k = 1:N
        % Finds all points where phase falls within a given bin
        index = find (fP_phase >= PhaseBins(k) & fP_phase < PhaseBins(k)+(360/N));

        % mean amplitude within the designated phase bins
        Bin(k).amp = mean(fA_amp(index));
    end

    %noramlize mean ampliude
    S =[]; Bin(k).NormAmp =[];
    S = sum([Bin.amp]);
    for k = 1:N
        Bin(k).NormAmp = Bin(k).amp / S;
    end

    %Kullback-Leibler distance: quantifies the deviation of P from the uniform distribution
    %D(P,Q)= sum [ P(j) * log(P(j)/Q(j))]for each bin of distribution
    P = [];
    Q(1:N) = 1/N;
    P = [Bin.NormAmp];
    KLdist = sum (P.*(log(P./Q)));

    %Modulation Index (MI) = KL distance divided by the log of the number of bins
    MI_alcohol(sessioni) = KLdist/log(N);

    disp(['Completed with MI ',num2str(sessioni),'/',num2str(length(lfpAE))])
end

%SI rats
MI_sham= [];
for sessioni= 1:length(lfpSI)
    pfc =[]; hpc = [];
    for triali= 1:length(lfpSI{sessioni})
        pfc{triali} = detrend(lfpSI{sessioni}{triali}(1,:),3);
        hpc{triali} = detrend(lfpSI{sessioni}{triali}(2,:),3);
    end

    %concatenate lfp across trials
    pfc_session=[]; hpc_session= [];
    pfc_session = horzcat(pfc{:});
    hpc_session = horzcat(hpc{:});

    %only include sessions w/ at least x seconds of data
    if length(pfc_session) < 60000
        continue
    end

    % Filter phase and amplitude signals
    fP_filt = []; fA_filt = [];
    fP_filt = skaggs_filter_var(hpc_session,signal_data.phase_bandpass(:,1),signal_data.phase_bandpass(:,2),signal_data.srate);
    fA_filt = skaggs_filter_var(pfc_session,signal_data.amplitude_bandpass(:,1),signal_data.amplitude_bandpass(:,2),signal_data.srate);

    %Hilbert transform -- extract phase and amplitude filtered signals
    fP_hilbert =[]; fP_phase =[]; fP_amp = []; fA_hilbert = []; fA_phase = []; fA_amp = [];
    fP_phase= hilbertPhase(fP_filt);
    fP_hilbert = hilbert(fP_filt);
    fP_amp = abs(fP_hilbert);

    fA_phase= hilbertPhase(fA_filt);
    fA_hilbert = hilbert(fA_filt);
    fA_amp = abs(fA_hilbert);

    %Bin by phase and find mean amplitude in each bin
    N = 18;  %number of bins (Tort 2010)
    PhaseBins = 0 : 360/N : 360-(360/N);
    Bin= []; index = [];
    for k = 1:N
        % Finds all points where phase falls within a given bin
        index = find (fP_phase >= PhaseBins(k) & fP_phase < PhaseBins(k)+(360/N));

        % mean amplitude within the designated phase bins
        Bin(k).amp = mean(fA_amp(index));
    end

    %noramlize mean ampliude
    S =[]; Bin(k).NormAmp =[];
    S = sum([Bin.amp]);
    for k = 1:N
        Bin(k).NormAmp = Bin(k).amp / S;
    end

    %Kullback-Leibler distance: quantifies the deviation of P from the uniform distribution
    %D(P,Q)= sum [ P(j) * log(P(j)/Q(j))]for each bin of distribution
    P = [];
    Q(1:N) = 1/N;
    P = [Bin.NormAmp];
    KLdist = sum (P.*(log(P./Q)));

    %Modulation Index (MI) = KL distance divided by the log of the number of bins
    MI_sham(sessioni) = KLdist/log(N);

    disp(['Completed with MI ',num2str(sessioni),'/',num2str(length(lfpSI))])
end

%generate bar graph modulation index for AE vs SI
data2plot = []; data2plot{1} = MI_sham; data2plot{2} = MI_alcohol;
figure('Color', 'w');
multiBarPlot(data2plot, [{'SI'}, {'AE'}], 'Modulation Index (min 30s LFP)')
[h,p] = ttest2(MI_sham, MI_alcohol)

%Generate plot of normalized gamma amplitude as a function of phase (two full cycles)
M.PhaseAxis = 360/N/2 : 360/N : 360-(360/N/2);

figure1=figure('Color', 'w');
axes1 = axes('Parent',figure1,'YTick',[0 0.02 0.04 0.06 0.08],'XTick',[0 180 360 540 720]);
xlim(axes1, [0 720]);
ylim(axes1, [0 0.08]);
hold(axes1,'all')

% phase=[M.PhaseAxis 360+M.PhaseAxis];
% ampSI=[Bin.NormAmp Bin.NormAmp]; %max MI from session 55
% ampAE=[Bin.NormAmp Bin.NormAmp]; %max MI from session 14

subplot 211; bar(phase,ampSI,'BarWidth',1, 'FaceColor', 'b');
xlabel('Theta phase')
ylabel('Normalized gamma amplitude')
title('PAC Sham')

subplot 212; bar(phase,ampAE,'BarWidth',1, 'FaceColor', 'r');
xlabel('Theta phase')
ylabel('Normalized gamma amplitude')
title('PAC Alcohol')


%phase comodologram: function creates a co-modulogram of normalized amplitude values
%   across phase bins

%testing -- find example signal - session w/ highest MI
maxAE  = max(MI_alcohol);
sessAE = find(MI_alcohol == maxAE);
maxSI  = max(MI_sham);
sessSI = find(MI_sham == maxSI);

%set up for MI calculations, co modulogram
signal_data.srate              = 2000;
signal_data.phase_extraction   = 2;    % 0= Hilbert transform, 2 = Morlet wavelet
phase_bins                     = 18;

amplitude_freq_bins = 1;
phase_freq_bins     = 1;
plot                = 1;

avg_comodAE = [];
%AE: get phase signal, amplitude signal, and timestamps, then get average co-modulogram 
% values across all frequencies for phase for each session 
for sessioni= 1:length(lfpAE)
    pfc =[]; hpc = []; time = [];
    for triali= 1:length(lfpAE{sessioni})
        pfc{triali} = detrend(lfpAE{sessioni}{triali}(1,:),3);
        hpc{triali} = detrend(lfpAE{sessioni}{triali}(2,:),3);
        time{triali}= lfpAE{sessioni}{triali}(3,:);
    end

    %concatenate lfp and timestamps across trials
    pfc_session =[]; hpc_session = []; timestamps = [];
    pfc_session = horzcat(pfc{:});
    hpc_session = horzcat(hpc{:});
    timestamps  = horzcat(time{:});
    if length(pfc_session) < 8000
        continue
    end

    %%make for loop to get MI info from each sesssion
    signal_data.timestamps = []; signal_data.phase_EEG = []; signal_data.amplitude_EEG = [];
    signal_data.phase_bandpass =[]; signal_data.amplitude_bandpass =[];

    signal_data.timestamps         = timestamps;  % array of timestamp values
    signal_data.phase_EEG          = hpc_session; % data from signal containing phase information
    signal_data.amplitude_EEG      = pfc_session; % data from signal containing amplitude information
    signal_data.phase_bandpass     = ([6 12]);    % Frequency band for filtering the phase signal
    signal_data.amplitude_bandpass = ([15 100]);  % Frequency band for filtering the amplitude signal
    
    phase_map_var = [];
    [phase_map_var,M,amplitude_highpass] = phase_comodgram(signal_data, phase_bins, amplitude_freq_bins, phase_freq_bins, plot);
    
    %store average co-modulogram values from each session
    avg_comodAE(:,:, end+1) = phase_map_var;
    disp(['Completed with MI ',num2str(sessioni),'/',num2str(length(lfpAE))])
end

%SI: get phase signal, amplitude signal, and timestamps 
for sessioni= 1:length(lfpSI)
    pfc =[]; hpc = []; time = [];
    for triali= 1:length(lfpSI{sessioni})
        pfc{triali} = detrend(lfpSI{sessioni}{triali}(1,:),3);
        hpc{triali} = detrend(lfpSI{sessioni}{triali}(2,:),3);
        time{triali}= lfpSI{sessioni}{triali}(3,:);
    end

    %concatenate lfp and timestamps across trials
    pfc_session =[]; hpc_session = []; timestamps = [];
    pfc_session = horzcat(pfc{:});
    hpc_session = horzcat(hpc{:});
    timestamps  = horzcat(time{:});
    if length(pfc_session) < 8000
        continue
    end

    %%make for loop to get MI info from each sesssion
    signal_data.timestamps = []; signal_data.phase_EEG = []; signal_data.amplitude_EEG = [];
    signal_data.phase_bandpass =[]; signal_data.amplitude_bandpass =[];

    signal_data.timestamps         = timestamps;  % array of timestamp values
    signal_data.phase_EEG          = hpc_session; % data from signal containing phase information
    signal_data.amplitude_EEG      = pfc_session; % data from signal containing amplitude information
    signal_data.phase_bandpass     = ([6 12]);    % Frequency band for filtering the phase signal
    signal_data.amplitude_bandpass = ([15 100]);  % Frequency band for filtering the amplitude signal
    
    phase_map_var = [];
    [phase_map_var,M,amplitude_highpass] = phase_comodgram(signal_data, phase_bins, amplitude_freq_bins, phase_freq_bins, plot);
    
    %store average co-modulogram values from each session
    avg_comodSI(:,:, end+1) = phase_map_var;
    disp(['Completed with MI ',num2str(sessioni),'/',num2str(length(lfpSI))])
end

%get average co-mod values across sessions
phase_mapAE = [];
phase_mapAE = mean(avg_comodAE,3);
phase_mapSI = [];
phase_mapSI = mean(avg_comodSI,3);
    
%This is useful in setting the limits of the colorbar, scaled to SI
bottom = min(min(phase_mapSI));
top    = max(max(phase_mapSI));
%bottom = min(min(min(phase_mapAE)),min(min(phase_mapSI)));
%top    = max(max(max(phase_mapAE)),max(max(phase_mapSI)));

%Create phase co-modulogram
figure('Color', 'w')
tiledlayout(1,2)

nexttile
pcolor(M.PhaseAxis,amplitude_highpass,phase_mapAE); hold on
colormap(jet)
shading 'interp'
ylabel('mPFC Frequency for Amplitude (Hz)')
xlabel('HPC Theta Phase')
title('AE')
c = colorbar;
caxis manual
caxis([bottom top]);
c.Label.String = 'Modulation Index';

nexttile
pcolor(M.PhaseAxis,amplitude_highpass,phase_mapSI)
colormap(jet)
shading 'interp'
ylabel('Frequency for Amplitude (Hz)')
xlabel('HPC Theta Phase')
title('SI')
c = colorbar;
caxis manual
caxis([bottom top]);
c.Label.String = 'Modulation Index';

 