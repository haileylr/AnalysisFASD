%% plot coherence and power
% generate coherence and power figures, show significant frequencies (2 sample t test-> fdr_bh)
%
% INPUTS
% (2) Structure with 3 fields containing: (outputs of coh_pow.m for two groups)...
%      group1.coh  (nomalized coherence by session)
%      group1.pHPC (nomalized HPC power by session)
%      group1.pPFC (normalized PFC power by session)
%
%      group2.coh  (nomalized coherence by session)
%      group2.pHPC (nomalized HPC power by session)
%      group2.pPFC (normalized PFC power by session)
%
% f = frequencies as specified in coh_pow.m
% titleCoh = title for coherence figure 
% titlePow = title for power figure
%
% OUTPUTS
% coherence and power figures, significance noted


function cohPowFigures(group1, group2, f, titleCoh, titlePow)

pValue=[]; adj_p = []; idxsig =[];
for freq = 1:length(f)
    %2 sample t test
    [~,p] =ttest2(group1.coh(:,freq),group2.coh(:,freq));
    pValue(freq) =p;
end
[~,~,~,adj_p] = fdr_bh(pValue); %find adjusted p value

%index significant frequencies and prepare to plot
idxsig = adj_p < 0.05;
freqSig = f(idxsig);
stars = ones(1,length(freqSig));

% Plotting coherence
figure('color','w'); hold on 
shadedErrorBar(f,mean(group1.coh,'omitnan'),stderr(group1.coh,1),'r',0);
shadedErrorBar(f,mean(group2.coh,'omitnan'),stderr(group2.coh,1),'b',0);
plot(freqSig, stars, '*', 'MarkerSize',5)
%xlim([4 12])
xlabel('Frequency (Hz)');
ylabel('Normalized Coherence'); 
title(titleCoh)

%alpha correction power normalized
pValue =[]; adj_p= []; idxsig =[];
for freq = 1:length(f)
    %2 sample t test
    [~,p] =ttest2(group1.pPFC(:,freq),group2.pPFC(:,freq));
    pValue(freq) =p;
end
[~,~,~,adj_p] = fdr_bh(pValue); %find adjusted p value

%index significant frequencies and prepare to plot
idxsig = adj_p < 0.05;
freqSigPFC = f(idxsig);
starsPFC = ones(1,length(freqSigPFC));

pValue =[]; adj_p=[]; idxsig =[];
for freq = 1:length(f)
    %2 sample t test
    [~,p] =ttest2(group1.pHPC(:,freq),group2.pHPC(:,freq));
    pValue(freq) =p;
end
[~,~,~,adj_p] = fdr_bh(pValue); %find adjusted p value

%index significant frequencies and prepare to plot
idxsig = adj_p < 0.05;
freqSigHPC = f(idxsig);
starsHPC = ones(1,length(freqSigHPC));

% Plotting power
figure('color','w'); 
subplot 211; hold on;
shadedErrorBar(f,mean(group1.pPFC,'omitnan'),stderr(group1.pPFC,1),'r',0);
shadedErrorBar(f,mean(group2.pPFC,'omitnan'),stderr(group2.pPFC,1),'b',0);
plot(freqSigPFC, starsPFC, '*', 'MarkerSize',5)
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title(titlePow)

subplot 212; hold on; 
shadedErrorBar(f,mean(group1.pHPC,'omitnan'),stderr(group1.pHPC,1),'r',0);
shadedErrorBar(f,mean(group2.pHPC,'omitnan'),stderr(group2.pHPC,1),'b',0);
plot(freqSigHPC, starsHPC, '*', 'MarkerSize',5)
xlabel('Frequency (Hz)');
ylabel('Normalized Power');

end