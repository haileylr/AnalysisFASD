%% coherence and power
% calculate coherence by trial or concatenated across sessions
% specify minimum data length to be used for analysis 
% currently using mscohere and pwelch, will need to edit if want to use
% cohherencyc/mtspectrumc
% 
% INPUTS: 
% lfp            = x by 1 cell containing trial data from x sessions of lfp (cleaned)
% f              = frequencies specified as a row or column vector with at least two elements
% srate          = sampling rate
% trialORsession = calculate by trial or concatenate across trials (0 =
%                  trial, 1 = session)
% datalength     = minimum length of lfp data for analysis by session (ex if you
%                  want at least 30 seconds of data per session, then
%                  datalength = 60000). Sessions below this length will be excluded 
%
% 
% OUTPUTS:
% structure containing:
%        coh  (nomalized coherence by session)
%        pHPC (nomalized HPC power by session)
%        pPFC (normalized PFC power by session)

function [output] = coh_pow(lfp, f, srate, trialORsession, datalength)

coh = []; powerPFC = []; powerHPC =[];
if trialORsession == 0
    for sessioni= 1:length(lfp)
        for triali= 1:length(lfp{sessioni})
            pfc = []; hpc = [];
            pfc = detrend(lfp{sessioni}{triali}(1,:),3);
            hpc = detrend(lfp{sessioni}{triali}(2,:),3);

            %only include signals w/ at least 1 sec data
            if length(pfc)<2000
                continue
            end

            %calculate coherence
            coh{sessioni}(triali,:)  = mscohere(pfc,hpc,[],[],f,srate);
            %coh{sessioni}(triali,:) =coherencyc(pfc,hpc,params);

            %calculate power
            S1 = []; S2 = [];
            S1 = pwelch(pfc,[],[],f,srate);
            S2 = pwelch(hpc,[],[],f,srate);
            % S1      = mtspectrumc(pfc,params);
            % [S2,Sf] = mtspectrumc(hpc,params);

            % log10 transform
            powerPFC{sessioni}(triali,:) = log10(S1);
            powerHPC{sessioni}(triali,:) = log10(S2);
        end
    end

    %loop over sessions and get coherence and power averages (get
    %averages at each frequency over trials in sessions), then concatenate
    %across sessions
    coh_sess = cellfun2(coh,'mean',{'1'});
    coh_mat  = vertcat(coh_sess{:});

    pPFC_sess = cellfun2(powerPFC,'mean',{'1'});
    pPFC_mat  = vertcat(pPFC_sess{:});

    pHPC_sess = cellfun2(powerHPC,'mean',{'1'});
    pHPC_mat  = vertcat(pHPC_sess{:});

    %normalize
    cohNor  = normalize(coh_mat, 2, 'range');
    pPFCnor = normalize(pPFC_mat, 2, 'range');
    pHPCnor = normalize(pHPC_mat, 2, 'range');

elseif trialORsession == 1
    for sessioni= 1:length(lfp)
        pfc =[]; hpc = [];
        for triali= 1:length(lfp{sessioni})
            pfc{triali} = detrend(lfp{sessioni}{triali}(1,:),3);
            hpc{triali} = detrend(lfp{sessioni}{triali}(2,:),3);
        end

        %concatenate lfp across trials
        pfc_session = []; hpc_session= [];
        pfc_session = horzcat(pfc{:});
        hpc_session = horzcat(hpc{:});

        %only include sessions w/ at least datalength/srate seconds of data
        if length(pfc_session) < datalength
            continue
        end

        %calculate coherence
        coh{end+1} = mscohere(pfc_session, hpc_session,[],[],f,srate);

        % calculate power
        S1 = []; S2 = [];
        S1 = pwelch(pfc_session,[],[],f,srate);
        S2 = pwelch(hpc_session,[],[],f,srate);

        % % log10 transform
        powerPFC{end+1} = log10(S1);
        powerHPC{end+1} = log10(S2);
    end

    %concatenate across sessions
    coh_mat  = vertcat(coh{:});
    pPFC_mat = vertcat(powerPFC{:});
    pHPC_mat = vertcat(powerHPC{:});

    %normalized coherence and power
    cohNor  = normalize(coh_mat, 2, 'range');
    pHPCnor = normalize(pHPC_mat, 2, 'range');
    pPFCnor = normalize(pPFC_mat, 2, 'range');

end

output.coh  = cohNor;
output.pHPC = pHPCnor;
output.pPFC = pPFCnor;

end




