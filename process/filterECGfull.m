% Author: M.J. Stam
% Date: 27-08-2024Linenoise_ECG
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to manually apply the Singular Value Decomposition ECG 
% suppression method by B.C.M. van Wijk, as described in Stam et al. 2023 
% (see reference below) to the BrainSenseTimeDomain Streaming data. Plots
% intermediate steps and allows to manually select R-peaks and SVD
% components for ECG suppression. 
% 
% Note for peak-detection: this script does not check whether cleaning 
% is already performed for specific recording (in order to use (and not 
% lose) the previously done (manual) R-peak detection).
%
% Stam MJ, van Wijk BCM, Sharma P, Beudel M, Pina-Fuentes DA, de Bie RMA, 
% Schuurman PR, Neumann W-J, Buijink AWG. A comparison of methods to 
% suppress electrocardiographic artifacts in local field potential 
% recordings. Clin Neurophysiol. 2023 Feb:146:147-161. DOI: 
% 10.1016/j.clinph.2022.11.011
%
% INPUT
%   info            =   struct containing metadata of recording
%   dataRec         =   struct containing info and BrainSenseTimeDomain
%                       Streaming time array, and LFP data in raw form of
%                       one recording.
%   rWindow         =   1x2 array defining window around R-peak [before 
%                       after] for calculating SVD components, in seconds.
%   savename_json   =   string of name under which to save figures
%   c               =   int indicating number of recording
%
% OUTPUT
%   dataRec         =   struct containing info and BrainSenseTimeDomain
%                       Streaming time array, and LFP data in raw, notch, 
%                       ECG and notch + ECG filtered form

function dataRec = filterECGfull(info, dataRec, rWindow, savename_json, c)

    % Location of folder where all figures are saved:
    savepathfig = [info.OriginalFolder [savename_json '_unpacked'] filesep 'Figures' filesep 'CleaningFigures'];
    if ~exist(savepathfig, 'dir')
        mkdir(savepathfig)
    end
    
    % Pwelch for plotting PSDs cannot be applied to NaNs, so this finds the start 
    % and end indices of segments of the Streaming data without missing values
    i_seg = [1; find(diff(any(isnan(dataRec.Data.LFP(:,1)),2)))+1; length(dataRec.Data.LFP)];
    
    % NaNs if ECG is not suppressed:
    dataRec.Data.LFP_Ecg = NaN(size(dataRec.Data.LFP,1), size(dataRec.Data.LFP,2));
    dataRec.Data.LFP_LinenoiseEcg = NaN(size(dataRec.Data.LFP,1), size(dataRec.Data.LFP,2));
    time = dataRec.Data.Time;
    
    % Initiate segment number:
    s = 0;
    
    % Loop over segments
    for int = 1:2:length(i_seg)
    
        % Find indices of particular segment
        % Because of function "diff" middle segments run to one sample too much
        if int == length(i_seg)
            ix = i_seg(int):1:i_seg(int+1);
        else
            ix = i_seg(int):1:i_seg(int+1)-1;
        end
    
        % To save info about the specific segment:
        s = s+1;
    
        % Check if recording contains data for left and right hemisphere
        chL=find(contains(dataRec.Data.Channel,'LEFT')==1);
        if ~isempty(chL)
            LFP_L{s}=dataRec.Data.LFP(ix,chL);
            LFP_L_linenoise{s} = dataRec.Data.LFP_Linenoise(ix,chL);
        end
    
        chR=find(contains(dataRec.Data.Channel,'RIGHT')==1);
        if ~isempty(chR)
            LFP_R{s}=dataRec.Data.LFP(ix,chR);
            LFP_R_linenoise{s} = dataRec.Data.LFP_Linenoise(ix,chR);
        end
    end
    
    %%
    
    % --- SVD-ECG Suppression to raw LFPs (output: lfp_ecg)
    
    % Initiate the general settings of SVD method
    
    % Generate SVD components of specific time window surrounding R-peaks
    % Default PQRST window is 250 ms before and 400 ms after R-peak
    % Alternatively use smaller QRS time window (200 ms before and after peaks)
    pre_R_time = rWindow(1); 
    post_R_time = rWindow(2); 
    
    % Set interactive = 1, to allow manual peak detection:
    settings.Interactive            = 1; 
    settings.ArtTimeBeforePeak      = pre_R_time;
    settings.ArtTimeAfterPeak       = post_R_time;
    settings.SamplingFrequency      = dataRec.Settings.SamplingFrequencyTimeDomain;
    settings.ShowFigs               = 1;
    settings.SaveFigs               = 1;
    
    % SVD ECG suppression is applied to LFP_L and LFP_R separately
    % For left hemisphere:
    if ~isempty(chL)
    
        % ECG suppression starts with an R-peak detection algorithm
    
        % Set threshold for minimal ECG R-peak height
        % Threshold is set according to the part of the signal before a 
        % potential artifact (> 5 std deviation) occurs (mx):
        upthres = mean(dataRec.Data.LFP(:,chL), "omitmissing")+(5*std(dataRec.Data.LFP(:,chL), "omitmissing"));
        lowthres = mean(dataRec.Data.LFP(:,chL), "omitmissing")-(5*std(dataRec.Data.LFP(:,chL), "omitmissing"));
    
        % Convert to string of zeros and ones
        Y = dataRec.Data.LFP(:,chL)>lowthres&dataRec.Data.LFP(:,chL)<upthres;
        idx=find(Y<1);
        if ~isempty(find(Y<1))
            idrange_full = [0; idx; length(dataRec.Data.LFP(:,chL))];
            idrange_start = find(diff(idrange_full)==max(diff(idrange_full)));
            mx = idrange_full(idrange_start)+1:idrange_full(idrange_start+1)-1;
        else
            mx = 1:length(dataRec.Data.LFP(:,chL));
        end
    
        % Initiate positive polarity of the R-peaks
        flagpeak_left = 1;

        % Normalize and find peaks for LFP_L
        LFPnorm_left = normalize(dataRec.Data.LFP(:,chL));

        % Threshold of 2 * sd (instead of 2.5 * sd as described in Stam et al. (2023) 
        % https://doi.org/10.1016/j.clinph.2022.11.011) for peak detection:
        [Rpeak_left, locs_Rwave_left] = findpeaks(LFPnorm_left, 'MinPeakHeight', (2 * std(LFPnorm_left(mx))), 'MinPeakDistance', dataRec.Settings.SamplingFrequencyTimeDomain / 2);
        [Speak_left, locs_Swave_left] = findpeaks(-LFPnorm_left, 'MinPeakHeight', (2 * std(-LFPnorm_left(mx))), 'MinPeakDistance', dataRec.Settings.SamplingFrequencyTimeDomain / 2);
        if isempty(Rpeak_left)
            Rpeak_left = 0;
        end
        if isempty(Speak_left)
            Speak_left = 0;
        end

        % Flip polarity if the mean height and the number of detected R-peaks
        % are higher with negative polarity
        if mean(Speak_left) > mean(Rpeak_left) && length(Speak_left) >= length(Rpeak_left)
            locs_Rwave_left = locs_Swave_left;
            locs_Rwave_left = unique(locs_Rwave_left); locs_Rwave_left = locs_Rwave_left(find(locs_Rwave_left>0));
            LFPnorm_left = -LFPnorm_left;
            flagpeak_left = -1;
        end
    
        % Check ECG artifact consistency for the left hemisphere and decide on ECG cleaning y/n:
        if isempty(locs_Rwave_left)
            key = 'n';
        else
            % Make sure locs_Rwave is an horizontal array for manual peak selection
            if size(locs_Rwave_left,1)>size(locs_Rwave_left,2)
                locs_Rwave_left = locs_Rwave_left';
            end
    
            % Plot LFP signal and detected R-peaks
            fig = figure('Name','Detected R-peaks and PSD','NumberTitle','off', 'visible', 'off');clf();
            fig.WindowState = 'maximized';
            sgtitle([savename_json ' LFP_L Raw - Recording ' num2str(c)'], 'Interpreter', 'none')
            subplot(2,s,1:s)
            plot(time,LFPnorm_left,  'k'), hold on,plot(time(locs_Rwave_left), LFPnorm_left(locs_Rwave_left),'r*', 'markersize', 15)
            xlabel('Time [s]', 'FontSize',12)
            box off
            legend('LFPs', ['Detected R-peaks: ' num2str(length(locs_Rwave_left))], 'FontSize',11, 'location', 'northwest')
            legend('boxoff')
    
            % A few pre-requisites before the SVD ECG suppression method is run:
                % If there are no R-peaks detected at all:
            if isempty(locs_Rwave_left)
                title('Automatically detected R-peaks: no obvious R-peaks', 'Interpreter', 'none', 'FontSize', 14);
                % If heartbeats are missed: no ECG artifact is detected:
            elseif any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                title('Automatically detected R-peaks: no R-peak every three seconds', 'Interpreter', 'none', 'FontSize', 14);
                % If heart rate is too low:
            elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain)
                title('Automatically detected R-peaks: unreliable heart rate (< 40 bpm)', 'Interpreter', 'none', 'FontSize', 14);
            else
                title(['Automatically detected R-peaks: consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'], 'Interpreter', 'none', 'FontSize', 14);
            end
    
            for k = 1:s
                if length(LFP_L{k}) > dataRec.Settings.SamplingFrequencyTimeDomain-1
                    [P,F]=pwelch(LFP_L{k},[],[],[],250);
                    [P1s,F1s]=pwelch(LFP_L{k},hamming(1*250),[],[],250);
                    subplot(2,s,s+k)
                    plot(F,P);hold on,plot(F1s,P1s,'r','linewidth',2);ylim([0 40]),xlim([0 130]);
                    xlabel('Frequency [Hz]', 'FontSize', 12); legend('PSD large window','PSD 1s window', 'FontSize', 11)
                    if k == 1
                        ylabel('Power', 'FontSize', 12)
                    end
                    title(['Segment ' num2str(k)], 'Interpreter', 'none', 'FontSize', 14);
                    box off
                else
                    subplot(2,s,s+k)
                    title(['Segment ' num2str(k) ' < 250 samples'], 'FontSize', 14')
                end
            end
    
            fig.Visible = 'on';
            savefig([savepathfig filesep [savename_json '_LFP_L_raw_rec' num2str(c) '_PSD_precleaning.fig']])
            print(gcf,'-djpeg', [savepathfig filesep [savename_json '_LFP_L_raw_rec' num2str(c) '_PSD_precleaning.jpg']]);
    
            key=input('Inspect detected R-peaks and PSDs of LFP-Left: Press "y" to apply ECG cleaning, press "n" to proceed without ECG cleaning:','s');
    
        end
    
        if key == 'n'
            if isempty(locs_Rwave_left)
                dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No obvious R-peaks';
                disp('ECG not applied to LFP_L: No obvious R-peaks');
            elseif any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No R-peak every three seconds';
                disp('ECG not applied to LFP_L: No R-peak every three seconds');
                close('Detected R-peaks and PSD')
            elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
                disp('ECG not applied to LFP_L: Unreliable heart rate (< 40 bpm)');
                close('Detected R-peaks and PSD')
            else
                dataRec.Data.EcgArtifactStatus{1,chL} = ['Data not ECG-suppressed despite consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'];
                disp(['ECG not applied to LFP_L despite consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)']);
                close('Detected R-peaks and PSD')
            end
            dataRec.Data.EcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.EcgArtifactStatus{3,chL} = [];
            dataRec.Data.EcgArtifactStatus{4,chL} = [];

        elseif key == 'y'
            close('Detected R-peaks and PSD')
            if any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.EcgArtifactStatus{1,chL} = 'Data ECG-suppressed despite no R-peak every three seconds';
            elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.EcgArtifactStatus{1,chL} = 'Data ECG-suppressed despite unreliable heart rate (< 40 bpm)';
            else
                dataRec.Data.EcgArtifactStatus{1,chL} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'];
            end
    
            dataRec.Data.EcgArtifactStatus{2,chL} = locs_Rwave_left;
    
            % Initiate LFP-specific settings of SVD method:
            % Initiated using 2 components, but manually select
            % specific components according to plot
            settings.Polarity               = flagpeak_left;
            settings.nComp                  = 2;
            settings.nCompFinal             = 1;

            % Location and name to save figures:
            settings.Folder                 = savepathfig;
            settings.Label                  = [savename_json '_LFP_L_raw_rec' num2str(c)];
            settings.Title                  = [savename_json ' LFP_L Raw - Recording ' num2str(c)];
    
            % Apply the SVD ECG suppression method to LFP_L and save filtered signal in LFP_Ecg:
            [dataRec.Data.LFP_Ecg(:,chL), settings.ProjOut, settings.FinalComp, locs_Rwave_left_final] = SVD_ECG_Filter(dataRec.Data.LFP(:,chL),settings,[],locs_Rwave_left,[]);

            % Save info about detected R-peaks to the status struct field
            if isempty(locs_Rwave_left_final)

                % If only one or two R-peaks are detected (very short LFP
                % signal) which are too close to the start and/or end, the 
                % (P)QRS(T) time window can not be created:
                dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: R-peaks too close to signal edges';
                dataRec.Data.EcgArtifactStatus{3,chL} = [];
                dataRec.Data.EcgArtifactStatus{4,chL} = [];
            else
                dataRec.Data.EcgArtifactStatus{3,chL} = locs_Rwave_left_final;
                dataRec.Data.EcgArtifactStatus{4,chL} = settings;
            end
            disp('ECG applied to LFP_L');
        end
    end
    
    % For right hemisphere:
    if ~isempty(chR)
    
        % ECG suppression starts with an R-peak detection algorithm
    
        % Set threshold for minimal ECG R-peak height
        % Threshold is set according to the part of the signal before a 
        % potential artifact (> 5 std deviation) occurs (mx):
        upthres = mean(dataRec.Data.LFP(:,chR), "omitmissing")+(5*std(dataRec.Data.LFP(:,chR), "omitmissing"));
        lowthres = mean(dataRec.Data.LFP(:,chR), "omitmissing")-(5*std(dataRec.Data.LFP(:,chR), "omitmissing"));
    
        % Convert to string of zeros and ones
        Y = dataRec.Data.LFP(:,chR)>lowthres&dataRec.Data.LFP(:,chR)<upthres; %// convert to string of zeros and ones
        idx=find(Y<1);
        if ~isempty(find(Y<1))
            idrange_full = [0; idx; length(dataRec.Data.LFP(:,chR))];
            idrange_start = find(diff(idrange_full)==max(diff(idrange_full)));
            mx = idrange_full(idrange_start)+1:idrange_full(idrange_start+1)-1;
        else
            mx = 1:length(dataRec.Data.LFP(:,chR));
        end
    
        % Initiate positive polarity of the R-peaks
        flagpeak_right = 1;

        % Normalize and find peaks for LFP_R
        LFPnorm_right = normalize(dataRec.Data.LFP(:,chR));

        % Threshold of 2 * sd (instead of 2.5 * sd as described in Stam et al. (2023)
        % https://doi.org/10.1016/j.clinph.2022.11.011) for peak detection:
        [Rpeak_right, locs_Rwave_right] = findpeaks(LFPnorm_right, 'MinPeakHeight', (2 * std(LFPnorm_right(mx))), 'MinPeakDistance', dataRec.Settings.SamplingFrequencyTimeDomain / 2);
        [Speak_right, locs_Swave_right] = findpeaks(-LFPnorm_right, 'MinPeakHeight', (2 * std(-LFPnorm_right(mx))), 'MinPeakDistance', dataRec.Settings.SamplingFrequencyTimeDomain / 2);
        if isempty(Rpeak_right)
            Rpeak_right = 0;
        end
        if isempty(Speak_right)
            Speak_right = 0;
        end

        % Flip polarity if the mean height and the number of detected R-peaks
        % are higher with negative polarity
        if mean(Speak_right) > mean(Rpeak_right) && length(Speak_right) >= length(Rpeak_right)
            locs_Rwave_right = locs_Swave_right;
            locs_Rwave_right = unique(locs_Rwave_right); locs_Rwave_right = locs_Rwave_right(find(locs_Rwave_right>0));
            LFPnorm_right = -LFPnorm_right;
            flagpeak_right = -1;
        end
    
        % Check ECG artifact consistency for the right hemisphere and decide on ECG cleaning y/n:
        if isempty(locs_Rwave_right)
            key = 'n';
        else
            % Make sure locs_Rwave is an horizontal array for manual peak selection
            if size(locs_Rwave_right,1)>size(locs_Rwave_right,2)
                locs_Rwave_right = locs_Rwave_right';
            end
            
            % Plot LFP signal and detected R-peaks
            fig = figure('Name','Detected R-peaks and PSD','NumberTitle','off', 'visible', 'off');clf();
            fig.WindowState = 'maximized';
            sgtitle([savename_json ' LFP_R Raw - Recording ' num2str(c)'], 'Interpreter', 'none')
            subplot(2,s,1:s)
            plot(time,LFPnorm_right,  'k'), hold on,plot(time(locs_Rwave_right), LFPnorm_right(locs_Rwave_right),'r*', 'markersize', 15)
            xlabel('Time [s]', 'FontSize',12)
            box off
            legend('LFPs', ['Detected R-peaks: ' num2str(length(locs_Rwave_right))], 'FontSize',11, 'location', 'northwest')
            legend('boxoff')
    
            % A few pre-requisites before the SVD ECG suppression method is run:
                % If there are no R-peaks detected at all:
            if isempty(locs_Rwave_right)
                title('Automatically detected R-peaks: no obvious R-peaks', 'Interpreter', 'none', 'FontSize', 14);
                % If heartbeats are missed: no ECG artifact is detected:
            elseif any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                title('Automatically detected R-peaks: no R-peak every three seconds', 'Interpreter', 'none', 'FontSize', 14);
                % If heart rate is too low:
            elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
                title('Automatically detected R-peaks: unreliable heart rate (< 40 bpm)', 'Interpreter', 'none', 'FontSize', 14);
            else
                title(['Automatically detected R-peaks: consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'], 'Interpreter', 'none', 'FontSize', 14);
            end
    
            for k = 1:s
                if length(LFP_R{k}) > dataRec.Settings.SamplingFrequencyTimeDomain-1
                    [P,F]=pwelch(LFP_R{k},[],[],[],250);
                    [P1s,F1s]=pwelch(LFP_R{k},hamming(1*250),[],[],250);
                    subplot(2,s,s+k)
                    plot(F,P);hold on,plot(F1s,P1s,'r','linewidth',2);ylim([0 40]),xlim([0 130]);
                    title(['Segment ' num2str(k)], 'Interpreter', 'none', 'FontSize', 14);
                    xlabel('Frequency [Hz]', 'FontSize', 12);
                    if k == 1
                        ylabel('Power', 'FontSize', 12)
                    end
                    legend('PSD large window','PSD 1s window', 'FontSize', 11)
                    box off
                else
                    subplot(2,s,s+k)
                    title(['Segment ' num2str(k) ' < 250 samples'], 'FontSize', 14')
                end
            end
    
            fig.Visible = 'on';
            savefig([savepathfig filesep [savename_json '_LFP_R_raw_rec' num2str(c) '_PSD_precleaning.fig']])
            print(gcf,'-djpeg', [savepathfig filesep [savename_json '_LFP_R_raw_rec' num2str(c) '_PSD_precleaning.jpg']]);
            
            key=input('Inspect detected R-peaks and PSDs of LFP-Right: Press "y" to apply ECG cleaning, press "n" to proceed without ECG cleaning:','s');
        end
    
        if key == 'n'
            if isempty(locs_Rwave_right)
                dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No obvious R-peaks';
                disp('ECG not applied to LFP_R: No obvious R-peaks');
            elseif any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No R-peak every three seconds';
                disp('ECG not applied to LFP_R: No R-peak every three seconds');
                close('Detected R-peaks and PSD')
            elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
                disp('ECG not applied to LFP_R: Unreliable heart rate (< 40 bpm)');
                close('Detected R-peaks and PSD')
            else
                dataRec.Data.EcgArtifactStatus{1,chR} = ['Data not ECG-suppressed despite consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'];
                disp(['ECG not applied to LFP_R despite consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)']);
                close('Detected R-peaks and PSD')
            end
            dataRec.Data.EcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.EcgArtifactStatus{3,chR} = [];
            dataRec.Data.EcgArtifactStatus{4,chR} = [];

        elseif key == 'y'
            close('Detected R-peaks and PSD')
            if any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.EcgArtifactStatus{1,chR} = 'Data ECG-suppressed despite no R-peak every three seconds';
            elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.EcgArtifactStatus{1,chR} = 'Data ECG-suppressed despite unreliable heart rate (< 40 bpm)';
            else
                dataRec.Data.EcgArtifactStatus{1,chR} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'];
            end
    
            dataRec.Data.EcgArtifactStatus{2,chR} = locs_Rwave_right;
    
            % Initiate LFP-specific settings of SVD method:
            % Initiated using 2 components, but manually select
            % specific components according to plot
            settings.Polarity               = flagpeak_right;
            settings.nComp                  = 2;
            settings.nCompFinal             = 1;

            % Location and name to save figures:
            settings.Folder                 = savepathfig;
            settings.Label                  = [savename_json '_LFP_R_raw_rec' num2str(c)];
            settings.Title                  = [savename_json ' LFP_R Raw - Recording ' num2str(c)];
    
            % Apply the SVD ECG suppression method to LFP_R and save filtered signal in LFP_Ecg:
            [dataRec.Data.LFP_Ecg(:,chR), settings.ProjOut, settings.FinalComp, locs_Rwave_right_final] = SVD_ECG_Filter(dataRec.Data.LFP(:,chR),settings,[],locs_Rwave_right,[]);

            % Save info about detected R-peaks to the status struct field
            if isempty(locs_Rwave_right_final)

                % If only one or two R-peaks are detected (very short LFP
                % signal) which are too close to the start and/or end, the 
                % (P)QRS(T) time window can not be created:
                dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: R-peaks too close to signal edges';
                dataRec.Data.EcgArtifactStatus{3,chR} = [];
                dataRec.Data.EcgArtifactStatus{4,chR} = [];
            else
                dataRec.Data.EcgArtifactStatus{3,chR} = locs_Rwave_right_final;
                dataRec.Data.EcgArtifactStatus{4,chR} = settings;
            end
            disp('ECG applied to LFP_R');
        end
    end
    
    %%
    
    % --- SVD-ECG Suppression to notch-filtered LFPs (output: LFP_LinenoiseEcg)
    
    % Generate SVD components of specific time window surrounding R-peaks
    % Here PQRST window (250 ms before and 400 ms after R-peak) is used
    % Alternatively use smaller QRS time window (200 ms before and after peaks)
    pre_R_time = rWindow(1); 
    post_R_time = rWindow(2); 
    
    % Set interactive = 1, to allow manual peak detection:
    settings.Interactive            = 1;
    settings.ArtTimeBeforePeak      = pre_R_time;
    settings.ArtTimeAfterPeak       = post_R_time;
    settings.SamplingFrequency      = dataRec.Settings.SamplingFrequencyTimeDomain;
    settings.ShowFigs               = 1;
    settings.SaveFigs               = 1;
    
    % SVD ECG suppression is applied to LFP_L and LFP_R separately
    % For left hemisphere:
    if ~isempty(chL)
    
        % ECG suppression starts with an R-peak detection algorithm
    
        % Normalize and find peaks for LFP_L
        LFPnorm_left = normalize(dataRec.Data.LFP_Linenoise(:,chL));
        LFPnorm_left = LFPnorm_left * flagpeak_left;
    
        % If SVD ECG suppression is applied to raw LFPs, use those final
        % R-peak indices for this SVD ECG suppression, otherwise use the
        % automatically detected R-peaks in raw LFP signal for ECG suppression
        % of notch-filtered LFP signal:
        if ~isempty(dataRec.Data.EcgArtifactStatus{3,chL})
            disp(['R-peak locations of previous cleaning (to raw signal) inserted as initial locations of current cleaning: ' num2str(length(dataRec.Data.EcgArtifactStatus{3,chL})) ' peaks']);
            locs_Rwave_left = dataRec.Data.EcgArtifactStatus{3,chL};
        else
            disp(['Automatic detection of R-peak locations in raw signal inserted as initial locations of current cleaning: ' num2str(length(dataRec.Data.EcgArtifactStatus{2,chL})) ' peaks']);
            locs_Rwave_left = dataRec.Data.EcgArtifactStatus{2,chL};
        end
    
        % Check ECG artifact consistency for the left hemisphere and decide on ECG cleaning y/n:
        if isempty(locs_Rwave_left)
            key = 'n';
        else
            % Make sure locs_Rwave is an horizontal array for manual peak selection
            if size(locs_Rwave_left,1)>size(locs_Rwave_left,2)
                locs_Rwave_left = locs_Rwave_left';
            end
            
            % Plot LFP signal and detected R-peaks
            fig = figure('Name','Detected R-peaks and PSD','NumberTitle','off', 'visible', 'off');clf();
            fig.WindowState = 'maximized';
            sgtitle([savename_json ' LFP_L - Linenoise filtered - Recording ' num2str(c)'], 'Interpreter', 'none')
            subplot(2,s,1:s)
            plot(time,LFPnorm_left,  'k'), hold on,plot(time(locs_Rwave_left), LFPnorm_left(locs_Rwave_left),'r*', 'markersize', 15)
            xlabel('Time [s]', 'FontSize',12)
            box off
            legend('LFPs', ['Detected R-peaks: ' num2str(length(locs_Rwave_left))], 'FontSize',11, 'location', 'northwest')
            legend('boxoff')
    
            % A few pre-requisites before the SVD ECG suppression method is run:
                % If there are no R-peaks detected at all:
            if isempty(locs_Rwave_left)
                title('Detected R-peaks: no obvious R-peaks', 'Interpreter', 'none', 'FontSize', 14);
                % If heartbeats are missed: no ECG artifact is detected:
            elseif any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                title('Detected R-peaks: no R-peak every three seconds', 'Interpreter', 'none', 'FontSize', 14);
                % If heart rate is too low:
            elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain)
                title('Detected R-peaks: unreliable heart rate (< 40 bpm)', 'Interpreter', 'none', 'FontSize', 14);
            else
                title(['Detected R-peaks: consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'], 'Interpreter', 'none', 'FontSize', 14);
            end
    
            for k = 1:s
                if length(LFP_L_linenoise{k}) > dataRec.Settings.SamplingFrequencyTimeDomain-1
                    [P,F]=pwelch(LFP_L_linenoise{k},[],[],[],250);
                    [P1s,F1s]=pwelch(LFP_L_linenoise{k},hamming(1*250),[],[],250);
                    subplot(2,s,s+k)
                    plot(F,P);hold on,plot(F1s,P1s,'r','linewidth',2);ylim([0 40]),xlim([0 130]);
                    title(['Segment ' num2str(k)], 'Interpreter', 'none', 'FontSize', 14);
                    xlabel('Frequency [Hz]', 'FontSize', 12);
                    if k == 1
                        ylabel('Power', 'FontSize', 12)
                    end
                    legend('PSD large window','PSD 1s window', 'FontSize', 11)
                    box off
                else
                    subplot(2,s,s+k)
                    title(['Segment ' num2str(k) ' < 250 samples'], 'FontSize', 14')
                end
            end
    
            fig.Visible = 'on';
            savefig([savepathfig filesep [savename_json '_LFP_L_linenoise_rec' num2str(c) '_PSD_precleaning.fig']])
            print(gcf,'-djpeg', [savepathfig filesep [savename_json '_LFP_L_linenoise_rec' num2str(c) '_PSD_precleaning.jpg']]);
    
            key=input('Inspect detected R-peaks and PSDs of LFP-Left after Notch: Press "y" to apply ECG cleaning, press "n" to proceed without ECG cleaning:','s');
        end
    
        if key == 'n'
            if isempty(locs_Rwave_left)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No obvious R-peaks';
                disp('ECG not applied to LFP_L_linenoise: No obvious R-peaks');
            elseif any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No R-peak every three seconds';
                disp('ECG not applied to LFP_L_linenoise: No R-peak every three seconds');
                close('Detected R-peaks and PSD')
            elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
                disp('ECG not applied to LFP_L_linenoise: Unreliable heart rate (< 40 bpm)');
                close('Detected R-peaks and PSD')
            else
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = ['Data not ECG-suppressed despite consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'];
                disp(['ECG not applied to LFP_L_linenoise despite consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)']);
                close('Detected R-peaks and PSD')
            end
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chL} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chL} = [];

        elseif key == 'y'
            close('Detected R-peaks and PSD')
            if any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data ECG-suppressed despite no R-peak every three seconds';
            elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data ECG-suppressed despite unreliable heart rate (< 40 bpm)';
            else
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'];
            end
    
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chL} = locs_Rwave_left;
    
            % Initiate LFP-specific settings of SVD method:
            % Initiated using 2 components, but manually select
            % specific components according to plot
            settings.Polarity               = flagpeak_left;
            settings.nComp                  = 2;
            settings.nCompFinal             = 1;

            % Location and name to save figures:
            settings.Folder                 = savepathfig;
            settings.Label                  = [savename_json '_LFP_L_linenoise_rec' num2str(c)];
            settings.Title                  = [savename_json ' LFP_L - Linenoise filtered - Recording ' num2str(c)];
    
            % If SVD ECG suppression is applied to raw LFPs, display which
            % (number of) components are chosen for ECG suppression:
            if ~isempty(dataRec.Data.EcgArtifactStatus{4,chL})
                disp(['FYI: for ECG cleaning of raw signal component(s): ' num2str(dataRec.Data.EcgArtifactStatus{4,chL}.FinalComp) ' chosen']);
            end
    
            % Apply the SVD ECG suppression method to notch-filtered LFP_L and save filtered signal in LFP_LinenoiseEcg:
            [dataRec.Data.LFP_LinenoiseEcg(:,chL), settings.ProjOut, settings.FinalComp, locs_Rwave_left_final] = SVD_ECG_Filter(dataRec.Data.LFP_Linenoise(:,chL),settings,[],locs_Rwave_left,[]);

            % Save info about detected R-peaks to the status struct field
            if isempty(locs_Rwave_left_final)

                % If only one or two R-peaks are detected (very short LFP
                % signal) which are too close to the start and/or end, the 
                % (P)QRS(T) time window can not be created:
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: R-peaks too close to signal edges';
                dataRec.Data.LinenoiseEcgArtifactStatus{3,chL} = [];
                dataRec.Data.LinenoiseEcgArtifactStatus{4,chL} = [];
            else
                dataRec.Data.LinenoiseEcgArtifactStatus{3,chL} = locs_Rwave_left_final;
                dataRec.Data.LinenoiseEcgArtifactStatus{4,chL} = settings;
            end
            disp('ECG applied to LFP_L_linenoise');
        end
    end
    
    % For right hemisphere:
    if ~isempty(chR)
    
        % ECG suppression starts with an R-peak detection algorithm
    
        % Normalize and find peaks for LFP_R
        LFPnorm_right = normalize(dataRec.Data.LFP_Linenoise(:,chR));
        LFPnorm_right = LFPnorm_right * flagpeak_right;
    
        % If SVD ECG suppression is applied to raw LFPs, use those final
        % R-peak indices for this SVD ECG suppressionm, otherwise use the
        % automatically detected R-peaks in raw LFP signal for ECG suppression
        % of notch-filtered LFP signal:
        if ~isempty(dataRec.Data.EcgArtifactStatus{3,chR})
            disp(['R-peak locations of previous cleaning (to raw signal) inserted as initial locations of current cleaning: ' num2str(length(dataRec.Data.EcgArtifactStatus{3,chR})) ' peaks']);
            locs_Rwave_right = dataRec.Data.EcgArtifactStatus{3,chR};
        else
            disp(['Automatic detection of R-peak locations in raw signal inserted as initial locations of current cleaning: ' num2str(length(dataRec.Data.EcgArtifactStatus{2,chR})) ' peaks']);
            locs_Rwave_right = dataRec.Data.EcgArtifactStatus{2,chR};
        end
    
        % Check ECG artifact consistency for the right hemisphere and decide on ECG cleaning y/n:
        if isempty(locs_Rwave_right)
            key = 'n';
        else
            % Make sure locs_Rwave is an horizontal array for manual peak selection
            if size(locs_Rwave_right,1)>size(locs_Rwave_right,2)
                locs_Rwave_right = locs_Rwave_right';
            end
            
            % Plot LFP signal and detected R-peaks
            fig = figure('Name','Detected R-peaks and PSD','NumberTitle','off', 'Visible','off');clf();
            fig.WindowState = 'maximized';
            sgtitle([savename_json ' LFP_R - Linenoise filtered - Recording ' num2str(c)'], 'Interpreter', 'none')
            subplot(2,s,1:s)
            plot(time,LFPnorm_right,  'k'), hold on,plot(time(locs_Rwave_right), LFPnorm_right(locs_Rwave_right),'r*', 'markersize', 15)
            xlabel('Time [s]', 'FontSize',12)
            box off
            legend('LFPs', ['Detected R-peaks: ' num2str(length(locs_Rwave_right))], 'FontSize',11, 'location', 'northwest')
            legend('boxoff')
    
            % A few pre-requisites before the SVD ECG suppression method is run:
                % If there are no R-peaks detected at all:
            if isempty(locs_Rwave_right)
                title('Detected R-peaks: no obvious R-peaks', 'Interpreter', 'none', 'FontSize', 14);
                % If heartbeats are missed: no ECG artifact is detected:
            elseif any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                title('Detected R-peaks: no R-peak every three seconds', 'Interpreter', 'none', 'FontSize', 14);
                % If heart rate is too low:
            elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
                title('Detected R-peaks: unreliable heart rate (< 40 bpm)', 'Interpreter', 'none', 'FontSize', 14);
            else
                title(['Detected R-peaks: consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'], 'Interpreter', 'none', 'FontSize', 14);
            end
    
            for k = 1:s
                if length(LFP_R_linenoise{k}) > dataRec.Settings.SamplingFrequencyTimeDomain-1
                    [P,F]=pwelch(LFP_R_linenoise{k},[],[],[],250);
                    [P1s,F1s]=pwelch(LFP_R_linenoise{k},hamming(1*250),[],[],250);
                    subplot(2,s,s+k)
                    plot(F,P);hold on,plot(F1s,P1s,'r','linewidth',2);ylim([0 40]),xlim([0 130]);
                    title(['Segment ' num2str(k)], 'Interpreter', 'none', 'FontSize', 14);
                    xlabel('Frequency [Hz]', 'FontSize', 12); legend('PSD large window','PSD 1s window', 'FontSize', 11)
                    if k == 1
                        ylabel('Power', 'FontSize', 12)
                    end
                    box off
                else
                    subplot(2,s,s+k)
                    title(['Segment ' num2str(k) ' < 250 samples'], 'FontSize', 14')
                end
            end
    
            fig.Visible = 'on';
            savefig([savepathfig filesep [savename_json '_LFP_R_linenoise_rec' num2str(c) '_PSD_precleaning.fig']])
            print(gcf,'-djpeg', [savepathfig filesep [savename_json '_LFP_R_linenoise_rec' num2str(c) '_PSD_precleaning.jpg']]);
    
            key=input('Inspect detected R-peaks and PSDs of LFP-Right after Notch: Press "y" to apply ECG cleaning, press "n" to proceed without ECG cleaning:','s');
    
        end
    
        if key == 'n'
            if isempty(locs_Rwave_right)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No obvious R-peaks';
                disp('ECG not applied to LFP_R_linenoise: No obvious R-peaks');
            elseif any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No R-peak every three seconds';
                disp('ECG not applied to LFP_R_linenoise: No R-peak every three seconds');
                close('Detected R-peaks and PSD')
            elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
                disp('ECG not applied to LFP_R_linenoise: Unreliable heart rate (< 40 bpm)');
                close('Detected R-peaks and PSD')
            else
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = ['Data not ECG-suppressed despite consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'];
                disp(['ECG not applied to LFP_R_linenoise despite consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)']);
                close('Detected R-peaks and PSD')
            end
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chR} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chR} = [];

        elseif key == 'y'
            close('Detected R-peaks and PSD')
            if any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data ECG-suppressed despite no R-peak every three seconds';
            elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data ECG-suppressed despite unreliable heart rate (< 40 bpm)';
            else
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'];
            end
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chR} = locs_Rwave_right;
    
            % Initiate LFP-specific settings of SVD method:
            % Initiated using 2 components, but manually select
            % specific components according to plot
            settings.Polarity               = flagpeak_right;
            settings.nComp                  = 2;
            settings.nCompFinal             = 1;

            % Location and name to save figures:
            settings.Folder                 = savepathfig;
            settings.Label                  = [savename_json '_LFP_R_Linenoise_rec' num2str(c)];
            settings.Title                  = [savename_json ' LFP_R - Linenoise filtered - Recording ' num2str(c)];
    
            % If SVD ECG suppression is applied to raw LFPs, display which
            % (number of) components are chosen for ECG suppression:
            if ~isempty(dataRec.Data.EcgArtifactStatus{4,chR})
                disp(['FYI: for ECG cleaning of raw signal component(s): ' num2str(dataRec.Data.EcgArtifactStatus{4,chR}.finalcomp) ' chosen']);
            end
            
            % Apply the SVD ECG suppression method to notch-filtered LFP_R and save filtered signal in LFP_LinenoiseEcg:
            [dataRec.Data.LFP_LinenoiseEcg(:,chR), settings.ProjOut, settings.FinalComp, locs_Rwave_right_final] = SVD_ECG_Filter(dataRec.Data.LFP_Linenoise(:,chR),settings,[],locs_Rwave_right,[]);

            % Save info about detected R-peaks to the status struct field
            if isempty(locs_Rwave_right_final)

                % If only one or two R-peaks are detected (very short LFP
                % signal) which are too close to the start and/or end, the 
                % (P)QRS(T) time window can not be created:
                dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: R-peaks too close to signal edges';
                dataRec.Data.LinenoiseEcgArtifactStatus{3,chR} = [];
                dataRec.Data.LinenoiseEcgArtifactStatus{4,chR} = [];
            else
                dataRec.Data.LinenoiseEcgArtifactStatus{3,chR} = locs_Rwave_right_final;
                dataRec.Data.LinenoiseEcgArtifactStatus{4,chR} = settings;
            end
            disp('ECG applied to LFP_R_linenoise');
        end
    end
    
    disp('Inspect ECG cleaning result. Press any key to continue.');
    pause
    close all
    
end

