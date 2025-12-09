% Author: M.J. Stam
% Date: 27-08-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to automatically apply the Singular Value Decomposition ECG 
% suppression method by B.C.M. van Wijk, as described in Stam et al. 2023 
% (see reference below) to the BrainSenseTimeDomain Streaming data.
%
% Stam MJ, van Wijk BCM, Sharma P, Beudel M, Pina-Fuentes DA, de Bie RMA, 
% Schuurman PR, Neumann W-J, Buijink AWG. A comparison of methods to 
% suppress electrocardiographic artifacts in local field potential 
% recordings. Clin Neurophysiol. 2023 Feb:146:147-161. DOI: 
% 10.1016/j.clinph.2022.11.011
%
% INPUT
%   dataRec         =   struct containing info and BrainSenseTimeDomain
%                       Streaming time array, and LFP data in raw form of
%                       one recording.
%   rWindow         =   1x2 array defining window around R-peak [before 
%                       after] for calculating SVD components, in seconds.
%
% OUTPUT
%   dataRec         =   struct containing info and BrainSenseTimeDomain
%                       Streaming time array, and LFP data in raw, 
%                       ECG-suppressed and notch + ECG-suppressed form                       

function dataRec = filterECG(dataRec, rWindow)

    % NaNs if ECG is not suppressed:
    dataRec.Data.LFP_Ecg = NaN(size(dataRec.Data.LFP,1), size(dataRec.Data.LFP,2));
    dataRec.Data.LFP_LinenoiseEcg = NaN(size(dataRec.Data.LFP,1), size(dataRec.Data.LFP,2));
    
    % Check if recording contains data for left and right hemisphere:
    chL=find(contains(dataRec.Data.Channel,'LEFT')==1);
    chR=find(contains(dataRec.Data.Channel,'RIGHT')==1);
    
    %%
    
    % Initiate the general settings of SVD method
    
    % Generate SVD components of specific time window surrounding R-peaks
    % Default PQRST window is 250 ms before and 400 ms after R-peak
    % Alternatively use smaller QRS time window (200 ms before and after peaks)
    pre_R_time = rWindow(1);
    post_R_time = rWindow(2);
    
    settings.Interactive            = 0;
    settings.ArtTimeBeforePeak      = pre_R_time;
    settings.ArtTimeAfterPeak       = post_R_time;
    settings.SamplingFrequency      = dataRec.Settings.SamplingFrequencyTimeDomain;
    settings.ShowFigs               = 0;
    settings.SaveFigs               = 0;
    
    % SVD ECG suppression is applied to Left and Right LFPs separately
    
    % For left hemisphere:
    if ~isempty(chL)
    
        % ECG suppression starts with an R-peak detection algorithm
        % disp('Performing R-peak detection for LFP_L');
    
        % Set threshold for minimal ECG R-peak height
        % Threshold is set according to the part of the signal before a 
        % potential artifact (> 5 std deviation) occurs (mx):
        upthres = mean(dataRec.Data.LFP(:,chL), "omitmissing")+(5*std(dataRec.Data.LFP(:,chL), "omitmissing"));
        lowthres = mean(dataRec.Data.LFP(:,chL), "omitmissing")-(5*std(dataRec.Data.LFP(:,chL), "omitmissing"));
    
        % Convert to string of zeros and ones
        Y = dataRec.Data.LFP(:,chL)>lowthres & dataRec.Data.LFP(:,chL)<upthres; 
        idx=find(Y<1);
        if ~isempty(find(Y<1,1))
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
    
        % A few pre-requisites before the SVD ECG suppression method is run:
           % If there are no R-peaks detected at all:
        if isempty(locs_Rwave_left)
            dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No obvious R-peaks';
            dataRec.Data.EcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.EcgArtifactStatus{3,chL} = [];
            dataRec.Data.EcgArtifactStatus{4,chL} = [];
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No obvious R-peaks';
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chL} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chL} = [];
            % If heartbeats are missed: no ECG artifact is detected:
        elseif any(diff(locs_Rwave_left) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
            dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No R-peak every three seconds';
            dataRec.Data.EcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.EcgArtifactStatus{3,chL} = [];
            dataRec.Data.EcgArtifactStatus{4,chL} = [];
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: No R-peak every three seconds';
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chL} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chL} = [];
            % If heart rate is too low:
        elseif length(locs_Rwave_left) < (40 / 60) * (length(LFPnorm_left) / dataRec.Settings.SamplingFrequencyTimeDomain) 
            dataRec.Data.EcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
            dataRec.Data.EcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.EcgArtifactStatus{3,chL} = [];
            dataRec.Data.EcgArtifactStatus{4,chL} = [];
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chL} = locs_Rwave_left;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chL} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chL} = [];
        else
            dataRec.Data.EcgArtifactStatus{1,chL} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'];
            dataRec.Data.EcgArtifactStatus{2,chL} = locs_Rwave_left;
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chL} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_left)) ' peaks detected)'];
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chL} = locs_Rwave_left;
    
            % Initiate LFP-specific settings of SVD method:
            % For automatic SVD: use components 1 and 2, do not manually select
            % specific components
            settings.Polarity               = flagpeak_left;
            settings.nComp                  = 2; 
            settings.nCompFinal             = 0;
    
            % Apply the SVD ECG suppression method to LFP_L and save filtered signal in LFP_Ecg:
            disp('Apply the SVD ECG suppression method to LFP_L and save filtered signal in LFP_Ecg');
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
    
            clear locs_Rwave_left_final;
    
            % Apply the SVD ECG Suppression to LFP_L_linenoise and save filtered signal in LFP_LinenoiseEcg:
            disp('Apply the SVD ECG suppression method to LFP_L_linenoise and save filtered signal in LFP_LinenoiseEcg');
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
    
        end
    end
    
    % For right hemisphere:
    if ~isempty(chR)
    
        % ECG suppression starts with an R-peak detection algorithm
        % disp('Performing R-peak detection for LFP_R');
    
        % Set threshold for minimal ECG R-peak height
        % Threshold is set according to the part of the signal before a 
        % potential artifact (> 5 std deviation) occurs (mx):
        upthres = mean(dataRec.Data.LFP(:,chR), "omitmissing")+(5*std(dataRec.Data.LFP(:,chR), "omitmissing"));
        lowthres = mean(dataRec.Data.LFP(:,chR), "omitmissing")-(5*std(dataRec.Data.LFP(:,chR), "omitmissing"));
    
        % Convert to string of zeros and ones
        Y = dataRec.Data.LFP(:,chR)>lowthres&dataRec.Data.LFP(:,chR)<upthres; 
        idx=find(Y<1);
        if ~isempty(find(Y<1,1))
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
    
        % A few pre-requisites before the SVD ECG suppression method is run:
           % If there are no R-peaks detected at all:
        if isempty(locs_Rwave_right)
            dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No obvious R-peaks';
            dataRec.Data.EcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.EcgArtifactStatus{3,chR} = [];
            dataRec.Data.EcgArtifactStatus{4,chR} = [];
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No obvious R-peaks';
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chR} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chR} = [];
            % If heartbeats are missed: no ECG artifact is detected:
        elseif any(diff(locs_Rwave_right) / dataRec.Settings.SamplingFrequencyTimeDomain > 3)
            dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No R-peak every three seconds';
            dataRec.Data.EcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.EcgArtifactStatus{3,chR} = [];
            dataRec.Data.EcgArtifactStatus{4,chR} = [];
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: No R-peak every three seconds';
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chR} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chR} = [];
            % If heart rate is too low:
        elseif length(locs_Rwave_right) < (40 / 60) * (length(LFPnorm_right) / dataRec.Settings.SamplingFrequencyTimeDomain)
            dataRec.Data.EcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
            dataRec.Data.EcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.EcgArtifactStatus{3,chR} = [];
            dataRec.Data.EcgArtifactStatus{4,chR} = [];
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = 'Data not ECG-suppressed: Unreliable heart rate (< 40 bpm)';
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chR} = locs_Rwave_right;
            dataRec.Data.LinenoiseEcgArtifactStatus{3,chR} = [];
            dataRec.Data.LinenoiseEcgArtifactStatus{4,chR} = [];
        else
            dataRec.Data.EcgArtifactStatus{1,chR} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'];
            dataRec.Data.EcgArtifactStatus{2,chR} = locs_Rwave_right;
    
            dataRec.Data.LinenoiseEcgArtifactStatus{1,chR} = ['Data ECG-suppressed: Consistent ECG artifact (' num2str(length(locs_Rwave_right)) ' peaks detected)'];
            dataRec.Data.LinenoiseEcgArtifactStatus{2,chR} = locs_Rwave_right;
    
            % Initiate LFP-specific settings of SVD method:
            % For automatic SVD: use components 1 and 2, do not manually select
            % specific components
            settings.Polarity               = flagpeak_right;
            settings.nComp                  = 2; 
            settings.nCompFinal             = 0;
    
            % Apply the SVD ECG suppression method to LFP_R and save filtered signal in LFP_Ecg:
            disp('Apply the SVD ECG suppression method to LFP_R and save filtered signal in LFP_Ecg');
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
    
            clear locs_Rwave_right_final; 
    
            % Apply the SVD ECG suppression method to LFP_R_linenoise and save filtered signal in lfp_LinenoiseEcg:
            disp('Apply the SVD ECG suppression method to LFP_R_linenoise and save filtered signal in LFP_LinenoiseEcg');
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
        end
    end
end

