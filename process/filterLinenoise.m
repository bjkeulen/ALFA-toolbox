% Author: M.J. Stam
% Date: 27-08-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to apply a power line noise filter to Streaming recordings. The
% central frequency and bandwidth of the filter are given by the input
% variable linenoise.
%
% INPUT
%   dataRec         =   struct containing info and BrainSenseTimeDomain
%                       data from one recording
%                       Streaming time array, and LFP data in raw form
%   linenoise       =   Frequency and bandwidth [Hz] of line-noise to 
%                       remove
%
% OUTPUT
%   dataRec         =   struct containing info and BrainSenseTimeDomain
%                       Streaming time array, and LFP data in raw and  
%                       notch filtered form

function dataRec = filterLinenoise(dataRec, linenoise)

% --- Linenoise notch filter to raw LFPs (output: lfp_notch)

% Cannot be applied to NaNs, so this finds the start and end indices of 
% segments of the Streaming data without missing values
i_seg = [1; find(diff(isnan(dataRec.Data.LFP(:,1))))+1; length(dataRec.Data.LFP)];

% Initiate such that if present, NaNs are correctly copied:
dataRec.Data.LFP_Linenoise = dataRec.Data.LFP;

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

    % Check if recording contains data for left and right hemisphere:
    chL=find(contains(dataRec.Data.Channel,'LEFT')==1);
    if ~isempty(chL)
        LFP_L=dataRec.Data.LFP(ix,chL);
    end

    chR=find(contains(dataRec.Data.Channel,'RIGHT')==1);
    if ~isempty(chR)
        LFP_R=dataRec.Data.LFP(ix,chR);
    end

    % Set parameters for notch filter
    NotchFreq = linenoise(1);
    NotchWidth = linenoise(2);

    if NotchWidth
        % Apply notch to Streaming segment and save in struct field: lfp_notch
        for j=1:length(NotchFreq)
            [B,A]=butter(4,[NotchFreq(j)-NotchWidth(j) NotchFreq(j)+NotchWidth(j)]/(.5*dataRec.Settings.SamplingFrequencyTimeDomain),'stop');
            if ~isempty(chL)
                LFP_L_notch=filtfilt(B,A,LFP_L);
                dataRec.Data.LFP_Linenoise(ix,chL) = LFP_L_notch;
                % disp(['Notch applied to segment ' num2str(s) ' of LFP_L']);
            end

            if ~isempty(chR)
                LFP_R_notch=filtfilt(B,A,LFP_R);
                dataRec.Data.LFP_Linenoise(ix,chR) = LFP_R_notch;
                % disp(['Notch applied to segment ' num2str(s) ' of LFP_R']);
            end
        end
    end
end
