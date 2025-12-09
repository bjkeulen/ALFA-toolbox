% Author: B.J. Keulen
% Date: 19-09-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function for retrieving stimulation and sensing parameters for a given
% stimulation group.
% 
% INPUT
%   group           =   struct containing group info, stored within a
%                       JSON report under Groups.Initial or Groups.Final
%
% OUTPUT
%   rate            =   1x2 double containing the stimulation frequency
%   pulseWidth      =   1x2 double containing the pulse width
%   chan            =   1x2 cell array containing the sensing channel
%   freq            =   1x2 double containing the sensing frequency
%   ampLower        =   1x2 double containing the lower stimulation
%                       amplitude
%   ampUpper        =   1x2 double containing the upper stimulation
%                       amplitude
%   lfpLower        =   1x2 double containing the lower LFP threshold
%   lfpUpper        =   1x2 double containing the upper LFP threshold 

function [rate, pulsewidth, chan, freq, ampLower, ampUpper, lfpLower, lfpUpper] = getGroupInfo(group)

    rate = [nan nan];
    pulsewidth = [nan nan];
    chan = cell(1,2);
    freq = [nan nan];
    ampLower = [nan nan];
    ampUpper = [nan nan];
    lfpLower = [nan nan];
    lfpUpper = [nan nan];

    if isfield(group, 'ProgramSettings')
        if isfield(group.ProgramSettings, 'SensingChannel')
            for s = 1:length(group.ProgramSettings.SensingChannel)
    
                % Check type of field SensingChannel
                try
                    sensing = group.ProgramSettings.SensingChannel{s};
                catch
                    sensing = group.ProgramSettings.SensingChannel(s);
                end

                % Get left frequency and channel
                if contains(sensing.HemisphereLocation, 'Left')
                    freq(1) = sensing.SensingSetup.FrequencyInHertz;
                    pulsewidth(1) = sensing.PulseWidthInMicroSecond;
                    rate(1) = sensing.RateInHertz;

                    chanL = split(sensing.Channel,'.');
                    if isfield(sensing, 'GangedToHemisphere')
                        if contains(sensing.GangedToHemisphere, 'Right')
                            chan{1} = strcat('RIGHT_',chanL{end});
                        else
                            chan{1} = strcat('LEFT_',chanL{end});
                        end
                    else
                        chan{1} = strcat('LEFT_',chanL{end});
                    end

                    % Check if adaptive therapy, collect parameters if so
                    if isfield(sensing, 'AdaptiveTherapyStatus')
                        if ~contains(sensing.AdaptiveTherapyStatus, 'NOT_CONFIGURED')
                            lfpLower(1) = sensing.LowerLfpThreshold;
                            lfpUpper(1) = sensing.UpperLfpThreshold;
                            if isfield(sensing, 'LowerLimitInMilliAmps')
                                ampLower(1) = sensing.LowerLimitInMilliAmps;
                                ampUpper(1) = sensing.UpperLimitInMilliAmps;
                            end
                        end
                    end                          

                % Get right frequency and channel
                elseif contains(sensing.HemisphereLocation, 'Right')
                    freq(2) = sensing.SensingSetup.FrequencyInHertz;
                    pulsewidth(2) = sensing.PulseWidthInMicroSecond;
                    rate(2) = sensing.RateInHertz;

                    chanR = split(sensing.Channel,'.');
                    if isfield(sensing, 'GangedToHemisphere')
                        if contains(sensing.GangedToHemisphere, 'Left')
                            chan{2} = strcat('LEFT_',chanR{end});
                        else
                            chan{2} = strcat('RIGHT_',chanR{end});
                        end
                    else
                        chan{2} = strcat('RIGHT_',chanR{end});
                    end

                    % Check if adaptive therapy, collect parameters if so
                    if isfield(sensing, 'AdaptiveTherapyStatus')
                        if ~contains(sensing.AdaptiveTherapyStatus, 'NOT_CONFIGURED')
                            lfpLower(2) = sensing.LowerLfpThreshold;
                            lfpUpper(2) = sensing.UpperLfpThreshold;
                            if isfield(sensing, 'LowerLimitInMilliAmps')
                                ampLower(2) = sensing.LowerLimitInMilliAmps;
                                ampUpper(2) = sensing.UpperLimitInMilliAmps;
                            end
                        end
                    end       
                end
            end
        else

            % Get left stimulation parameters
            if isfield(group.ProgramSettings,'LeftHemisphere')
                if isfield(group.ProgramSettings.LeftHemisphere, 'Programs')
                    pulsewidth(1) = group.ProgramSettings.LeftHemisphere.Programs.PulseWidthInMicroSecond;
                    if isfield(group.ProgramSettings, 'RateInHertz')
                        rate(1) = group.ProgramSettings.RateInHertz; 
                    elseif isfield(group.ProgramSettings.LeftHemisphere.Programs, 'RateInHertz')
                        rate(1) = group.ProgramSettings.LeftHemisphere.Programs.RateInHertz; 
                    end
                    if isfield(group.ProgramSettings.LeftHemisphere.Programs, 'LowerLimitInMilliAmps')
                        ampLower(1) = group.ProgramSettings.LeftHemisphere.Programs.LowerLimitInMilliAmps;
                        ampUpper(1) = group.ProgramSettings.LeftHemisphere.Programs.UpperLimitInMilliAmps;
                    end
                end
            end
    
            % Get right stimulation parameters
            if isfield(group.ProgramSettings,'RightHemisphere')
                if isfield(group.ProgramSettings.RightHemisphere, 'Programs')
                    pulsewidth(2) = group.ProgramSettings.RightHemisphere.Programs.PulseWidthInMicroSecond;
                    if isfield(group.ProgramSettings, 'RateInHertz')
                        rate(2) = group.ProgramSettings.RateInHertz; 
                    elseif isfield(group.ProgramSettings.RightHemisphere.Programs, 'RateInHertz')
                        rate(2) = group.ProgramSettings.RightHemisphere.Programs.RateInHertz; 
                    end
                    if isfield(group.ProgramSettings.RightHemisphere.Programs, 'LowerLimitInMilliAmps')
                        ampLower(2) = group.ProgramSettings.RightHemisphere.Programs.LowerLimitInMilliAmps;
                        ampUpper(2) = group.ProgramSettings.RightHemisphere.Programs.UpperLimitInMilliAmps;
                    end
                end   
            end
        end
    end
end