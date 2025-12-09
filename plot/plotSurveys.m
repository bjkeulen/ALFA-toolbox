% Author: B.J. Keulen
% Date: 25-10-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function for plotting data from dataSurvey and dataIdentifier, the output
% from the getSurvey and getIdentifier functions of the Observe toolbox. 
% Four different subplots are created: for each combination of hemisphere 
% and ring/segments. The PSDs are calculated from the raw LFP signal using
% a Hamming window of 256 data points with 50% overlap, resulting in a
% frequency resolution of ~0.98 Hz which is also used in the Percept.
% 
% INPUT
%   dataSurveys     =   struct containing BrainSense Electrode Survey or
%                       Electrode Identifier data
%   savepath        =   string with path into which to save figure
%   savename        =   string of name under which to save figure
%   showFig         =   int to determine whether to show the figure or not
%
% OUTPUT
%   none

function plotSurveys(dataSurveys, savepath, savename, showFig)

    arguments
        dataSurveys
        savepath
        savename
        showFig = 0
    end

    % Define colors
    colors = [orderedcolors("gem"); orderedcolors("glow")];

    % Loop over number of runs
    if ~isfield(dataSurveys, 'nRuns')
        dataSurveys.nRuns = 1;
    end
    for r = 1:dataSurveys.nRuns

        % Create figure
        if strcmp(dataSurveys.Mode, 'ElectrodeSurvey')
            fig = figure('Name','Electrode Survey','visible','off','WindowState','maximized');
            sgtitle({dataSurveys.Info.OriginalFile(1:end-5); ['BrainSense Electrode Survey - run ' num2str(r) '/' num2str(dataSurveys.nRuns)]}, ...
                    'Interpreter','none')
        else
            fig = figure('Name','Electrode Identifier','visible','off','WindowState','maximized');
            sgtitle({dataSurveys.Info.OriginalFile(1:end-5); ['BrainSense Electrode Identifier - run ' num2str(r) '/' num2str(dataSurveys.nRuns)]}, ...
                'Interpreter','none')
        end
        isub = [1 3; 2 4];
        ymax = 0;
    
        % Loop over hemispheres
        sides = fieldnames(dataSurveys.Data);
        for s = 1:length(sides)
    
            % Loop over modes
            modes = fieldnames(dataSurveys.Data.(sides{s}));
            for m = 1:length(modes)
    
                % Determine index of subplot
                if length(sides) == 1
                    iplot = m;
                elseif length(modes) == 1
                    iplot = s;
                else
                    iplot = isub(s,m);
                end
    
                % Create subfigure
                subplot(length(modes),length(sides),iplot)

                % Find indices of channels from current run and loop over 
                labels = {};
                idx = find([dataSurveys.Data.(sides{s}).(modes{m}).Run] == r);
                for c = 1:length(idx)

                    % Retrieve data, calculate PSD
                    data = dataSurveys.Data.(sides{s}).(modes{m})(idx(c));
                    [psd,freq] = pwelch(data.LFP, 256, [], 0:1/(256/dataSurveys.SamplingFrequency):100, dataSurveys.SamplingFrequency);
                    plot(round(freq,2), psd,'LineWidth',1,'Color',colors(c,:)); hold on;

                    % Add channel to labels
                    labels = [labels, {data.Channel}];

                    % Check for max y-value
                    if max(psd) > ymax
                        ymax = max(psd);
                    end
                end

                % Settings
                box off
                xlabel('Frequency [Hz]','FontSize', 12)
                ylabel('PSD [\muV^{2}/Hz]','FontSize', 12)
                title(sprintf('%s %ss', sides{m}, lower(modes{m})),'FontSize', 14)
                legend(labels,'Interpreter','none','FontSize', 10)
                legend('boxoff')

            end
        end
    
        % Set y-limit
        ax = findall(fig,'type','axes');
        ylim(ax, [0 ceil(ymax)])
        linkaxes(ax)

        % Create directory if needed
        if ~exist(savepath, 'dir')
            mkdir(savepath)
        end
    
        % Save figure and keep visible if requested
        if strcmp(dataSurveys.Mode, 'ElectrodeSurvey')
            savename_fig = [savename '_Survey_run' num2str(r)];
        else
            savename_fig = [savename '_Identifier_run' num2str(r)];
        end
        fig.Visible = 'on';
        saveas(fig, [savepath filesep savename_fig])
        frame = getframe(fig);
        imwrite(frame.cdata, [savepath filesep savename_fig '.png'])
        if ~showFig
            close(fig)
        end
    end
end