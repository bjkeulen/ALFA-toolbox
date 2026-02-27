% Author: B.J. Keulen
% Date: 25-10-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for plotting data from dataSetup, the output from the getSetup
% function of the Observe toolbox. A subplot is created for each
% hemisphere. The PSDs are calculated from the raw LFP signal using
% a Hamming window of 256 data points with 50% overlap, resulting in a
% frequency resolution of ~0.98 Hz which is also used in the Percept.
% 
% INPUT
%   dataSetup       =   struct containing BrainSense Setup data
%   savepath        =   string with path into which to save figure
%   savename        =   string of name under which to save figure
%   showFig         =   int to determine whether to show the figure or not
%
% OUTPUT
%   none

function plotSetup(dataSetup, savepath, savename, showFig)

    arguments
        dataSetup
        savepath
        savename
        showFig = 0
    end


    % Define colors and linestyles
    colors = orderedcolors("gem");

    % Loop over number of runs
    for r = 1:dataSetup.nRuns

        % Create figure
        fig = figure('Name','Setup','visible','off','WindowState','maximized');
        fig.WindowState = 'maximized';
        sgtitle({dataSetup.Info.OriginalFile(1:end-5); ['BrainSense Setup - run ' num2str(r) '/' num2str(dataSetup.nRuns)]},'Interpreter','none')
        ymax = 0;
    
        % Loop over hemispheres and channels, plot data
        sides = fieldnames(dataSetup.Data);
        for s = 1:length(sides)
            
            % Create subplot
            subplot(1,length(sides),s)

            % Find indices of channels from current run and loop over 
            labels = {};
            idx = find([dataSetup.Data.(sides{s}).Run] == r);
            for c = 1:length(idx)

                % Retrieve data, calculate PSD
                data = dataSetup.Data.(sides{s})(idx(c));
                [psd,freq] = pwelch(data.LFP, 256, [], 0:1/(256/dataSetup.SamplingFrequency):100, dataSetup.SamplingFrequency);
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
            title([sides{s} ' electrode'],'FontSize', 14)
            legend(labels,'Interpreter','none','FontSize', 10)
            legend('boxoff')

        end
    
        % Set y limit
        ax = findall(fig,'type','axes');
        ylim(ax, [0 ceil(ymax)])
        linkaxes(ax)

        % Create directory if needed
        if ~exist(savepath, 'dir')
            mkdir(savepath)
        end
    
        % Save figure and keep visible if requested
        fig.Visible = 'on';
        saveas(fig, [savepath filesep savename '_Setup_run' num2str(r) '.fig'])
        frame = getframe(fig);
        imwrite(frame.cdata, [savepath filesep savename '_Setup_run' num2str(r) '.png'])
        if  ~showFig
            close(fig)
        end
    end
end
