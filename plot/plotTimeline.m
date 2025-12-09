% Author: B.J. Keulen
% Date: 25-10-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function for plotting data from dataTimeline, the output from the
% getTimeline function of the Observe toolbox. The timeline data, 
% stimulation amplitude and if applicable, LFP thresholds and stimulation
% limits are plotted. Sensing channel and sensing frequency are annotated 
% with a green vertical line at the beginning of the dataset and whenever 
% one of the parameters has changed. Sessions are shown as light green
% patches.
% 
% INPUT
%   dataTimeline    =   struct containing BrainSense Timeline data
%   savepath        =   string with path into which to save figure
%   savename        =   string of name under which to save figure
%   showFig         =   int to determine whether to show the figure or not
%
% OUTPUT
%   none

function plotTimeline(dataTimeline, savepath, savename, showFig)

    arguments
        dataTimeline
        savepath
        savename
        showFig = 0
    end

    % Define colors
    colors = orderedcolors("gem");

    % Get time, sensing channel and sensing frequency data outside sessions
    iS = nan([length(dataTimeline.Data.DateTime),height(dataTimeline.History.Sessions)]);
    for s = 1:height(dataTimeline.History.Sessions)
        iS(:,s) = isbetween(dataTimeline.Data.DateTime, dataTimeline.History.Sessions{s,"SessionStart"}, dataTimeline.History.Sessions{s,"SessionEnd"});
    end
    time = dataTimeline.Data.DateTime(~any(iS,2));
    freq = dataTimeline.Data.SensingFrequency(~any(iS,2),:);
    freq(isnan(freq)) = 0;
    chan = dataTimeline.Data.SensingChannel(~any(iS,2),:);
    chan(cellfun(@isempty,chan)) = {char()};

    % Find changes in sensing frequency or channel
    idxL = unique([1; find(diff(freq(:,1)))+1; find(~strcmp([chan(1,1); chan(1:end-1,1)], chan(:,1)))]);
    idxR = unique([1; find(diff(freq(:,2)))+1; find(~strcmp([chan(1,2); chan(1:end-1,2)], chan(:,2)))]);

    % Define labels for changes in frequency or channel on left lead
    textL = strings(length(idxL),1);
    for i = 1:length(idxL)
        fi = freq(idxL(i),1);
        ci = chan{idxL(i),1};

        if ~isempty(ci) && fi ~= 0
            textL(i) = "\color{black}" + sprintf("%s | %.2f Hz",ci,fi);
        elseif strcmp(dataTimeline.History.Sessions{find(time(idxL(i)) < dataTimeline.History.Sessions{:,'SessionStart'},1,'first'),'SessionType'},'Missing Session')
            textL(i) = "\color{black}" + sprintf("Info missing");
        else
            textL(i) = "\color{black}" + sprintf("No sensing");
        end
    end
    textL = replace(textL,{'GROUP_','_AND_','ZERO','ONE','TWO','THREE','LEFT_','RIGHT_'},{'','','0','1','2','3','L-','R-'});

    % Define labels for changes in frequency or channel on right lead
    textR = strings(length(idxR),1);
    for i = 1:length(idxR)
        fi = freq(idxR(i),2);
        ci = chan{idxR(i),2};

        if ~isempty(ci) && fi ~= 0
            textR(i) = "\color{black}" + sprintf("%s | %.2f Hz",ci,fi);
        elseif strcmp(dataTimeline.History.Sessions{find(time(idxR(i)) < dataTimeline.History.Sessions{:,'SessionStart'},1,'first'),'SessionType'},'Missing Session')
            textR(i) = "\color{black}" + sprintf("Info missing");
        else
            textR(i) = "\color{black}" + sprintf("No sensing");
        end
    end
    textR = replace(textR,{'GROUP_','_AND_','ZERO','ONE','TWO','THREE','LEFT_','RIGHT_'},{'','','0','1','2','3','L-','R-'});

    % Define maximal y-value while ignoring outliers
    mask = abs(dataTimeline.Data.LfpPower) <= mean(dataTimeline.Data.LfpPower,'omitnan') + 5*std(dataTimeline.Data.LfpPower,'omitnan');
    ymax = [1.2*max(dataTimeline.Data.StimulationAmplitude(mask),[],'all') 1.1*max(dataTimeline.Data.LfpPower(mask),[],'all')];
    if ymax(1) < 1
        ymax(1) = 1;
    end

    % Create figure, plot data
    fig = figure('Name','Timeline','visible','off','WindowState','maximized'); hold on
    sgtitle([savename ' - BrainSense Timeline'],'Interpreter','none')
    titles = {'Left electrode', 'Right electrode'};
    for i = 1:2

        % Create subplot
        subplot(2,1,i); hold on
        set(gca,'color', [0.9 0.9 0.9]);
        title(titles{i},'FontSize', 14);
        xlabel('DateTime','FontSize', 12)

        % Plot stimulation amplitudes
        yyaxis right; 
        plot(dataTimeline.Data.DateTime, dataTimeline.Data.StimulationAmplitude(:,i),'Color',colors(2,:))
        plot(dataTimeline.Data.DateTime, dataTimeline.Data.LowerStimulationLimit(:,i),'LineStyle',':','LineWidth',1,'Color',colors(2,:))
        plot(dataTimeline.Data.DateTime, dataTimeline.Data.UpperStimulationLimit(:,i),'LineStyle',':','LineWidth',1,'Color',colors(2,:))
        ylabel('Stimulation amplitude [mA]','FontSize',12,'Rotation',270)
        ylim([-0.01 ymax(1)]);

        % Plot LFP values
        yyaxis left;
        plot(dataTimeline.Data.DateTime, dataTimeline.Data.LfpPower(:,i),'Color',colors(1,:))
        plot(dataTimeline.Data.DateTime, dataTimeline.Data.LowerLfpThreshold(:,i),'LineStyle','--','LineWidth',1.5,'Color',colors(3,:))
        plot(dataTimeline.Data.DateTime, dataTimeline.Data.UpperLfpThreshold(:,i),'LineStyle','--','LineWidth',1.5,'Color',colors(3,:))
        ylabel('LFP power','FontSize',12)
        ylim([0 ymax(2)]);

        % Plot changes in channel or frequency
        if i == 1
            xline(time(idxL), '-', textL, 'LineWidth', 1.5, 'FontSize', 8, 'Color', colors(5,:), ...
                  'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'center')
        else
            xline(time(idxR), '-', textR, 'LineWidth', 1.5, 'FontSize', 8, 'Color', colors(5,:), ...
                  'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'center')
        end
        fill([dataTimeline.History.Sessions{:,2:3} dataTimeline.History.Sessions{:,[3 2]}]', ...
             [zeros(height(dataTimeline.History.Sessions),2) ones(height(dataTimeline.History.Sessions),2).*ymax(2)]',...
             colors(5,:), 'FaceAlpha', 0.2, 'EdgeColor', colors(5,:), 'EdgeAlpha', 0.2, 'LineStyle', '-', 'Marker', 'none')
        box off
    end

    % Link axes
    linkaxes(findall(fig,'type','axes'),'xy')

    % Create directory if needed
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end

    % Save figure and keep visible if requested
    fig.Visible = 'on';    
    saveas(fig, [savepath filesep savename '_Timeline.fig'])
    frame = getframe(fig);
    imwrite(frame.cdata, [savepath filesep savename '_Timeline.png'])
    if ~showFig
        close(fig)
    end
end