% Author: B.J. Keulen
% Date: 25-10-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for plotting data from dataEvents, the output from the
% mainExtract function of the Observe toolbox. Each type of event is
% plotted as the mean PSD +- standard deviation.
% 
% INPUT
%   dataEvents      =   struct containing BrainSense Events
%   savepath        =   string with path into which to save figure
%   savename        =   string of name under which to save figure
%   showFig         =   int to determine whether to show the figure or not
%
% OUTPUT
%   none

function plotEvents(dataEvents, savepath, savename, showFig)

    arguments
        dataEvents
        savepath
        savename
        showFig = 0
    end

    % Check presence of PSD data
    if isempty([dataEvents.Data{:,'Magnitude'}{:}])
        warning('no PSD data available in %s.\nFigure will not be created.', [savepath filesep savename '_Events'])
        return
    end

    % Retrieve unique events
    events = unique(dataEvents.Data{:,'EventName'});

    % Create array for saving variables for plotting
    avg_l = nan([length(events),100]);
    sd_l = nan([length(events),100]);
    avg_r = nan([length(events),100]);
    sd_r = nan([length(events),100]);

    % Initiate array for counting valid recordings
    num = nan(length(events),2);

    % Loop over events. Retrieve, check and plot data
    iData = 1:length(events);
    for e = 1:length(events)
        data = dataEvents.Data(strcmp(dataEvents.Data{:,'EventName'}, events{iData(e)}),:);
        
        % Check channels and values
        if isempty([data{:,'Magnitude'}{:,1}]) && isempty([data{:,'Magnitude'}{:,2}])
            warning('no PSD data available for event "%s".\n Figure: %s.', events{iData(e)}, [savepath filesep savename '_Events'])
            iData(e) = nan;
            continue
        end

        % Collect data
        left = nan([height(data),100]);
        right = left;
        for i = 1:height(data)
            if ~isempty(data{i,'Magnitude'}{1})
                if length(data{i,'Magnitude'}{1}) == 100
                    left(i,:) = data{i,'Magnitude'}{1}';
                end
            end
            if ~isempty(data{i,'Magnitude'}{2})
                if length(data{i,'Magnitude'}{2}) == 100
                    right(i,:) = data{i,'Magnitude'}{2}';
                end
            end
        end

        % Calculate mean and standard deviation
        avg_l(e,:) = mean(left,1,"omitmissing");
        sd_l(e,:) = std(left,0,1,"omitmissing");
        avg_r(e,:) = mean(right,1,"omitmissing");
        sd_r(e,:) = std(right,0,1,"omitmissing");

        % Count valid recordings
        num(e,1) = sum(~isnan(left(:,1)));
        num(e,2) = sum(~isnan(right(:,1)));

        % Get frequencies
        j = 1;
        while isempty(data{j,'Frequency'}{:}) || length(data{j,'Frequency'}{:}) ~= 100
            j = j + 1;
        end
        freq = data{j,'Frequency'}{:};
    end

    % Plot all events in one figure
    fig = figure('Name','Events','visible','off','WindowState','maximized'); hold on
    sgtitle([savename ' - BrainSense Events'],'Interpreter','none')
    colors = [orderedcolors("gem"); orderedcolors("glow")];
    labelsL = {}; labelsR = {};

    % Loop over events
    iData = iData(~isnan(iData));
    for e = 1:length(iData)

        % Plot event, left hemisphere
        subplot(2,1,1); hold on
        plot(freq, avg_l(iData(e),:), 'Color', colors(e,:), 'LineWidth', 1.5, 'DisplayName', events{iData(e)})
        fill([freq', fliplr(freq')], [avg_l(iData(e),:)-sd_l(iData(e),:), fliplr(avg_l(iData(e),:)+sd_l(iData(e),:))], colors(e,:), 'EdgeColor', 'none')
        alpha(0.10)
        box off
        title('Left electrode','FontSize', 14)
        xlabel('Frequency [Hz]','FontSize', 12)
        ylabel('Magnitude [\muVp]','FontSize', 12)

        % Set legend
        if all(isnan(avg_l(iData(e),:)))
            labelsL = [labelsL {'', ''}];
        else
            labelsL = [labelsL {sprintf('%s (n=%i)',events{iData(e)},num(iData(e),1)),''}];
        end
        if e == length(iData)
            legend(labelsL,'FontSize', 11)
            legend('boxoff')
        end

        % Plot event, right hemisphere
        subplot(2,1,2); hold on
        plot(freq, avg_r(iData(e),:), 'Color', colors(e,:), 'LineWidth', 1.5, 'DisplayName', events{iData(e)})
        fill([freq', fliplr(freq')], [avg_r(iData(e),:)-sd_r(iData(e),:), fliplr(avg_r(iData(e),:)+sd_r(iData(e),:))], colors(e,:), 'EdgeColor', 'none')
        alpha(0.10)
        box off
        title('Right electrode','FontSize', 14)
        xlabel('Frequency [Hz]','FontSize', 12)
        ylabel('Magnitude [\muVp]','FontSize', 12)

        % Set legend
        if all(isnan(avg_r(iData(e),:)))
            labelsR = [labelsR {'', ''}];
        else
            labelsR = [labelsR {sprintf('%s (n=%i)',events{iData(e)},num(iData(e),2)),''}];
        end
        if e == length(iData)
            legend(labelsR,'FontSize', 11)
            legend('boxoff')
        end
    end

    % Set y-limit
    ax = findall(fig,'type','axes');
    ylim(ax, [0 ceil(max([avg_l+sd_l; avg_r+sd_r],[],'all'))])
    linkaxes(ax)

    % Create directory if needed
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end

    % Save figure and keep visible if requested
    fig.Visible = 'on';
    saveas(fig, [savepath filesep savename '_Events.fig'])
    frame = getframe(fig);
    imwrite(frame.cdata, [savepath filesep savename '_Events.png'])
    if ~showFig
        close(fig)
    end
end