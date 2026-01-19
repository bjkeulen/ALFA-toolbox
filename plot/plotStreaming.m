% Author: B.J. Keulen
% Date: 21-10-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function for plotting data from dataStreaming, the output from the
% getStreaming function of the Observe toolbox. Four different subplots are 
% created of time domain data: one with raw data from the original json 
% file and up to three for each type of processing (linenoise filtered, ECG 
% suppressed and both linenoise filtered and ECG suppressed), if performed. 
% Next to that, spectrograms are made of the linenoise filtered and ECG 
% suppressed signal if available, otherwise the linesnoise filtered signal
% is used.
%
% The spectrogram is calculated using a 1s window with 50% overlap and a 
% frequency resolution of ~0.98 Hz. Under the spectrogram the LFP power 
% values are plotted of the selected frequency band defined by the 
% frequency of interest (FOI) and below that the stimulation amplitude is 
% plotted. The LFP thresholds and stimulation limits are also shown if 
% aDBS has been configured.
% 
% INPUT
%   dataRec         =   struct containing one recording of BrainSense 
%                       Streaming data
%   nRuns           =   number (int) of Streaming recordings within JSON
%   ecgMethod       =   int indicating ECG suppression method
%   savepath        =   string with path into which to save figure
%   savename        =   string of name under which to save figure
%   showFig         =   int to determine whether to show the figure or not
%
% OUTPUT
%   none

function plotStreaming(dataRec, nRuns, ecgMethod, savepath, savename, showFig)

    arguments
        dataRec
        nRuns
        ecgMethod
        savepath
        savename
        showFig = 0
    end

    %% Set basic variables
    if ecgMethod == 0
        plots = {'LFP', 'LFP_Linenoise'};
        cols = 1;
    else
        plots = {'LFP', 'LFP_Linenoise', 'LFP_Ecg', 'LFP_LinenoiseEcg'};
        cols = 2;
    end
    titles = {'Raw data', 'Line noise filtered', 'ECG suppressed', 'Line noise filtered and ECG suppressed'};

    %% Time domain plots
    % Create figure
    fig1 = figure('Name','Streaming: LFP','visible','off','WindowState','maximized');
    sgtitle([dataRec.Info.OriginalFile(1:end-5) sprintf('\nBrainSense Streaming - run %i/%i',dataRec.Run,nRuns)],'Interpreter','none')

    % Determine change value of y-axis to separate left and right within plot
    maxvals = zeros([1,4]);
    for p = 1:length(plots)
        maxvals(p) = max(abs(dataRec.Data.(plots{p})),[],'all');
    end
    ychange = max(maxvals);

    % Loop over the different arrays with LFP data
    for p = 1:length(plots)
        subplot(2,cols,p); hold on
        if length(dataRec.Data.Channel) == 2
    
            % Plot data from left electrode
            yyaxis left; plot(dataRec.Data.Time, dataRec.Data.(plots{p})(:,1),'Color',[0 0.4470 0.7410])
            ylim([-3*ychange, 1.2*ychange])
            ylabel('LFP [\muV]','FontSize',12)
    
            % Plot data from right electrode
            yyaxis right; plot(dataRec.Data.Time, dataRec.Data.(plots{p})(:,2),'Color',[0.8500 0.3250 0.0980])
            ylim([-1.2*ychange, 3*ychange])
            ylabel('LFP [\muV]','FontSize',12,'Rotation',270)
    
            % General settings
            title(titles{p},'FontSize',14,'Interpreter','none')
            xlabel('Time [s]','FontSize',12)
            legend(dataRec.Data.Channel,'Interpreter','none','FontSize',10)
            legend('boxoff')

        else
    
            % Plot data from electrode
            plot(dataRec.Data.Time, dataRec.Data.(plots{p})(:,1),'Color',[0 0.4470 0.7410])
            ylabel('LFP [\muV]','FontSize',12)
    
            % General settings
            title(titles{p},'FontSize',14,'Interpreter','none')
            xlabel('Time [s]','FontSize',12)
            legend(dataRec.Data.Channel,'Interpreter','none','FontSize',10)
        end
    end

    % Link axes
    linkaxes(findall(fig1,'type','axes'))

    % Create directory if needed
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end

    % Save figure and keep visible if requested
    fig1.Visible = 'on';
    saveas(fig1, [savepath filesep savename '_LFP.fig'])
    frame1 = getframe(fig1);
    imwrite(frame1.cdata, [savepath filesep savename '_LFP.png'])
    if ~showFig
        close(fig1)
    end

    %% Spectrogram
    % Settings
    fs = dataRec.Settings.SamplingFrequencyTimeDomain;
    window = 1;
    overlap = 0.5;
    colors = orderedcolors("gem");

    % Check recorded and used channels
    if length(dataRec.Data.Channel) == 1
        iChan = 1;
        if any(strcmp(dataRec.Settings.AdaptiveStatus,'NOT_CONFIGURED'))
            if any(contains(dataRec.Settings.SensingChannel,'LEFT'))
                iInfo = 1;
            else
                iInfo = 2;
            end
        else
            iInfo = [1 2];
        end
    else
        iChan = [1 2];
        iInfo = [1 2];
    end

    % Create figure
    ax = gobjects(5,length(iInfo));
    fig2 = figure('Name','Streaming: spectrogram','visible','off','WindowState','maximized');
    sgtitle([dataRec.Info.OriginalFile(1:end-5) sprintf('\nBrainSense Streaming - run %i/%i',dataRec.Run,nRuns)],'Interpreter','none')

    barLims = nan(length(iInfo),2);
    for h = 1:length(iInfo)
        ax(1,h) = subplot(5,length(iInfo),[h h+length(iInfo)]);
        t_empty = [];
        allTime = [];
        allPSD = [];

        if length(iChan) == 1
            iC = 1;
        else
            iC = iChan(h);
        end

        % If LFP_LinenoiseEcg is empty, use LFP_Linenoise
        if all(isnan(dataRec.Data.LFP_LinenoiseEcg(:,iC)))

            % Find segments
            i_seg = [0; find(diff(isnan(dataRec.Data.LFP_Linenoise(:,iC)))); length(dataRec.Data.LFP_Linenoise(:,iC))];

            % Calculate spectrogram for each segment
            for s = 1:2:length(i_seg)
                if i_seg(s+1)-i_seg(s)+1 >= fs*window
                    [~,f,t,psd] = spectrogram(dataRec.Data.LFP_Linenoise(i_seg(s)+1:i_seg(s+1),iC), window*fs, overlap*fs, 0:1/(256/fs):100, fs);
                    if ~isempty(allTime)
                        t_empty = [t_empty allTime(end)+(overlap*window):(overlap*window):dataRec.Data.Time(i_seg(s)+1)];
                    end
                    allTime = [allTime t + dataRec.Data.Time(i_seg(s)+1)];
                    allPSD = [allPSD psd];
                end
            end

            % Check whether there is data present
            if isempty(allPSD)
                continue
            end

            % Insert empty data for NaN segments
            psd_empty = nan(length(f),round(sum(diff(allTime)-(overlap*window))/(overlap*window)));
            [allTime, i_sort] = sort([allTime t_empty]);
            allPSD = [allPSD psd_empty];
            allPSD = allPSD(:,i_sort);
    
            % Plot figure
            barLims(h,:) = [min(pow2db(allPSD),[],'all') max(pow2db(allPSD),[],'all')];
            imagesc(allTime,f,pow2db(allPSD))
            set(gca,'YDir','normal')

            % Set title and y-label
            if length(iInfo) == 2
                if h == 1
                    lead = 'Left electrode';
                else
                    lead = 'Right electrode';
                end
            else
                if contains(dataRec.Settings.SensingChannel{1}, 'LEFT')
                    lead = 'Left electrode';
                else
                    lead = 'Right electrode';
                end
            end
            title(sprintf('%s | Sensing channel %s\n%s | FOI = %.2f Hz', lead, dataRec.Settings.SensingChannel{iInfo(h)}, ...
                  titles{2}, dataRec.Settings.SensingFrequency(iInfo(h))),'FontSize',14,'Interpreter','none')
            ylabel('Freqency [Hz]','FontSize',12)

            % Set colorbar
            cb = colorbar(ax(1,h),'Limits',[min(barLims(:,1)) max(barLims(:,2))]);
            ylabel(cb,'Power [dB]','FontSize',12,'Rotation',270)

            % Set positions of spectrograms to align subplots
            if length(iInfo) == 2
                if h == 1
                    ax(1,1).Position = [0.1300,0.5927,0.3347,0.2803];
                else
                    ax(1,2).Position = [0.5703,0.5927,0.3347,0.2803];
                end
            else
                ax(1,1).Position = [0.1300,0.5927,0.7750,0.2803];
            end

        % Use LFP_LinenoiseEcg
        else

            % Find segments
            i_seg = [0; find(diff(isnan(dataRec.Data.LFP_LinenoiseEcg(:,iC)))); length(dataRec.Data.LFP_LinenoiseEcg(:,iC))];

            % Calculate spectrogram for each segment
            for s = 1:2:length(i_seg)
                if i_seg(s+1)-i_seg(s)+1 >= fs*window
                    [~,f,t,psd] = spectrogram(dataRec.Data.LFP_LinenoiseEcg(i_seg(s)+1:i_seg(s+1),iC), window*fs, overlap*fs, 0:1/(256/fs):100, fs);
                    if ~isempty(allTime)
                        t_empty = [t_empty allTime(end)+(overlap*window):(overlap*window):dataRec.Data.Time(i_seg(s)+1)];
                    end
                    allTime = [allTime t + dataRec.Data.Time(i_seg(s)+1)];
                    allPSD = [allPSD psd];
                end
            end

            % Check whether there is data present
            if isempty(allPSD)
                continue
            end

            % Insert empty data for NaN segments
            psd_empty = nan(length(f),round(sum(diff(allTime)-(overlap*window))/(overlap*window)));
            [allTime, i_sort] = sort([allTime t_empty]);
            allPSD = [allPSD psd_empty];
            allPSD = allPSD(:,i_sort);
    
            % Plot figure
            barLims(h,:) = [min(pow2db(allPSD),[],'all') max(pow2db(allPSD),[],'all')];
            imagesc(allTime,f,pow2db(allPSD))
            set(gca,'YDir','normal')

            % Set title and y-label
            if length(iInfo) == 2
                if h == 1
                    lead = 'Left electrode';
                else
                    lead = 'Right electrode';
                end
            else
                if contains(dataRec.Settings.SensingChannel{1}, 'LEFT')
                    lead = 'Left electrode';
                else
                    lead = 'Right electrode';
                end
            end
            title(sprintf('%s | Sensing channel %s\n%s | FOI = %.2f Hz', lead, dataRec.Settings.SensingChannel{iInfo(h)}, ...
                  titles{4}, dataRec.Settings.SensingFrequency(iInfo(h))),'FontSize',14,'Interpreter','none')
            ylabel('Freqency [Hz]','FontSize',12) 

            % Set colorbar
            cb = colorbar(ax(1,h),'Limits',[min(barLims(:,1)) max(barLims(:,2))]);
            ylabel(cb,'Power [dB]','FontSize',12,'Rotation',270)

            % Set positions of spectrograms to align subplots
            if length(iInfo) == 2
                if h == 1
                    ax(1,1).Position = [0.1300,0.5927,0.3347,0.2803];
                else
                    ax(1,2).Position = [0.5703,0.5927,0.3347,0.2803];
                end
            else
                ax(1,1).Position = [0.1300,0.5927,0.7750,0.2803];
            end

        end

        % Plot power values, use mask to reduce influence of outliers on y-limit
        ax(2,h) = subplot(5,length(iInfo),[h+2*length(iInfo) h+2*length(iInfo)+length(iInfo)]);
        plot(dataRec.Data.LfpPower.Time, dataRec.Data.LfpPower.Value(:,iC),'Color',colors(1,:)); hold on;

        dataMask = dataRec.Data.LfpPower.Value(dataRec.Data.LfpPower.Value <= mean(dataRec.Data.LfpPower.Value,'all','omitnan') + 5*std(dataRec.Data.LfpPower.Value,[],'all','omitnan'));
        ylim([0 1.1*max(dataMask,[],'all')])
        xlim([0 max(dataRec.Data.LfpPower.Time)])
        ylabel('LFP power [a.u.]','FontSize',12)

        if ~strcmp(dataRec.Settings.AdaptiveStatus{iInfo(h)},'NOT_CONFIGURED')
            yline(dataRec.Settings.LowerLfpThreshold(:,iInfo(h)),'LineStyle','-.','LineWidth',1.5,'Color',colors(3,:))
            yline(dataRec.Settings.UpperLfpThreshold(:,iInfo(h)),'LineStyle','-.','LineWidth',1.5,'Color',colors(3,:))
        end

        % Set position to align subplots
        if length(iInfo) == 1
            ax(2,1).Position = [0.1300,0.2668,0.7750,0.2803];
        end

        % Plot stimulation amplitudes
        ax(3,h) = subplot(5,length(iInfo),h+4*length(iInfo));
        plot(dataRec.Data.Stimulation.Time, dataRec.Data.Stimulation.Amplitude(:,iInfo(h)),'Color',colors(2,:))
        xlim([0 max(dataRec.Data.Stimulation.Time)])
        xlabel('Time [s]','FontSize',12)
        ylabel({'Stimulation';'amplitude [mA]'},'FontSize',12)

        if ~strcmp(dataRec.Settings.AdaptiveStatus{iInfo(h)},'NOT_CONFIGURED')
            yline(dataRec.Settings.LowerStimulationLimit(:,iInfo(h)),'LineStyle',':','Color',colors(2,:))
            yline(dataRec.Settings.UpperStimulationLimit(:,iInfo(h)),'LineStyle',':','Color',colors(2,:))
            ylim([dataRec.Settings.LowerStimulationLimit(:,iInfo(h))-0.2 dataRec.Settings.UpperStimulationLimit(:,iInfo(h))+0.2])
        else
            ylim([min(dataRec.Data.Stimulation.Amplitude(:,iInfo(h)))-0.2 max(dataRec.Data.Stimulation.Amplitude(:,iInfo(h)))+0.2])
        end
    end

    % Check presence of data
    if isempty(findobj(fig2,'-property','XData'))
        warning('Insufficient data in %s for calculation of spectrogram.',savename)
    else
    
        % Link axes
        warning('off','MATLAB:linkaxes:RequireDataAxes')
        linkaxes(ax,'x')
        linkaxes(ax(1,:),'yz')
        linkaxes(ax(2,:),'y')
    
        % Save figure and keep visible if requested
        fig2.Visible = 'on';
        saveas(fig2, [savepath filesep savename '_spectrogram.fig'])
        frame2 = getframe(fig2);
        imwrite(frame2.cdata, [savepath filesep savename '_spectrogram.png'])
        if ~showFig
            close(fig2)
        end
    end
end