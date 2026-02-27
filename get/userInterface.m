% Author: B.J. Keulen
% Date: 10-11-2025
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% User interface to collect settings for LFP processing within the ALFA
% toolbox. Setting regard the type of dataset, ECG suppression method, line
% noise filter, timezone, whether data should be plotted and whether 
% figures should be shown.
% 
% INPUT
%   none
%
% OUTPUT
%   settings        =   struct containing settings

function settings = userInterface()

    % Create GUI window
    ui = uifigure('Name','Settings','WindowStyle','alwaysontop','Resize','off','Position',[600 500 420 490]);
    ui.CloseRequestFcn = @(src,event) cancel();

    % Dataset type
    bg_dataset = uibuttongroup(ui,'Position',[30 365 160 105],'Title','Type of dataset','FontWeight','bold','BorderType','none');
    b1_dataset = uiradiobutton(bg_dataset,'Position',[10 55 150 25],'Text','Single file');
    b2_dataset = uiradiobutton(bg_dataset,'Position',[10 30 150 25],'Text','Folder');
    b3_dataset = uiradiobutton(bg_dataset,'Position',[10 5 150 25],'Text','Set of folders');

    % Figure settings
    bg_fig = uibuttongroup(ui,'Position',[30 250 160 105],'Title','Figure settings','FontWeight','bold','BorderType','none');
    b1_fig = uicheckbox(bg_fig,'Position',[10 55 150 25],'Text',' Create figures','Value',1);
    b2_fig = uicheckbox(bg_fig,'Position',[10 30 150 25],'Text',' Show figures','Value',1);

    % Timezone
    bg_time = uibuttongroup(ui,'Position',[30 90 160 170],'Title','Timezone','FontWeight','bold','BorderType','none','SelectionChangedFcn',@(src, event) enableZoneSelector(event));
    b1_time = uiradiobutton(bg_time,'Position',[10 120 150 25],'Text','Local (system)','Value',true);
    b2_time = uiradiobutton(bg_time,'Position',[10 95 150 25],'Text','Other','Value',false);
    area_time = uidropdown(bg_time,'Position',[25 65 130 25],'Enable','off','Items',unique(timezones().Area),'ValueChangedFcn',@(src,event) updateTimeZoneArea(src));
    zone_time = uidropdown(bg_time,'Position',[25 35 130 25],'Enable','off','Items',timezones('Africa').Name);
    b3_time = uiradiobutton(bg_time,'Position',[10 5 150 25],'Text','JSON default (UTC)','Value',false);

    % ECG suppression
    bg_ecg = uibuttongroup(ui,'Position',[230 285 160 185],'Title','ECG suppression method','FontWeight','bold','BorderType','none');
    b1_ecg = uiradiobutton(bg_ecg,'Position',[10 135 150 25],'Text','None','Value',false);
    b2_ecg = uiradiobutton(bg_ecg,'Position',[10 110 150 25],'Text','Automatic','Value',true);
    b3_ecg = uiradiobutton(bg_ecg,'Position',[10 85 150 25],'Text','Manual','Value',false);

    uilabel(bg_ecg,'Position',[10 50 120 25],'Text','R-peak window [ms]','HorizontalAlignment','center');
    uilabel(bg_ecg,'Position',[15 5 50 20],'Text','Before','FontAngle','italic','HorizontalAlignment','center');
    b4_ecg = uieditfield(bg_ecg,'numeric','Position',[15 25 50 25],'Value',250,'Limits',[50 500]);
    uilabel(bg_ecg,'Position',[75 5 50 20],'Text','After','FontAngle','italic','HorizontalAlignment','center');
    b5_ecg = uieditfield(bg_ecg,'numeric','Position',[75 25 50 25],'Value',400,'Limits',[50 500]);

    % Line noise filter
    bg_line = uibuttongroup(ui,'Position',[230 90 160 170],'Title','Line noise filter','FontWeight','bold','BorderType','none');
    uilabel(bg_line,'Position',[10 120 120 25],'Text','Center frequency [Hz]');
    b1_line = uiradiobutton(bg_line,'Position',[10 95 50 25],'Text','50 Hz','Value',true);
    b2_line = uiradiobutton(bg_line,'Position',[80 95 50 25],'Text','60 Hz','Value',false);
    uilabel(bg_line,'Position',[10 55 120 25],'Text','Bandwidth [Hz]');
    bw_line = uislider(bg_line,'Position',[10 40 130 3],'Value',0.2,'Limits',[0 1],'MajorTicks',0:0.1:1,'MajorTickLabels',{'0.0','','0.2','','','0.5','','','','','1.0'},'MinorTicks',[]);

    function enableZoneSelector(event)
        if strcmp(event.NewValue.Text,'Other')
            area_time.Enable = true;
            zone_time.Enable = true;
        else
            area_time.Enable = false;
            zone_time.Enable = false;
        end
    end

    function updateTimeZoneArea(src)
        zone = src.Value;
        zone_time.Items = timezones(zone).Name;
    end

    % Create cancel and confirm button
    bg_done = uibuttongroup(ui,'Position',[30 15 360 75],'Title',' ','BorderType','none');
    uibutton(bg_done,'push','Text','Cancel','Position',[90 5 70 35],'FontWeight','bold','ButtonPushedFcn',@(btn,event) cancel);
    uibutton(bg_done,'push','Text','Confirm','Position',[200 5 70 35],'FontWeight','bold','ButtonPushedFcn',@(btn,event) confirm);

    % Wait to give output
    waitfor(ui)

    % In case window was closed or input canceled
    function cancel()
        settings = struct;
        closereq()
    end

    % If input was confirmed
    function confirm()

        % Initialise struct
        settings = struct;

        % Dataset
        if b1_dataset.Value
            settings.dataset = 0;
        elseif b2_dataset.Value
            settings.dataset = 1;
        elseif b3_dataset.Value
            settings.dataset = 2;
        end

        % Figures
        settings.plotData = b1_fig.Value;
        settings.showFig = b2_fig.Value;

        % Timezone
        if b1_time.Value
            settings.tZone = datetime().SystemTimeZone;
        elseif b2_time.Value
            settings.tZone = zone_time.Value;
        end

        % ECG method
        settings.rWindow = [nan nan];
        if b1_ecg.Value
            settings.ecgMethod = 0;
        else
            if b2_ecg.Value
                settings.ecgMethod = 1;
            elseif b3_ecg.Value
                settings.ecgMethod = 2;
            end
            settings.rWindow(1) = b4_ecg.Value/1000;
            settings.rWindow(2) = b5_ecg.Value/1000;
        end

        % Line noise
        if b1_line.Value
            settings.linenoise = [50 bw_line.Value];
        elseif b2_line.Value
            settings.linenoise = [60 bw_line.Value];
        end

        % Close window
        closereq()
    end
end