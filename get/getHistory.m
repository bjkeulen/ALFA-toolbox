% Author: B.J. Keulen and J.T. Boonstra
% Date: 17-02-2025
%
% Copyright 2025 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function to retrieve and overview of all sessions and group changes
% within the time period of dataTimeline. First, all known sessions which
% are already stored within dataTimeline.Info are put within a single table
% containing all stimulation and sensing parameters with the exception of
% stimulation amplitudes and LFP values, which are already within
% dataTimeline. Then, all sessions of which no JSON is available (missing
% sessions) and all group changes are collected from the EventLogs. The
% initial and final stimulation and sensing parameters of missing sessions
% are assumed to be equal to the session before and after each session, 
% respectively.
% 
% INPUT
%   dataTimeline    =   struct with Timeline data
%
% OUTPUT
%   sessions        =   table containing initial and final stimulation
%                       and sensing parameters of all sessions
%   changes         =   table containing all group changes

function [sessions, changes] = getHistory(dataTimeline)

    % Collect all general group information from info
    sessions = repmat(struct(),1,length(dataTimeline.Info));
    for i = 1:length(dataTimeline.Info)
        sessions(i).SessionType = {'JSON Session Report'};
        sessions(i).SessionStart = datetime(dataTimeline.Info(i).SessionStartDateUtc,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
        sessions(i).SessionEnd = datetime(dataTimeline.Info(i).SessionEndDateUtc,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
        sessions(i).Initial = dataTimeline.Info(i).GroupsInitial;
        sessions(i).Final = dataTimeline.Info(i).GroupsFinal;
    end

    % Get unique sessions and sort
    sessions = struct2table(sessions, "AsArray", true);
    [~, idx, ~] = unique(sessions.SessionStart);
    sessions = sortrows(sessions(idx,:),'SessionStart');

    % Retrieve intial and final values of channels and frequencies
    for h = 1:height(sessions)

        % Initiate structures for initial values
        sessions.InitialActive(h) = "";
        sessions.InitialID(h) = {nan};
        sessions.InitialSensingChannel(h) = {cell(4,2)};
        sessions.InitialSensingFrequency(h) = {nan([4,2])};
        sessions.InitialLowerLfpThreshold(h) = {nan([4,2])};
        sessions.InitialUpperLfpThreshold(h) = {nan([4,2])};
        sessions.InitialLowerAmplitudeLimit(h) = {nan([4,2])};
        sessions.InitialUpperAmplitudeLimit(h) = {nan([4,2])};
        sessions.InitialStimulationFrequency(h) = {nan([4,2])};
        sessions.InitialPulseWidth(h) = {nan([4,2])};

        % Collect initial data 
        try
            initial = sessions.Initial{h};
        catch
            initial = sessions.Initial(h);
        end

        % Loop over groups and retrieve data
        for g = 1:length(initial)

            % Retrieve group info
            [rateI, pulseI, chanI, freqI, ampLI, ampUI, lfpLI, lfpUI] = getGroupInfo(initial(g));

            % Get group ID and store data
            grId = str2double(replace(initial(g).GroupId,{'GroupIdDef.GROUP_','A','B','C','D'},{'','1','2','3','4'}));
            sessions.InitialSensingChannel{h}(grId,:) = chanI;
            sessions.InitialSensingFrequency{h}(grId,:) = freqI;
            sessions.InitialLowerLfpThreshold{h}(grId,:) = lfpLI;
            sessions.InitialUpperLfpThreshold{h}(grId,:) = lfpUI;
            sessions.InitialLowerAmplitudeLimit{h}(grId,:) = ampLI;
            sessions.InitialUpperAmplitudeLimit{h}(grId,:) = ampUI;
            sessions.InitialStimulationFrequency{h}(grId,:) = rateI;
            sessions.InitialPulseWidth{h}(grId,:) = pulseI;

            % Check if active group, store if so
            if initial(g).ActiveGroup
                sessions.InitialActive{h} = initial(g).GroupId(end-6:end);
                sessions.InitialID{h} = grId;
            end
        end

        % Initiate structures for final values
        sessions.FinalActive(h) = "";
        sessions.FinalID(h) = {nan};
        sessions.FinalSensingChannel(h) = {cell(4,2)};
        sessions.FinalSensingFrequency(h) = {nan([4,2])};
        sessions.FinalLowerLfpThreshold(h) = {nan([4,2])};
        sessions.FinalUpperLfpThreshold(h) = {nan([4,2])};
        sessions.FinalLowerAmplitudeLimit(h) = {nan([4,2])};
        sessions.FinalUpperAmplitudeLimit(h) = {nan([4,2])};
        sessions.FinalStimulationFrequency(h) = {nan([4,2])};
        sessions.FinalPulseWidth(h) = {nan([4,2])};

        % Collect final data
        try
            final = sessions.Final{h};
        catch
            final = sessions.Final(h);
        end

        % Loop over groups and retrieve data
        for g = 1:length(final)
            
            % Retrieve group info
            [rateF, pulseF, chanF, freqF, ampLF, ampUF, lfpLF, lfpUF] = getGroupInfo(final(g));

            % Get group ID and store data
            grId = str2double(replace(final(g).GroupId,{'GroupIdDef.GROUP_','A','B','C','D'},{'','1','2','3','4'}));
            sessions.FinalSensingChannel{h}(grId,:) = chanF;
            sessions.FinalSensingFrequency{h}(grId,:) = freqF;
            sessions.FinalLowerLfpThreshold{h}(grId,:) = lfpLF;
            sessions.FinalUpperLfpThreshold{h}(grId,:) = lfpUF;
            sessions.FinalLowerAmplitudeLimit{h}(grId,:) = ampLF;
            sessions.FinalUpperAmplitudeLimit{h}(grId,:) = ampUF;
            sessions.FinalStimulationFrequency{h}(grId,:) = rateF;
            sessions.FinalPulseWidth{h}(grId,:) = pulseF;

            % Check if active group, store if so
            if final(g).ActiveGroup
                sessions.FinalActive{h} = final(g).GroupId(end-6:end);
                sessions.FinalID{h} = grId;
            end
        end
    end

    % Remove initial and final
    sessions(:,[4 5]) = [];

    % Define structures for missing sessions and changes
    missLogs = table('Size',[0,2],'VariableNames',{'SessionStart','SessionEnd'},'VariableTypes',{'datetime','datetime'});
    changes = table('Size',[0,3],'VariableNames',{'DateTime','OldGroup','NewGroup'},'VariableTypes',{'datetime','string','string'});

    % Loop over EventLogs
    for i = 1:height(dataTimeline.Info)
        for e = 1:length(dataTimeline.Info(i).EventLogs)
            try
                logsE = dataTimeline.Info(i).EventLogs{e};
            catch
                logsE = dataTimeline.Info(i).EventLogs{e};
            end

            % Get ends of sessions
            if isfield(logsE,'SessionType')
                if logsE.SessionType == "SessionStateDef.SessionStart"
                    stop = NaT(1);
                    for b = e:length(dataTimeline.Info(i).EventLogs)
                        try
                            logsB = dataTimeline.Info(i).EventLogs{b};
                        catch
                            logsB = dataTimeline.Info(i).EventLogs{b};
                        end
                        if isfield(logsB,'SessionType')
                            if logsB.SessionType == "SessionStateDef.SessionEnd"
                                if ~strcmp(logsE.DateTime, logsB.DateTime)
                                    start = datetime(logsE.DateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
                                    stop = datetime(logsB.DateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
                                    break
                                end
                            end
                        end
                    end
                    if ~any(max(missLogs{:,'SessionStart'}, start) < min(missLogs{:,'SessionEnd'}, stop))
                        missLogs = [missLogs; {start, stop}];
                    end
                end
            end

            % Get group changes
            if isfield(logsE,'OldGroupId')
                    if ~any(changes{:,'DateTime'} == datetime(logsE.DateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z'''))
                        changes = [changes; {datetime(logsE.DateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z'''), ...
                                   logsE.OldGroupId(end-6:end), logsE.NewGroupId(end-6:end)}];
                    end
            end
        end
    end

    % Remove sessions until most recent after start dataTimeline, and remove when already in history
    missLogs = sortrows(missLogs,'SessionStart');
    missLogs = missLogs(find(missLogs{:,"SessionStart"} > dataTimeline.Data{1,'DateTime'},1,'first'):end,:);
    for m = height(missLogs):-1:1
        if any(max(sessions{:,'SessionStart'}, missLogs{m, 'SessionStart'}) < min(sessions{:,'SessionEnd'}, missLogs{m, 'SessionEnd'}))
            missLogs(m,:) = [];
        end
    end

    % Find first and last indices of adjacent missing sessions
    iM = nan([height(missLogs),1]);
    for m = 1:height(missLogs)
        iM(m) = find(missLogs{m,"SessionEnd"} < sessions.SessionStart,1);
    end
    [iU, iF, ~] = unique(iM);
    [~, iL, ~] = unique(iM,'last');

    % Fill and sort table
    missing = table('Size',[height(missLogs),width(sessions)],'VariableNames',sessions.Properties.VariableNames,'VariableTypes',varfun(@class,sessions,'OutputFormat','cell'));
    missing(:,1) = {'Missing Session'};
    missing(:,2:3) = missLogs;
    for i = 1:height(iU)
        missing(iL(i),16:23) = sessions(iU(i),6:13);
        if iU(i) > 1
            missing(iF(i),6:13) = sessions(iU(i)-1,16:23);
        end
    end

    % Remove changes with the same datetime and different groups
    iC = true(height(changes),1);
    if ~isempty(changes)
        for b = 2:height(changes)
            if isequal(changes.DateTime(b), changes.DateTime(b-1))
                if b > 2
                    iC([b-1 b]) = strcmp(changes.NewGroup{b-2},{changes.OldGroup{b-1}; changes.OldGroup{b}});
                elseif b > 1
                    iC([b-1 b]) = strcmp(changes.OldGroup{b+1},{changes.NewGroup{b-1}; changes.NewGroup{b}});
                end
            end
        end
    end
    changes = changes(iC,:);

    % Assign groups and group IDs to missing sessions
    for m = 1:height(missing)
        if isempty(changes)
            im = find(missing{m,"SessionEnd"} < sessions{:,"SessionStart"}, 1);
            missing{m,"FinalActive"} = sessions{im,"InitialActive"};
            missing{m,"FinalID"} = sessions{im,"InitialID"};
            if im == 1
                break
            else
                im = find(missing{m,"SessionStart"} > sessions{:,"SessionEnd"}, 1, "last");
                if isempty(im)
                    break
                else
                    missing{m,"InitialActive"} = sessions{im,"FinalActive"};
                    missing{m,"InitialID"} = sessions{im,"FinalID"};
                end
            end
        else
            im = find(missing{m,"SessionStart"} > changes{:,"DateTime"}, 1, "last");
            if isempty(im)
                missing{m,"InitialActive"} = changes{1,"OldGroup"};
                missing{m,"InitialID"} = {str2double(replace(missing{m,"InitialActive"},{'GROUP_','A','B','C','D'},{'','1','2','3','4'}))};
                missing{m,"FinalActive"} = missing{m,"InitialActive"};
                missing{m,"FinalID"} = missing{m,"InitialID"};
            else
                missing{m,"InitialActive"} = changes{im,"NewGroup"};
                missing{m,"InitialID"} = {str2double(replace(missing{m,"InitialActive"},{'GROUP_','A','B','C','D'},{'','1','2','3','4'}))};

                im = find(missing{m,"SessionEnd"} > changes{:,"DateTime"}, 1, "last");
                missing{m,"FinalActive"} = changes{im,"NewGroup"};
                missing{m,"FinalID"} = {str2double(replace(missing{m,"FinalActive"},{'GROUP_','A','B','C','D'},{'','1','2','3','4'}))};
            end
        end
    end

    % Add missing to sessions and sort
    if ~isempty(missing)
        sessions = sortrows([sessions; missing],'SessionStart');
    end
end