% Author: B.C.M van Wijk
% Date: 02-02-2018
% Adjusted by: M.J. Stam
% Date: 01-10-2021
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to suppress ECG-artifact in LFP time series, described in-depth
% in Stam et al. 2023:
% 
% Stam MJ, van Wijk BCM, Sharma P, Beudel M, Pina-Fuentes DA, de Bie RMA, 
% Schuurman PR, Neumann W-J, Buijink AWG. A comparison of methods to 
% suppress electrocardiographic artifacts in local field potential 
% recordings. Clin Neurophysiol. 2023 Feb:146:147-161. DOI: 
% 10.1016/j.clinph.2022.11.011
% 
% Epochs the data around each heartbeat (R-peak) using Medtronic code. 
% Performs a singular value decomposition across epochs. Removes selected
% ncomp from each epoch. Corrects for sudden jumps at start and end of each 
% epoch. Output is the corrected time series, ECG artifact subtracted per 
% epoch, number of components chosen to generate ECG template, and index 
% numbers of selected peaks.
%
% INPUT (REQUIRED):
%   sigin                                 % LFP signal containing ECG artifact
%   settings.ArtTimeBeforePeak            % time in s before R peak (try to cover P)
%   settings.ArtTimeAfterPeak             % time in s after R peak (try to cover T)
%   settings.nComp                        % number of SVD components used to reconstruct artifact.
%                                         Could be increased if there is a lot of residual ECG left.
%                                         Note: ncomp > 20 is interpreted as a target percentage
%                                         of explained variance, e.g. 95 = 95%
%   settings.SamplingFrequency            % sampling rate in Hz
%   settings.Polarity                     % 1 if R peak points upwards, -1 if R peak points downwards
%
% INPUT (OPTIONAL):
%   ecg                                   % synced ECG timeseries for visualization
%   ecg_peak_indices                      % index numbers of pre-selected peaks of input LFP signal
%   ecg_pks                               % index numbers of pre-selected peaks of synced ECG timeseries
%   settings.ShowFigs                     % plot figures to monitor correction
%   settings.SaveFigs                     % save figures to disk
%   settings.Label                        % label used to save figure
%   settings.Title                        % title displayed as figure title
%   settings.Interactive                  % enable manual adjustments of detected peaks
%   settings.Finalcomp                    % enable manual adjustments of selected SVD components
%
% OUTPUT:
%   sigin                                 % LFP signal after suppressing ECG artifact
%   proj_out                              % projections used for supressing ECG artifact for each R-peak in input LFP signal
%   finalcomp                             % final components chosen to suppress ECG artifact
%   ecg_peak_indices                      % final index numbers of R-peaks of ECG artifact in input LFP signal
%
% -------------------------------------------------------------------------
% NOTE: perform this algorithm on a slighly longer time series as intended
% for analysis so first few seconds can be cut off due to filter
% -------------------------------------------------------------------------

function [sigout, proj_out, finalcomp, ecg_peak_indices] = SVD_ECG_Filter(sigin, settings, ecg,ecg_peak_indices, ecg_pks)

    if nargin == 2
        ecg = [];
        ecg_peak_indices = [];
    elseif nargin == 3
        ecg_peak_indices = [];
    end
    if size(sigin,2)>size(sigin,1)
        sigin = sigin';
    end
    
    % Mandatory fields
    Fs = settings.SamplingFrequency;
    polarity = settings.Polarity;
    art_time_b4_peak = settings.ArtTimeBeforePeak;
    art_time_after_peak = settings.ArtTimeAfterPeak;
    ncomp = settings.nComp;
    
    % Other fields
    if isfield(settings,'ShowFigs')
        showfigs = settings.ShowFigs;
    else
        showfigs = 1;
    end
    if isfield(settings,'SaveFigs')
        savefigs = settings.SaveFigs;
    else
        savefigs = 0;
    end
    if isfield(settings,'Interactive')
        interactive = settings.Interactive;
    else
        interactive = 1;
    end
    if isfield(settings,'Label')
        label = settings.Label;
    else
        label = 'test';
    end
    if isfield(settings,'Title')
        figname = settings.Title;
    else
        figname = 'test';
    end
    % Check whether co-recorded synced ECG timeseries is present:
    if isempty(ecg)
        ecg = [];
        plotecg = 0;
        if size(ecg,2)>size(ecg,1)
            ecg = ecg';
        end
    else
        plotecg = 1;
        scaleecg = std(ecg,'omitnan')/std(sigin,'omitnan');
        ecg = (ecg - mean(ecg,'omitnan'))./scaleecg;
    end
    if isempty(ecg_peak_indices)
        ecg_peak_indices = [];
    end
    
    % Conversion of the parameters in number of samples
    art_pts_b4_peak = round(art_time_b4_peak * Fs);
    art_pts_after_peak = round(art_time_after_peak * Fs);
    art_pts_plus_t_after_peak = art_pts_after_peak;
    art_width_pts = art_pts_b4_peak + art_pts_plus_t_after_peak + 1;
    
    % Initialization
    artifact_count = 0;
    startedgepeak = [];
    endedgepeak = [];
    NaNpeak = [];
    time = (1:length(sigin))/Fs;
    
    % If no pre-selected R-peaks are indicated, ECG suppression starts with R-peak detection
    % method described in Stam et al. 2023: https://doi.org/10.1016/j.clinph.2022.11.011
    if isempty(ecg_peak_indices)
        LFPecg = sigin;
    
        % Normalize and find peaks
        LFPnorm = normalize(LFPecg);
        [Rpeak,locs_Rwave] = findpeaks(LFPnorm,'MinPeakHeight',(2.5*std(LFPnorm)),...
            'MinPeakDistance',Fs/2);

        % Also search the flipped signal for R-peaks:
        [Speak,locs_Swave] = findpeaks(-LFPnorm,'MinPeakHeight',(2.5*std(-LFPnorm)),...
            'MinPeakDistance',Fs/2);
    
        % Determine which way the peaks are highest, corresponding to R-peaks:
        if isempty(Rpeak)
            Rpeak = 0;
        end
        if isempty(Speak)
            Speak = 0;
        end

        % Flip polarity if the mean height and the number of detected R-peaks
        % are higher with negative polarity
        if mean(Speak) > mean(Rpeak) && length(Speak) >= length(Rpeak)
            locs_Rwave = locs_Swave;
            polarity = -1;
            LFPnorm = -LFPnorm;
        end
        disp(['Number of peaks detected: ' num2str(length(locs_Rwave))]);
    
        % A few pre-requisites before the SVD ECG suppression method is run
            % If there are no R-peaks detected at all:
        if isempty(locs_Rwave)
            disp('Inconsistent ECG artifact: LFP signal has no obvious R-peaks')
            % If heartbeats are missed: no ECG artifact is detected:
        elseif any(diff(locs_Rwave)/Fs > 3)
            disp('Inconsistent ECG artifact: LFP signal has no R-peak every three seconds')
            % If heart rate is too low:
        elseif length(locs_Rwave) < (40/60)*(length(LFPecg)/Fs)
            disp('Inconsistent ECG artifact: LFP signal has unreliable heart rate (< 40 bpm)')
        % Else: continue with removal
        else
            disp('Consistent ECG artifact')
        end
    else

        % disp('Pre-defined ECG peaks used');
        LFPnorm = sigin;
        locs_Rwave = ecg_peak_indices;
    end
    
    % Make sure locs_Rwave is an horizontal array for manual peak selection 
    if size(locs_Rwave,1)>size(locs_Rwave,2)
        locs_Rwave = locs_Rwave';
    end
    
    % If manual adjustments of detected peaks is enabled:
    if interactive
    
        % Flip signal according to polarity for manual R-peak detection:
        LFPnorm = LFPnorm * polarity;
    
        % Plot LFP signal (with ECG if present) and detected R-peaks
        fig = figure('Name','Detected R-peaks','NumberTitle','off');clf();
        fig.WindowState = 'maximized';
        if plotecg
            plot(ecg,'g', 'LineWidth', 1.5); hold on
            if nargin>4
                plot(ecg_pks, ecg(ecg_pks),'b*', 'markersize', 15)
            end
        end

        plot(time,LFPnorm,  'k'), hold on,plot(time(locs_Rwave), LFPnorm(locs_Rwave),'r*', 'markersize', 15)
        if plotecg
            legend('ECG','LFPs', ['Detected R-peaks: ' num2str(length(locs_Rwave))], 'FontSize', 11)
        else
            legend('LFPs', ['Detected R-peaks: ' num2str(length(locs_Rwave))], 'FontSize', 11)
        end

        legend('boxoff')
        title('Do you want to add or remove peaks y/n?', 'FontSize', 14);
        xlabel('Time [s]','FontSize', 12)
        box off
        key = input('Inspect detected R-peaks: Press "y" to make manual adjustments, press "n" to continue filtering:','s');
        close('Detected R-peaks')
    
        % Loop over zoomed-in epochs to (de)select R-peaks manually:
        if key == 'y'
            figure('Name','Manual (de)select R-peaks','NumberTitle','off')
            fh = findobj( 'Type', 'Figure', 'Name', 'Manual (de)select R-peaks' ); clf(fh);
            if locs_Rwave == 1
                l_epoch = 10000;
            else
                l_epoch=2000; %1500;%10000;%2500;%5000 %1000 % adjust this value for bigger epoch figures to select R-peaks
            end

            nint=length(LFPnorm)/(l_epoch);
            j = 1;
            while j < nint+1
                fh.WindowState = 'maximized';
                plot(LFPnorm,'k'),hold on,plot(locs_Rwave, LFPnorm(locs_Rwave),'r*')
                xlim([(j-1)*l_epoch j*l_epoch])
                checkpeaks=1;
                while checkpeaks
                    figure(fh)
                    if plotecg
                        plot(ecg,'g'); hold on
                        if nargin>4
                            plot(ecg_pks, ecg(ecg_pks),'b*', 'markersize', 15)
                        end
                    end

                    plot(LFPnorm, 'k'),hold on,plot(locs_Rwave, LFPnorm(locs_Rwave),'r*')
                    title({['Epoch ' num2str(j) '/' num2str(ceil(nint))]; 'Press "r" to remove a peak, "a" to add a peak, "b" to go back one epoch, "f" to skip this step';'Press any other key to continue to next epoch'});
                    xlim([(j-1)*l_epoch j*l_epoch])
                    xlabel('Samples', 'FontSize', 12)
                    box off

                    key=input('Press "r" to remove a peak, "a" to add a peak, "b" to go back one epoch, "f" to skip this step.\nPress any other key to continue to next epoch:','s');
                    if key=='r' % remove peak on click
                        [X,~] = ginput(1);
                        [~,ix] = min(abs(locs_Rwave-X));
                        locs_Rwave(ix) = [];
                        clf

                    elseif key=='a' % add peak on click
                        [X,~] = ginput(1);
                        mn = round(X)-20;%round(X)-5; % adjust this value to select R-peaks more specifically
                        mx = round(X)+20;%round(X)+5;
                        if mx > length(LFPnorm)
                            mx = length(LFPnorm);
                        elseif mn <= 0
                            mn = 1;
                        end

                        x1 = find(LFPnorm == max(LFPnorm(mn:mx)));
                        if length(x1) > 1
                            x2 = find(x1>mn & x1<mx);
                            x1 = x1(max(x2));
                        end

                        locs_Rwave = [locs_Rwave, x1];
                        locs_Rwave = sort(locs_Rwave);

                    elseif key == 'b' % go back one epoch
                        j = j - 1;
                    elseif key == 'f' % finish interactive
                        j = nint+1;
                    % elseif key == 'x'
                    %     j = j + 10;
                    % elseif key=='t'
                    %     j = j - 10;
                    % elseif key == 'c'
                    %     j = j + 100;
                    % elseif key == 'h'
                    %     j = j - 100;
                    else
                        checkpeaks=0; % done with this epoch, go to the next
                        j = j + 1;
                    end
                end
            end
        close('Manual (de)select R-peaks')
        end
        j = 1;
        disp(['Number of peaks detected after manual adjustments: ' num2str(length(locs_Rwave))]);
    end
    
    % Save indices of R-peaks in art_indices:
    art_indicies = locs_Rwave; 
    artifact_count = length(locs_Rwave);
    
    % Average of artifacts
    cnt = 1;
    toremove = [];

    % Create epochs around R-peaks
    for i = 1:artifact_count
        start_index = art_indicies(i) - art_pts_b4_peak; end_index = art_indicies(i) + art_pts_plus_t_after_peak;
        if start_index > 0 && end_index <= length(sigin)

            % Remove R-peak if epoch contains NaNs:
            if any(isnan(sigin(start_index:end_index)))
                disp('Warning: window around peak contains NaN values which cannot be used in SVD')
                toremove = [toremove i];
                NaNpeak = art_indicies(i);
            else
                % Save epochs around R-peaks in allqrst:
                allqrst(cnt,:) = sigin(start_index:end_index);
                if plotecg
                    allqrst_ecg(cnt,:) = ecg(start_index:end_index);
                end
                cnt = cnt+1;
            end

        elseif start_index <= 0
            disp('Warning: first peak near start of the time window cannot be used in SVD')
            toremove = [toremove i];
            startedgepeak = art_indicies(i);
    
        elseif end_index > length(sigin)
            disp('Warning: last peak near end of the time window cannot be used in SVD')
            toremove = [toremove i];
            endedgepeak = art_indicies(i);
        end
    end
    
    art_indicies(toremove) = [];
    
    % Check to make sure at least one R-peak is left to perform the SVD:
    if exist('allqrst','var')
        disp(['Final number of peaks = ',num2str(size(allqrst,1))])

        if showfigs
            fig = figure('Name','Final R-peaks','NumberTitle','off');
            fig.WindowState = 'maximized';
            set(gcf,'color','w');
            if plotecg
                plot(time,ecg,'g'); hold on
            end

            siginplot = sigin  * polarity;
            plot(time,siginplot, 'k'),hold on,plot(time(art_indicies),siginplot(art_indicies),'r*')
            if ~isempty(startedgepeak)
                plot(time(startedgepeak),siginplot(startedgepeak),'m*');
            end
            if ~isempty(NaNpeak)
                plot(time(NaNpeak),siginplot(NaNpeak),'m*');
            end
            if ~isempty(endedgepeak)
                plot(time(endedgepeak),siginplot(endedgepeak),'m*');
            end

            xlabel('Time [s]','FontSize', 12);
            title('Final R-peaks used for SVD', 'FontSize', 14)
            box off
            if savefigs
                print(gcf,'-djpeg',[settings.Folder,filesep,[label,'_detected_R_peaks.jpg']]);
                savefig([settings.Folder,filesep,[label,'_detected_R_peaks']]);
            end
        end

        qrst = mean(allqrst,1);
        time_qrst = (1:length(qrst))/Fs;
    
        % Initiate output signal
        sigout = sigin;
    
        % Perform SVD
        [U,S,V] = svd(allqrst');
        var_ex=diag(S).^2/sum(diag(S).^2);
        if ncomp > 20 % take explained variance
            ind = find(cumsum(var_ex) >= ncomp/100);
            ncomp = ind(1);
        end
    
        if showfigs
            if length(var_ex)>9
                fig = figure('Name','SVD components','NumberTitle','off'); 
                set(gcf,'color','w','position', [ 43         150        1700         720]);
                fig.WindowState = 'maximized';
                sgtitle(figname, 'FontSize', 18, 'Interpreter', 'none')
                plt = [1 2 4 5 7 8];
                clrs1 = [[65 182 196]/256;[29 145 192]/256;[34 94 168]/256;[12 44 132]/256;[0 24 46]/256;[0 12 23]/256];
                labels = [];

                for num = 1:length(plt)
                    subplot(3,3,[plt(num)])
                    for num2 = 1:num
                        plot(time_qrst,U(:,num2),'linewidth',2, 'Color', clrs1(num2,:)); title(['SVD' num2str(num) ' - explained % energy = ',num2str(100*sum(var_ex(1:num))),' %'], 'Interpreter', 'none');
                        hold on;
                    end
                    labels = [labels;num2str(num)];
                    legend(labels, 'Location', 'southeast')
                    if plt(num) < 6
                        hAx = gca;
                        hAx.Position = hAx.Position.*[1 0.92 1 1];
                    else
                        xlabel('Time [s]','FontSize', 12);
                    end
                    ylim([min(min(U(:,1:4)))*1.05 max(max(U(:,1:4)))*1.05])
                    yticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6])
                    set(gca,'fontsize',10)
                end

                subplot(3,3,[3 6 9]),bar(var_ex(1:10)*100,'k');title('SVD components - explained % energy');
                hAx = gca;
                hAx.OuterPosition = hAx.OuterPosition.*[1 1 1 0.95];
                set(gca,'fontsize',10)
                xlabel('Component number')
            end
        end
    
        % Select number of components based on visual inspection:
        if settings.nCompFinal
            prs = input('Indicate which components have an ECG pattern (write as e.g.: "[1, 2, 4]"): ','s');
            ncomp = str2num(prs);
            x = 0;
    
            while x == 0
                try U(:,ncomp);
                    x = 1;
                catch
                    prs = input('Incorrect input. Write as e.g.: "[1, 2, 4]": ','s');
                    ncomp = str2num(prs);
                end
            end 
    
        else
            ncomp = 1:ncomp;
            disp(['SVD components chosen: ' num2str(ncomp)]);
        end
    
        if showfigs
            if savefigs
                saveas(fig, [settings.Folder,filesep,[label,'_SVD_components.jpg']]);
                saveas(fig, [settings.Folder,filesep,[label,'_SVD_components.fig']]);
            end
            close('Final R-peaks')
            close('SVD components')
        end
    
        % Save ECG projection per R-peak in "proj":
        proj = U(:,ncomp)*S(ncomp,:)*V';
    
        % Remove artifacts from data
        % Loop over each R-peak
        for i = 1:length(art_indicies)
    
            % Save original ECG projection per R-peak
            columnMeans = proj(:,i);
            columnMeans1 = columnMeans;
    
            % Avoid the introduction of artifact by correcting sudden jumps at start and end of each epoch
            % The first and last 5 ms of the template is searched for samples with the smallest difference
            if length(columnMeans1) >= round(0.05*Fs)
                Qwave = round(0.04*Fs);
                Swave = length(columnMeans1)-round(0.04*Fs);
            else
                Qwave = round((length(columnMeans1)/2)-2);
                Swave = length(columnMeans1)-round((length(columnMeans1)/2)-2);
            end
            firsthalf = columnMeans1(1:Qwave);
            secondhalf = columnMeans1(Swave:end);
            [firstzero, secondzero] = findSmallestDifference(firsthalf, secondhalf);
            if isnan(firstzero) % In this case, artifact corrupted the signal
                continue
            end

            % In case two samples of the same tail have the same value, choose samples to obtain
            % longest template (so firstzero(1) and secondzero(end)):
            if length(firstzero)>1
                firstzero = firstzero(1);
            end
            secondzero = (Swave-1) + secondzero;
            if length(secondzero)>1
                secondzero = secondzero(end);
            end
            % Let tails run smoothly to smallest value
            columnMeans(1:firstzero-1) = columnMeans1(firstzero);       
            columnMeans(secondzero+1:end) = columnMeans1(secondzero);
    
            % Important that both tails have same value:
            if columnMeans(end) < columnMeans(1)
                if firstzero > 1
                    columnMeans(1:firstzero-1) = columnMeans(end);
                elseif firstzero == 1
                    columnMeans(firstzero) = columnMeans(end);
                end
            elseif columnMeans(end) > columnMeans(1)
                if secondzero < length(columnMeans)
                    columnMeans(secondzero+1:end) = columnMeans(1);
                elseif secondzero == length(columnMeans)
                    columnMeans(secondzero) = columnMeans(1);
                end
            end
    
            % Optimisation of the ECG projection of each R-peak epoch with a linear least squares problem:
            start_index(i) = art_indicies(i) - art_pts_b4_peak;
            end_index(i) = art_indicies(i) + art_pts_plus_t_after_peak;
            y = sigout(start_index(i):end_index(i));
    
            % Initial Guess & Upper & Lower Bounds
            p0 = 2;
            lb = -10;
            ub = 10;
    
            % errorfunc
            func = @(p)fune(p,columnMeans,y);
    
            % b=mean(yu)
            parest = min(max(mean(y - columnMeans), lb), ub);
    
            % Estimate & Plot yest
            [sigout(start_index(i):end_index(i)), ymod2(i,:)] = fune(parest,columnMeans,y);
    
            % Save the projection of the ECG artifact of each R-peak epoch 
            proj_out(:,i) = ymod2(i,:);
            clear y columnMeans columnMeans1 tri
        end
    
        % Remove zero-projections - artifact
        proj_out( :, ~any(proj_out,1) ) = [];
    
        % Compare spectral profiles of input and output LFP signal
        if showfigs
    
            % Pwelch for plotting PSDs cannot be applied to NaNs, so this finds the start
            % and end indices of segments of the Streaming data without missing values
            i_seg = [1; find(diff(any(isnan(sigin(:,1)),2)))+1; length(sigin)];
    
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
                LFP_original{s} = sigin(ix);
                LFP_cleaned{s} = sigout(ix);
    
            end
    
            % Plot LFP signal and detected R-peaks
            fig = figure('Name','Original and ECG suppressed LFP signal','NumberTitle','off', 'visible', 'off');clf();
            fig.WindowState = 'maximized';
            sgtitle(figname, 'Interpreter', 'none')
            subplot(2,s,1:s)
            plot(time,sigin,  'r'),hold on, plot(time,sigout,'k')
            xlabel('Time [s]', 'FontSize',12)
            box off
            legend('Original LFP', 'ECG suppressed LFP', 'FontSize',11, 'location', 'northwest')
            legend('boxoff')
    
            % Plot PSD of original and ECG suppressed (segments of) LFP signal
            for k = 1:s
                if length(LFP_original{k}) > Fs-1
                    [P_in,~]=pwelch(LFP_original{k},[],[],[],Fs);
                    [P_out,F]=pwelch(LFP_cleaned{k},[],[],[],Fs);
                    subplot(2,s,s+k)
                    plot(F,P_in,'r'),hold on, plot(F,P_out,'k');
                    ylim([0 40]),xlim([0 130]);
                    xlabel('Frequency [Hz]', 'FontSize', 12)
                    if k == 1
                        ylabel('Power', 'FontSize', 12)
                    end
                    title(['Segment ' num2str(k)], 'Interpreter', 'none', 'FontSize', 14);
                    legend('Original','ECG suppressed', 'FontSize', 11)
                    box off
                else
                    subplot(2,s,s+k)
                    title(['Segment ' num2str(k) ' < 250 samples'], 'FontSize', 14')
                end
            end

            fig.Visible = 'on';
            if savefigs
                savefig(fig, [settings.Folder,filesep,[label,'_postcleaning_',num2str(ncomp),'comp.fig']]);
                print(fig,'-djpeg',[settings.Folder,filesep,[label,'_postcleaning_',num2str(ncomp),'comp.jpg']]);
            end
        end
    else
        allqrst = [];
        disp(['Final number of peaks = ',num2str(size(allqrst,1))])
    
        sigout = sigin;
        proj_out = [];
    end
    
    % Save only R-peak indices used for ECG suppression:
    ecg_peak_indices = sort(art_indicies);

    % Save selection of SVD components:
    finalcomp = ncomp;
end
    
% ---- Below are two functions used in this script:
    
function [idxA, idxB, result] = findSmallestDifference(A, B)
    
    % This function finds the indices and value of one sample of A and one 
    % sample of B that are closest to each other.
    
    % Sort both arrays using sort function
    A1 = sort(A);
    B1 = sort(B);
    
    m = length(A);
    n = length(B);
    
    % Initialize result as max value
    result = 70;
    
    % Scan Both Arrays up to size of of the arrays
    for a = 1:m
        for b = 1:n
            if (abs(A1(a) - B1(b)) < result)
                result = abs(A1(a) - B1(b));
                a1 = a;
                b1 = b;
            end
        end
    end
    if exist('a1','var')
        idxA = find(A == A1(a1));
        idxB = find(B == B1(b1));
    else
        idxA = NaN;
        idxB = NaN;
    end
    result;
end
    
% ----
    
function [e, yhat] = fune(p,u,y)
    
    % Error Function used for linear least square method
    b = p(1);
    
    % Define estimated y: yhat
    yhat = u + b;
    
    % Define error: e
    e = y - yhat;
    % J = e' * e;
end
