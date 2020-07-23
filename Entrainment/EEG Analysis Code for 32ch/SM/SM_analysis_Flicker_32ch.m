%open file and plot data EEG preprocessed format 



clear

input_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG Testing\CABI\EEG';
input_file = 'MGtestFlickerAV_f_1_200_butterord3_32ch'; %without .mat extension
file_name = 'MGtstoutput';
%%%comment/uncomment sections below to choose montage for analysis

% montage = 'common'; %'average', 'common', 'csd'
% channels = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
% common_montage_ch = 48; %channel used for common montage (48 is Cz)
% montage_plot_order = [0 0 1 0 33 0 34 0 0    2 0 3 0 37 0 36 0 35    7 6 5 4 38 39 40 41 42    8 9 10 11 47 46 45 44 43   15 14 13 12 48 49 50 51 52   16 17 18 19 32 56 55 54 53   23 22 21 20 31 57 58 59 60   25 0 26 0 30 0 63 0 62    0 0 27 0 29 0 64 0 0];
% %channels = {'F3' 'F4' 'C3' 'C4' 'O1' 'O2'}; %channels that will be analyzed 

montage = 'bipolar';
channels = {'Fpz_AFz' 'AFz_Fz' 'Fz_FCz' 'FCz_Cz' 'Fpz_F2' 'F2_FC2' 'FC2_C2' 'C2_CP2' 'CP2_P2' 'Fp2_AF4' 'AF4_F4' 'F4_FC4' 'FC4_C4' 'C4_CP4' 'CP4_P4' 'P4_PO4' 'PO4_O2' 'Fp2_F6' 'F6_FC6' ...
    'FC6_C6' 'C6_CP6' 'CP6_P6' 'P6_O2' 'Fp2_AF8' 'AF8_F8' 'F8_FT8' 'FT8_T8' 'T8_TP8' 'TP8_P8' 'P8_PO8' 'PO8_O2' 'F1_FC1' 'FC1_C1' 'C1_CP1' 'CP1_P1' 'Fp1_AF3' 'AF3_F3' 'F3_FC3' 'FC3_C3' 'C3_CP3' 'CP3_P3' 'P3_PO3' 'PO3_O1' 'Fp1_F5' 'F5_FC5' ...
    'FC5_C5' 'C5_CP5' 'CP5_P5' 'P5_O1' 'Fp1_AF7' 'AF7_F7' 'F7_FT7' 'FT7_T7' 'T7_TP7' 'TP7_P7' 'P7_PO7' 'PO7_O1' 'P1_Oz' 'CPz_Pz' 'Pz_POz' 'POz_Oz' 'Fpz_F1' 'Cz_CPz'  'P2_POz'};  
montage_plot_order = [50 0 36 0 1 0 10 0 24    51 44 37 62 2 5 11 18 25    52 45 38 32 3 6 12 19 26    53 46 39 33 4 7 13 20 27       54 47 40 34 63 8 14 21 28    55 48 41 35 59 9 15 22 29     56 49 42 58 60 64 16 23 30    57 0 43 0 61 0 17 0 31];
%channels = {'FC3_C3', 'FC4_C4', 'PO3_O1', 'PO4_O2'};


analysis_start = 0; %in seconds (set to 0 for beginning of file); applied after artifact rejection
analysis_end = -1; %in seconds (set to -1 for end of file); applied after artifact rejection
remove_artifact = 1; %1 to cut out artifact data; 0 not to


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data file
load([input_dir '\' input_file '.mat']);

% disp('PLOTTING ALL CHANNELS!!!')
% unused_eeg_channels = []; %tmp to see all channels

%assuming 512Hz sampling freq
WINDOW = 512;
NOOVERLAP = 256;
NFFT = 512;

%convert data labels in 'channels' to numerical indices. These are channels that will be analyzed
%if channel labels are bipolar (meaning they contain a dash), then channel_ind will be a 2-column matrix
channels_ind = [];  
for i = 1:length(channels)
    flag = 0; ch_label = channels{i};
    if (findstr(ch_label, '_'))  %bipolar montage (channel label contains underscore - DO NOT USE DASH BECAUSE IT CAN BE INTERPRETED AS A MINUS BY MATLAB)
        dash_index = strfind(ch_label, '_'); 
        ch1_label = ch_label(1:dash_index-1); flag1 = 0; ch1_ind = 0; 
        ch2_label = ch_label(dash_index+1:end); flag2 = 0; ch2_ind = 0; 
        for j = 1:length(eeg_channel_labels)
            if strcmp(eeg_channel_labels(j), ch1_label)
                ch1_ind = j; flag1 = 1; 
            end
            if strcmp(eeg_channel_labels(j), ch2_label)
                ch2_ind = j; flag2 = 1; 
            end
        end 
        if (flag1 == 1 && flag2 == 1)
            channels_ind = [channels_ind; [ch1_ind ch2_ind]];
            flag = 1;
        end    
    else    
        for j = 1:length(eeg_channel_labels)
            if strcmp(eeg_channel_labels(j), ch_label)
                channels_ind = [channels_ind j]; flag = 1; break
            end
        end 
    end    
    if (flag == 0)
        fprintf('WARNING: Could not find channel %s in eeg_channel_labels. Channel indexing will be inaccurate.\n', ch_label)
    end    
end 
if (isempty(channels_ind))
    fprintf('ERROR: Could not find any channels for analysis...exiting\n')
    return; 
end    

fprintf('File: %s ... original file is %f seconds long, sampling rate for analysis is %f Hz\n', [input_file '.mat'], size(eeg_data,2)/samp_freq, samp_freq)
fprintf('Unused channels: %d\n', length(unused_eeg_channels))
if (~isempty(unused_eeg_channels))
    eeg_channel_labels(unused_eeg_channels)
end    

%check if any artifact times, and if so, cut them out
if (~exist('rej_time', 'var') || isempty(rej_time))
    warning('No artifact rejections in this file')
end
if (exist('rej_time', 'var') && ~isempty(rej_time) && remove_artifact == 1)
        eeg_rejects = ones(1,size(eeg_data,2));
        for i = 1:size(rej_time,1) 
            eeg_rejects(1,floor(rej_time(i,1)*samp_freq)+1:floor(rej_time(i,2)*samp_freq)) = 0;  
        end    
        eeg_rejects = logical(eeg_rejects); 
        seconds_data_removed = (size(eeg_data,2) - sum(eeg_rejects))/samp_freq; 
        fprintf('%d points or %.2f seconds of data (%.1f percent) was manually cut from eeg\n', size(eeg_data,2) - sum(eeg_rejects), (size(eeg_data,2) - sum(eeg_rejects))/samp_freq, 100*((size(eeg_data,2) - sum(eeg_rejects))/size(eeg_data,2)));  
        eeg_data = eeg_data(:,eeg_rejects);
   
    if (~isempty(aux_data))
        aux_data = aux_data(:,eeg_rejects);
    end      
end

%convert analysis_start and _end from seconds to points, and select desired time interval
my_time_start = round(analysis_start*samp_freq);
my_time_end = round(analysis_end*samp_freq);
if (analysis_end == -1 || round(analysis_end*samp_freq) > size(eeg_data,2))
    my_time_end = size(eeg_data,2);
end    
if (analysis_start < 0 || my_time_start > size(eeg_data,2) || analysis_end < -1)
    disp(sprintf('ERROR: Impossible analysis_start...exiting'))
    return; 
end    
if (my_time_start == 0)
    my_time_start = 1;
end
eeg_data = eeg_data(:, my_time_start:my_time_end); %cut data to desired analysis length 
aux_data = aux_data(:, my_time_start:my_time_end); %cut data to desired analysis length 

%setup montage 
if (size(channels_ind,2) == 2 && size(channels_ind,1) > 1 && strcmp(montage, 'bipolar') == 0 )
    sprintf('WARNING: channel names contain two channels, but montage NOT set to BIPOLAR. Will ignore second channel')
end    
eeg = [];
if (strcmp(montage, 'none')) %no montage (technically not correct)
    for i = 1:length(channels_ind)
        eeg(i,:) = eeg_data(channels_ind(i),:); 
    end
elseif (strcmp(montage, 'bipolar'))
    if (size(channels_ind,2) ~= 2)
        error('ERROR: Bipolar montage not appropriately defined (needs two channels with underscore in between)...exiting\n')  %expecting channel labels to be 'C2-C3' so that channels_ind is a 2-column matrix
    else
        for i = 1:size(channels_ind,1)
            eeg(i,:) = eeg_data(channels_ind(i,1),:) - eeg_data(channels_ind(i,2),:); 
            if (ismember(channels_ind(i,1), unused_eeg_channels) || ismember(channels_ind(i,2), unused_eeg_channels))  
                eeg(i,:) = zeros(1,length(eeg_data));   %set unused channel to zero
            end    
        end  
    end    
elseif (strcmp(montage, 'average'))  %calculate average signal without unused channels
    avg_signal = mean(eeg_data(setdiff([1:size(eeg_data,1)],[unused_eeg_channels]),:));    
    for i = 1:length(channels_ind)
        eeg(i,:) = eeg_data(channels_ind(i),:) - avg_signal; 
        if (ismember(channels_ind(i), unused_eeg_channels))  
            eeg(i,:) = zeros(1,length(eeg_data));   %set unused channel back to zero (above line will set it to avg_signal)
        end    
    end  
elseif (strcmp(montage, 'common'))  %use single channel to reference (last one in channels list)
    if (ismember(channels_ind(common_montage_ch), unused_eeg_channels))
        fprintf('ERROR: Reference channel is unused channel...exiting\n')
    end    
    ref_signal = (eeg_data(channels_ind(common_montage_ch),:)); 
    for i = 1:length(channels_ind)
        eeg(i,:) = eeg_data(channels_ind(i),:) - ref_signal; 
        if (ismember(channels_ind(i), unused_eeg_channels))  
            eeg(i,:) = zeros(1,length(eeg_data));   %set unused channel back to zero (above line will set it to ref_signal)
        end 
    end      
elseif (strcmp(montage, 'csd')) %current source density
    load('C:\Users\Svjetlana\Box Sync\Research\MATLAB\SCRIPTS\helper_functions\CSD_GHmat_64ch_spline4.mat') %load G and H matrices for spatial laplacian
    eeg_dataCSD = CSD(eeg_data, G, H);  
    for i = 1:length(channels_ind)
        eeg(i,:) = eeg_dataCSD(channels_ind(i),:); 
        if (ismember(channels_ind(i), unused_eeg_channels))  
            eeg(i,:) = zeros(1,length(eeg_data));   %set unused channel back to zero 
        end
    end    
else
    error('Unknown montage')
end


%eeg=ref_data;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = (1/samp_freq)*[1:size(eeg,2)]; % time points

%eegplot(eeg,'srate', samp_freq, 'title', 'EEG - after artifact rejection and start-end analysis times');


%calculate and plot spectrogram 
if (length(channels) <= 10)
 figure;
 for i = 1:length(channels)
	signal = eeg(i,:);   
    [S,F,T,P] = spectrogram(signal,WINDOW,NOOVERLAP,NFFT,samp_freq);
    h = subplot(1,length(channels),i);
    surf(T,F,log10(P),'edgecolor','none');  
    axis tight; shading interp; view(0,90);
    %colormap; colorbar; 
    set(gca, 'YLim', [2 130]);  %freq WAS 
    set(gca, 'CLim', [-2 2]); %WAS [0 3]   [-2 1]
    if (size(channels_ind, 1) == 1)  %single channel name
        title(['Channel ', num2str(channels_ind(i)) ' ' char(eeg_channel_labels(channels_ind(i))), ' Time-freq'], 'FontSize', 7,'interpreter','none')
    else %bipolar channel name
        title(['Channel ', num2str(i) ' ' char(eeg_channel_labels(channels_ind(i,1))) '_' char(eeg_channel_labels(channels_ind(i,2))), ' Time-freq'], 'FontSize', 7,'interpreter','none')
    end    
   % tmpt=round((my_time_end-my_time_start)/samp_freq);
    %set(gca, 'XTick', [0:round(tmpt/10):tmpt], 'XTickLabel', [round(my_time_start/samp_freq):round(tmpt/10):round(my_time_start/samp_freq+tmpt)]);
    xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');
    set(gcf,'color','w');
 end
end

%calculate PSD 
Pxx = [];
for i = 1:length(channels)
	[Pxx(:,i), f] = pwelch(eeg(i,:),WINDOW,NOOVERLAP,NFFT,samp_freq);   
end
Pxx_log = log10(Pxx); 


%use PSD to define bad channels (those where gamma power is higher than alpha/beta power)
%bad channels will have channel name written in red
bad_channels = []; %list of channels where gamma power is too high 
low_band_freq = [8 13]; %start and stop values for defining beta band (in Hz), for calc bad channels  
high_band_freq = [70 90]; %start and stop values for defining gamma band (in Hz), for calc bad channels  
lowband_start_index = find(f >= low_band_freq(1),1);
lowband_stop_index = find(f >= low_band_freq(2),1)-1;
highband_start_index = find(f >= high_band_freq(1),1);
highband_stop_index = find(f >= high_band_freq(2),1)-1;
for i = 1:length(channels)
	avg_low_power = mean(Pxx(lowband_start_index:lowband_stop_index,i),1); 
	avg_high_power = mean(Pxx(highband_start_index:highband_stop_index,i),1);    
    if (avg_high_power > avg_low_power)
    	bad_channels = [bad_channels i];
    end    
end
bad_channels


%calculate amplitude of discrete peak (log), relative to baseline
peak_discrete_freq = 40; %freq at which to calculate peak amplitude
peak_discrete_index = find(f >= peak_discrete_freq,1);
amp_discrete_peak_log = [];
for i = 1:length(channels)
    amp_discrete_peak_log(i) = Pxx_log(peak_discrete_index,i) - (Pxx_log(peak_discrete_index+2,i)+Pxx_log(peak_discrete_index-2,i))/2; %peak height is diff between peak and baseline
    if (amp_discrete_peak_log(i) < 0)
        amp_discrete_peak_log(i) = 0; 
    end    
    if ismember(i, bad_channels)
        amp_discrete_peak_log(i) = 0; 
    end    
end

%calculate amplitude of discrete peak (absolute), 31-39Hz and 41-49Hz bands (no log)
%entrained channels are those where peak map is above mean+3STD of surrounding bands 
peak_discrete_freq = 40; %freq at which to calculate peak amplitude
peak_discrete_index = find(f >= peak_discrete_freq,1);
band1 = [31 39]; 
band2 = [41 49];
band1_start_index = find(f >= band1(1),1); band1_stop_index = find(f >= band1(2),1);
band2_start_index = find(f >= band2(1),1); band2_stop_index = find(f >= band2(2),1);
amp_discrete_peak = []; entrained_ch = [];std_above_mean = [];
for i = 1:length(channels)
    amp_discrete_peak(i) = Pxx(peak_discrete_index,i);   
    band1_mean(i) = mean(Pxx(band1_start_index:band1_stop_index,i));
    band2_mean(i) = mean(Pxx(band2_start_index:band2_stop_index,i));
    band1_std(i) = std(Pxx(band1_start_index:band1_stop_index,i));
    band2_std(i) = std(Pxx(band2_start_index:band2_stop_index,i));
    if (amp_discrete_peak(i) > band1_mean(i)+3*band1_std(i) && amp_discrete_peak(i) > band2_mean(i)+3*band2_std(i))
        entrained_ch = [entrained_ch i];
    end
    
    std_above_mean(i) = ...
    (amp_discrete_peak(i) - mean([band1_mean(i),band2_mean(i)]))/mean([band1_std(i),band2_std(i)]); 
end


fprintf('Number of channels entrained at %d Hz is %d (out of %d total)\n', peak_discrete_freq, length(entrained_ch), length(channels));


%plot log PSD
% if (length(channels) <= 10)
%     figure()
%     for i = 1:length(channels)
%         h2 = subplot(1,length(channels),i); 
%         plot(f,Pxx_log(:,i))  
%         xlim([2,100]); ylim([-2,2]) 
%         set(gcf,'color','w');
%         
%         if (ismember(i, bad_channels))
%             if (size(channels_ind, 1) == 1)  %single channel name
%                 title(['Log PSD Ch ', num2str(channels_ind(i)) ' ' char(eeg_channel_labels(channels_ind(i)))], 'FontSize', 7, 'color', 'r','interpreter','none')
%             else %bipolar channel name
%                 title(['Log PSD Ch ', num2str(i) ' ' char(eeg_channel_labels(channels_ind(i,1))) '_' char(eeg_channel_labels(channels_ind(i,2)))], 'FontSize', 7, 'color', 'r', 'Interpreter','none')
%             end  
%         else
%             if (size(channels_ind, 1) == 1)  %single channel name
%                 title(['Log PSD Ch ', num2str(channels_ind(i)) ' ' char(eeg_channel_labels(channels_ind(i)))], 'FontSize', 7,'interpreter','none')
%             else %bipolar channel name
%                 title(['Log PSD Ch ', num2str(i) ' ' char(eeg_channel_labels(channels_ind(i,1))) '_' char(eeg_channel_labels(channels_ind(i,2)))], 'FontSize', 7,'interpreter','none')
%             end              
%         end    
%     end
% else
%     figure() ;  clf ; set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] ) ;
%     nCol = 9;  nRow = round(length(montage_plot_order)/nCol);  %ASSUMES 9 ROWS OF PLOTS
%     rowH = 0.58 / nRow ;  colW = 0.7 / nCol ; colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ; rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
%     ct = 0; 
%     for c = montage_plot_order
%         ct=ct+1; 
%         rowId = ceil( ct / nCol ) ;
%         colId = ct - (rowId - 1) * nCol ;
%         axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
%         if (c == 0)
%             plot(1,1); continue
%         end
%         plot(f,Pxx_log(:,c),'b-'); %ylim([-2.5, -1])
%         xlim([2 100]);  ylim([-3,2]) 
% 
%         if (ismember(c, bad_channels))
%            if (size(channels_ind, 1) == 1)  %single channel name
%                 title(['Ch ', num2str(channels_ind(c)) ' ' char(eeg_channel_labels(channels_ind(c)))], 'FontSize', 7, 'color', 'r','interpreter','none')
%            else %bipolar channel name
%                 title(['Ch ', num2str(c) ' ' char(eeg_channel_labels(channels_ind(c,1))) '_' char(eeg_channel_labels(channels_ind(c,2)))], 'FontSize', 7, 'color', 'r','interpreter','none')
%            end  
%         else
%            if (size(channels_ind, 1) == 1)  %single channel name
%                 title(['Ch ', num2str(channels_ind(c)) ' ' char(eeg_channel_labels(channels_ind(c)))], 'FontSize', 7,'interpreter','none')
%            else %bipolar channel name
%                 title(['Ch ', num2str(c) ' ' char(eeg_channel_labels(channels_ind(c,1))) '_' char(eeg_channel_labels(channels_ind(c,2)))], 'FontSize', 7,'interpreter','none')
%            end  
%         end          
%     end
%     axes('Position', [0, 0.95, 1, 0.05]) ; set(gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White') ;
%     text(0.5, 0, 'Log PSD', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','interpreter','none') ;
% end    



% %plot other data
% if (length(channels) <= 10)
%  figure()
%  for i = 1:length(channels)
%     subplot(length(channels),1,i)
%     plot(T, eeg(i,:));
%     title('EEG channels (entire signal)')
%  end
% else
%     figure() ;  clf ; set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] ) ;
%     nCol = 9;  nRow = round(length(montage_plot_order)/nCol);
%     rowH = 0.58 / nRow ;  colW = 0.7 / nCol ; colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ; rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
%     ct = 0; 
%     for c = montage_plot_order
%         ct=ct+1; 
%         rowId = ceil( ct / nCol ) ;
%         colId = ct - (rowId - 1) * nCol ;
%         axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
%         if (c == 0)
%             plot(1,1); continue
%         end
%         plot(T,eeg(c,:),'b-'); axis tight;
%         if (ismember(c, bad_channels))
%             if (size(channels_ind, 1) == 1)  %single channel name
%                 title(['Ch ', num2str(channels_ind(c)) ' ' char(eeg_channel_labels(channels_ind(c)))], 'FontSize', 7, 'color', 'r','interpreter','none')
%             else %bipolar channel name
%                 title(['Ch ', num2str(c) ' ' char(eeg_channel_labels(channels_ind(c,1))) '_' char(eeg_channel_labels(channels_ind(c,2)))], 'FontSize', 7, 'color', 'r','interpreter','none')
%             end  
%         else
%             if (size(channels_ind, 1) == 1)  %single channel name
%                 title(['Ch ', num2str(channels_ind(c)) ' ' char(eeg_channel_labels(channels_ind(c)))], 'FontSize', 7,'interpreter','none')
%             else %bipolar channel name
%                 title(['Ch ', num2str(c) ' ' char(eeg_channel_labels(channels_ind(c,1))) '_' char(eeg_channel_labels(channels_ind(c,2)))], 'FontSize', 7,'interpreter','none')
%             end  
%         end    
%     end
%     axes('Position', [0, 0.95, 1, 0.05]) ; set(gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White') ;
%     text(0.5, 0, 'EEG ', 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','interpreter','none') ;
% end
% 
% 
% if (exist('aux_data','var') && ~isempty(aux_data))  
%     figure; 
%     for i = 1:size(aux_data,1)
%         subplot(size(aux_data,1),1,i)
%         plot(T, aux_data(i,:))
%         try
%             title(['AUX ' num2str(i) '-' char(aux_labels(i))])
%         end    
%     end    
% end

% if (length(channels) == 64)
%     my_limits = [min(amp_discrete_peak_log) max(amp_discrete_peak_log)]; 
%     if (strcmp(montage, 'bipolar'))
%         figure; topoplot(amp_discrete_peak_log, 'E:\EEG_Data\Standard-10-20-Cap64.locs', 'maplimits', my_limits, 'style', 'map');
%         colorbar; title(['Max log value for 40Hz peak ' num2str(max(amp_discrete_peak_log))])
%     else
%         figure; topoplot(amp_discrete_peak_log,'E:\EEG_Data\Standard-10-20-Cap64.locs', 'maplimits', my_limits); 
%         colorbar; title(['Max log value for 40Hz peak ' num2str(max(amp_discrete_peak_log))])
%     end    
% end


%calculate average power in various frequency bands
%normalization is optional (dividing entire PSD by power in normalization band)
% delta_band_freq = [2 4]; %start and stop values for defining frequency band
% theta_band_freq = [4 8];  
% alpha_band_freq = [8 12]; 
% beta_band_freq = [13 30]; 
% lobeta_band_freq = [13 20];
% hibeta_band_freq = [20 30];
% logamma_band_freq = [30 50];
% gamma_band_freq = [50 200];
% norm_band_freq = [3 50]; %freq band for normalization 
% 
% delta_start_index = find(f >= delta_band_freq(1),1); delta_stop_index = find(f >= delta_band_freq(2),1);
% theta_start_index = find(f >= theta_band_freq(1),1); theta_stop_index = find(f >= theta_band_freq(2),1);
% alpha_start_index = find(f >= alpha_band_freq(1),1); alpha_stop_index = find(f >= alpha_band_freq(2),1);
% beta_start_index = find(f >= beta_band_freq(1),1); beta_stop_index = find(f >= beta_band_freq(2),1);
% lobeta_start_index = find(f >= lobeta_band_freq(1),1); lobeta_stop_index = find(f >= lobeta_band_freq(2),1);
% hibeta_start_index = find(f >= hibeta_band_freq(1),1); hibeta_stop_index = find(f >= hibeta_band_freq(2),1);
% logamma_start_index = find(f >= logamma_band_freq(1),1); logamma_stop_index = find(f >= logamma_band_freq(2),1);
% gamma_start_index = find(f >= gamma_band_freq(1),1); gamma_stop_index = find(f >= gamma_band_freq(2),1);
% norm_start_index = find(f >= norm_band_freq(1),1); norm_stop_index = find(f >= norm_band_freq(2),1);
% 
% % %normalize psd before calculating average power
% norm_flag = 0; %1 to normalize, 0 not to 
% if (norm_flag == 1)
%     norm_power = mean(Pxx(norm_start_index:norm_stop_index,:),1); 
%     Pxx2 = Pxx./repmat(norm_power, size(Pxx,1),1);
% 
%     figure
%     plot(f, log10(Pxx(:,1)),'b-')
%     hold on
%     plot(f, log10(Pxx2(:,1)), 'r-'); 
%     xlim([0 100])
%     legend({'Original Log PSD', 'Normalized Log PSD'})
% else
%     Pxx2 = Pxx; 
% end
% 
% delta_power = mean(Pxx2(delta_start_index:delta_stop_index,:),1); 
% theta_power = mean(Pxx2(theta_start_index:theta_stop_index,:),1); 
% alpha_power = mean(Pxx2(alpha_start_index:alpha_stop_index,:),1); 
% beta_power = mean(Pxx2(beta_start_index:beta_stop_index,:),1); 
% lobeta_power = mean(Pxx2(lobeta_start_index:lobeta_stop_index,:),1); 
% hibeta_power = mean(Pxx2(hibeta_start_index:hibeta_stop_index,:),1); 
% logamma_power = mean(Pxx2(logamma_start_index:logamma_stop_index,:),1); 
% gamma_power = mean(Pxx2(gamma_start_index:gamma_stop_index,:),1); 

% fprintf('Channel\t Delta\t Theta\t Alpha\t Beta\t LoBeta\t Hibeta\t LoGamma\t Gamma\n')
% for c = 1:length(channels)
%     fprintf('%s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n', char(eeg_channel_labels(channels_ind(c))), delta_power(c), theta_power(c), alpha_power(c), beta_power(c), ...
%                                                                                  lobeta_power(c), hibeta_power(c), logamma_power(c), gamma_power(c))
% end

save([input_dir '\' file_name '.mat'], ...
'eeg', 'aux_data', 'samp_freq', 'unused_eeg_channels', 'eeg_channel_labels', ...
'aux_labels', 'rej_time', 'event_stop_times_pts', 'event_stop_times','std_above_mean' );

