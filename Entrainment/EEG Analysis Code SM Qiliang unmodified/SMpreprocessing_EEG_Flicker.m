%open .bdf EEG file, filter, organize channels and labels, and save in .mat format
clear

%user defined variables
input_dir = 'C:\Users\mattokaren3\Box\Project_FlickerBloodFlow\EEG\Emory BHC Recordings\DO NOT USE FOR FLICKER Subject MA\DO NOT USE FOR FLICKER Entrainment Testing 01OCT2019';
input_file = 'Subject MA TEST 24.4 V eyes open with 15 seconds post stim';  %without .bdf extension
output_dir = input_dir;
cap_locs = 'C:\Users\mattokaren3\Box\Project_FlickerBloodFlow\EEG Testing\EEG Cap\Standard-10-20-Cap64.locs';

% for config file 64+8
eeg_channels = [1:64];  
aux_channels = [73 74]; %auxiliary channels [73 74]
aux_labels = {'photodiode', 'microphone'};  %{'photodiode', 'microphone'}
unused_eeg_channels = []; %unused EEG channels indices (they may have been disconnected OR they are too noisy); vector; data in those channels will be set to 0

eeg_locutoff = 1;     %low cutoff frequency for EEG filtering, set to -1 if not using; 1 for all
eeg_hicutoff = 200;    %high cutoff frequency for EEG filtering, set to -1 if not using; set to 1000 for 2kHz files, 8000 for 8kHz files
eeg_filter_order = 3;  %filter order for eeg processing (both high/lo and notch if applicable); 4000 for fir1, 3 for butter
eeg_filter_type = 'butter'; %'fir1' or 'butter' (both high/lo and notch if applicable)
aux_locutoff = -1;     %low cutoff frequency for aux filtering, set to -1 if not using; 1 for all
aux_hicutoff = -1;    %high cutoff frequency for aux filtering, set to -1 if not using; set to 1000 for 2kHz files, 8000 for 8kHz files
aux_filter_order = 3; %filter order for aux processing (both high/lo and notch if applicable); 4000 for fir1, 3 for butter
aux_filter_type = 'butter'; %'fir1' or 'butter' (both high/lo and notch if applicable)

notch_freq_vector = []; %frequencies where to apply notch filter; leave blank if notch filtering not desired
notch_filter_on = 0; %1 to apply notch filter to saved preprocessed data, 0 not to apply. If set to 0, notch filter will be applied only to data used for artifact rejection.

disp_channel_labels = 1; % 1 to display channel labels in eegplot; 0 to display channel number

%%%END OF USER DEFINED VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read in EEG file
time_range = [1 120]; % only take the first 2 mins
down_fs = 512; %downsample frequency to 
EEG = pop_biosig([input_dir '\' input_file '.bdf'],'blockrange',time_range);
EEG = pop_resample( EEG, down_fs);

%extract channel labels 
channel_labels = cell(EEG.nbchan, 1);
for i = 1:EEG.nbchan
    channel_labels(i) = cellstr(EEG.chanlocs(i).labels);  %%extract channel names from eeglab EEG data structure and create cell array 
end    

fprintf('\nAnalyzing file: %s\n', [input_dir '\' input_file '.bdf'])
fprintf('Loaded file contains %f min of data, sampled at %d Hz; %d channels (%d unused)\n', EEG.pnts/EEG.srate/60, EEG.srate, EEG.nbchan, length(unused_eeg_channels))

%extract triggers if available (these is are in data points, not seconds)
event_stop_times_pts = []; 
if ~isempty(EEG.event)
    event_start_idx = find([EEG.event.type]==767); %F2 press
    for i = 1:length(event_start_idx)
        event_stop_times_pts(i) = EEG.event(event_start_idx(i)).latency;
    end
end
event_stop_times = event_stop_times_pts./EEG.srate %stop times in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove mean signal from each channel (DC offset)  
disp('Removing mean from each channel')
for chan = 1:EEG.nbchan
        EEG.data(chan,:) = double(EEG.data(chan,:)) - mean(double(EEG.data(chan,:))); %remove mean signal from each channel (DC shift)   
end

%remove all freqs in notch_freq_vector
if (~isempty(notch_freq_vector) && notch_filter_on)
    fprintf('Performing notch filtering on EEG data +/- 2Hz ...\n');
    for i = notch_freq_vector
        eval(['[b,a] = ' eeg_filter_type '(eeg_filter_order, [i-2 i+2]./(EEG.srate/2),''stop'');'])
        for chan = eeg_channels
            EEG.data(chan,:) = filtfilt(b,a, double(EEG.data(chan,:)));
        end
    end
    fprintf('Performing notch filtering on REF data (if any) +/- 2Hz ...\n');
    for i = notch_freq_vector
        eval(['[b,a] = ' eeg_filter_type '(eeg_filter_order, [i-2 i+2]./(EEG.srate/2),''stop'');'])
        for chan = ref_channels
            EEG.data(chan,:) = filtfilt(b,a, double(EEG.data(chan,:)));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filter data 

beeg = [];
if (eeg_hicutoff > -1 && eeg_locutoff > -1)
    eval(['[beeg,a] = ' eeg_filter_type '(eeg_filter_order, [eeg_locutoff, eeg_hicutoff]./(EEG.srate/2),''bandpass'');'])
elseif (eeg_locutoff > -1)
    eval(['[beeg,a] = ' eeg_filter_type '(eeg_filter_order, eeg_locutoff/(EEG.srate/2),''high'');'])
elseif (eeg_hicutoff > -1)
    eval(['[beeg,a] = ' eeg_filter_type '(eeg_filter_order, eeg_hicutoff/(EEG.srate/2),''low'');'])    
end

if (length(beeg) > 1)
    disp('Filtering EEG data')
    for chan = eeg_channels
        EEG.data(chan,:) = filtfilt(beeg, a, double(EEG.data(chan,:)));   %filtfilt (zero phase filtering);
    end 
    for chan = unused_eeg_channels
        EEG.data(chan,:) = zeros(1,EEG.pnts);  %set unused channels to 0  
    end
end

baux= [];
if (aux_hicutoff > -1 && aux_locutoff > -1)
    eval(['[baux,a] = ' aux_filter_type '(aux_filter_order, [aux_locutoff, aux_hicutoff]./(EEG.srate/2),''bandpass'');'])
elseif (aux_locutoff > -1)
    eval(['[baux,a] = ' aux_filter_type '(aux_filter_order, aux_locutoff/(EEG.srate/2),''high'');'])
elseif (aux_hicutoff > -1)
    eval(['[baux,a] = ' aux_filter_type '(aux_filter_order, aux_hicutoff/(EEG.srate/2),''low'');'])
end
if (length(baux) > 1)
    disp('Filtering AUX data')
    for chan = aux_channels
        EEG.data(chan,:) = filtfilt(baux, a, double(EEG.data(chan,:)));   
    end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process data which will be saved into a file for later use (in matlab format)
eeg_data = zeros(length(eeg_channels), EEG.pnts);
aux_data = zeros(length(aux_channels), EEG.pnts);

for i = 1:length(eeg_channels)
    eeg_data(i,:) = EEG.data(eeg_channels(i),:);  
end 
eeg_channel_labels = channel_labels(eeg_channels);

for i = 1:length(aux_channels)
    aux_data(i,:) = EEG.data(aux_channels(i),:);  
end 

samp_freq = EEG.srate; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process EEG for artifact removal; 

if (~exist('analysis_channels', 'var'))
    analysis_channels = channel_labels(eeg_channels);
end

channels_ind = [];  
for i = 1:length(analysis_channels)
    flag = 0; ch_label = analysis_channels{i};
    for j = 1:length(channel_labels)
            if strcmp(channel_labels(j), ch_label)
                channels_ind = [channels_ind j]; flag = 1; break
            end
    end 
    if (flag == 0)
        fprintf('WARNING: Could not find channel %s in channel_labels.\n', ch_label)
    end    
end 
if (isempty(channels_ind))
    fprintf('ERROR: Could not find any channels for artifact analysis...exiting\n')
    return; 
end    

%set 'average' montage
avg_signal = mean(eeg_data(setdiff([1:size(eeg_data,1)],unused_eeg_channels),:)); 
for i = 1:size(eeg_data,1)
	eeg_data_avg(i,:) = eeg_data(i,:) - avg_signal; 
	if (ismember(i, unused_eeg_channels))  
    	eeg_data_avg(i,:) = zeros(1,length(eeg_data));   %set unused channel back to zero (above line will set it to avg_signal)
	end    
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot data 

if (disp_channel_labels == 1) %display EEG channel names
	eegplot(eeg_data_avg, 'srate', samp_freq, 'eloc_file', cap_locs, 'title', 'Processed EEG data for artifact removal - average montage') 
else %do NOT display EEG channel names (just numbers)
	eegplot(eeg_data_avg, 'srate', samp_freq, 'title', 'Processed EEG data for artifact removal - average montage') 
end

% eegplot(eeg_data, 'srate', EEG.srate, 'title', 'Original EEG data - unmontaged')  
% eegplot(EEG.data(65:end,:), 'srate', EEG.srate, 'title', 'All non-eeg channels')

my_T = (1/EEG.srate)*[0:length(eeg_data)-1];

% if (~isempty(aux_channels > 0))
% figure()
% for i = 1:length(aux_channels)
%     subplot(length(aux_channels), 1, i)
%     plot(my_T, aux_data(i,:))
%     title(['AUX ' num2str(i) ' '  char(aux_labels(i))])
% end 
% end


%calculate and plot spectrogram for one channel (average montage)
% spect_channel = 13; %C3 is 13
% WINDOW = 1024; %for power spectral density calculation, and Time-Freq analysis
% NOOVERLAP = 512; %for power spectral density calculation, and Time-Freq analysis
% NFFT = 1024; %for power spectral density calculation, and Time-Freq analysis (freq_res = Fs/NFFT)
% figure;
% [S,F,T2,P] = spectrogram(eeg_data_avg(spect_channel,:),WINDOW,NOOVERLAP,NFFT,EEG.srate);
% surf(T2,F,log10(P),'edgecolor','none');  
% axis tight; shading interp; view(0,90);
% set(gca, 'YLim', [0 130]);  %freq
% set(gca, 'CLim', [-1 1]);
% colorbar
% xlabel('Time (Seconds)'); ylabel('Hz');
% title(['Channel ',  char(analysis_channels(spect_channel))], 'FontSize', 7,'interpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare variables
rej_orig_samp_freq = []; rej_time = []; rej = [];
if (~isempty(notch_freq_vector) && notch_filter_on == 1)
    add_label = '_notch'; 
    for i = 1:length(notch_freq_vector)
        add_label = strcat(add_label, num2str(notch_freq_vector(i)));
    end
else
    add_label = ''; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%manually select artifacts (select artifact segments visually, then choose 'File' > 'Accept and Close', and Rejections will be saved in variable TMPREJ)
error('SELECT ARTIFACT SEGMENTS MANUALLY, THEN GO TO FIGURE > ACCEPT&CLOSE, THEN EXECUTE REMAINDER OF THE SCRIPT FROM LINE 245')
%%
rej1 = TMPREJ(:,1:2);  %TMPREJ has 78 columns (not sure why), first two columns are start and stop points for the manually-chosen segments (in data points)
rej=sortrows(rej1); %artifact sections will be not be ordered temporally, if they were not chosen like that, so need to sort in ascending order
%make sure that one artifact segment does not contain another one because that causes problems with movement analysis (this is rare, but happens, not sure why)
i=1;
while (i <= size(rej,1)-1)
    if (rej(i,2) >= rej(i+1,1))
        rej(i+1,:) = [];
        disp('REMOVING ARTIFACT SEGMENT')
    else
        i = i+1; 
    end    
end    
rej_time = rej./samp_freq; %convert to seconds 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    EXECUTE TO SAVE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %use this to append new rejections only
% save([output_dir '\' input_file '_f_' num2str(eeg_locutoff) '_' num2str(eeg_hicutoff), '_' eeg_filter_type 'ord' num2str(eeg_filter_order) '_' num2str(length(eeg_channels)) 'ch' add_label '.mat'], 'rej_time', '-append');

%save processed file
disp('Saving preprocessed file...')
if (size(eeg_data,2) < 2*10^6)
	save([output_dir '\' input_file '_f_' num2str(eeg_locutoff) '_' num2str(eeg_hicutoff), '_' eeg_filter_type 'ord' num2str(eeg_filter_order) '_' num2str(length(eeg_channels)) 'ch' add_label '.mat'], ...
    'eeg_data', 'aux_data', 'samp_freq', 'unused_eeg_channels', 'eeg_channel_labels', 'aux_labels', 'rej_time', 'event_stop_times_pts', 'event_stop_times' );

else
	disp('Large eeg matrix...using 7.3 switch')
    save([output_dir '\' input_file '_f_' num2str(eeg_locutoff) '_' num2str(eeg_hicutoff), '_' eeg_filter_type 'ord' num2str(eeg_filter_order) '_' num2str(length(eeg_channels)) 'ch' add_label '.mat'], ...
    'eeg_data', 'aux_data', 'samp_freq', 'unused_eeg_channels', 'eeg_channel_labels', 'aux_labels', 'rej_time', 'event_stop_times_pts', 'event_stop_times','-v7.3' );
end   
disp('Done...')