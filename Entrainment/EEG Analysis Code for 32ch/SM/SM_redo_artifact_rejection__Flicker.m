
%this is for visualizing and adjusting preprocessed files. Script will:
%1. plot eeg data with rejections marked
%2. user can redo artifact rejections and resave data

clear all

input_dir = 'C:\Users\mattokaren3\Box\Project_FlickerBloodFlow\EEG Testing\CABI\EEG';
input_file = 'MGtestFlickerAV_f_1_200_butterord3_32ch'; %without .mat extension

disp_channel_labels = 1; % 1 to display channel labels in eegplot; 0 to display channel number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([input_dir '\' input_file '.mat']);

fprintf('File: %s ... file is %.2f seconds long, sampling rate for analysis is %d Hz\n', [input_file '.mat'], length(eeg_data)/samp_freq, samp_freq)

%%% remove artifact data and concatenate the rest
if (exist('rej_time', 'var') && ~isempty(rej_time))
    if (~isempty(eeg_data))
        eeg_rejects = ones(1,size(eeg_data,2));
        for i = 1:size(rej_time,1) 
            eeg_rejects(1,floor(rej_time(i,1)*samp_freq)+1:floor(rej_time(i,2)*samp_freq)) = 0;  
        end    
        eeg_rejects = logical(eeg_rejects); 
        fprintf('%d points or %.2f seconds (%.2f percent) of data was manually cut from data\n', size(eeg_data,2) - sum(eeg_rejects), (size(eeg_data,2) - sum(eeg_rejects))/samp_freq, 100*(size(eeg_data,2) - sum(eeg_rejects))/size(eeg_data,2));  
        eeg_without_artif = eeg_data(:,eeg_rejects);
    end    
    if (~isempty(aux_data))
        aux_without_artif = aux_data(:,eeg_rejects);
    end   
else
    disp('WARNING: no artifact rejections in this file')
    eeg_without_artif = eeg_data; 
    aux_without_artif = aux_data; 
end

%set 'average' montage
avg_signal = mean(eeg_data(setdiff([1:size(eeg_data,1)],unused_eeg_channels),:)); 
for i = 1:size(eeg_data,1)
	eeg_data_avg(i,:) = eeg_data(i,:) - avg_signal; 
	if (ismember(i, unused_eeg_channels))  
    	eeg_data_avg(i,:) = zeros(1,length(eeg_data));   %set unused channel back to zero (above line will set it to avg_signal)
	end    
end  

avg_signal2 = mean(eeg_without_artif(setdiff([1:size(eeg_without_artif,1)],unused_eeg_channels),:)); 
for i = 1:size(eeg_without_artif,1)
	eeg_without_artif_avg(i,:) = eeg_without_artif(i,:) - avg_signal2; 
	if (ismember(i, unused_eeg_channels))  
    	eeg_without_artif_avg(i,:) = zeros(1,length(eeg_without_artif));   %set unused channel back to zero (above line will set it to avg_signal)
	end    
end  


%plot eeg data with artifact rejections marked
if (disp_channel_labels == 1)
	eegplot(eeg_data_avg,'srate', samp_freq, 'winrej', rej_time*samp_freq, 'eloc_file', 'C:\Users\mattokaren3\Box\Project_FlickerBloodFlow\EEG Testing\Standard-10-20-Cap64.locs', 'title', 'EEG - all data with artifact rejections marked');
else
	eegplot(eeg_data_avg,'srate', samp_freq, 'winrej', rej_time*samp_freq, 'title', 'EEG - all data with artifact rejections marked');
end
 
%plot spectrogram for motor cortex channels
ch_to_spect = [12 13]; %which channels to calculate spectrogram for
WINDOW = 2048;
NOOVERLAP = 1024;
NFFT = 2048;

figure;
for i = 1:length(ch_to_spect) 
    [S,F,T2,P] = spectrogram(eeg_data_avg(ch_to_spect(i),:),WINDOW,NOOVERLAP,NFFT,samp_freq);
    h = subplot(1,2,i);
    surf(T2,F,log10(P),'edgecolor','none');  
    axis tight; shading interp; view(0,90);
    %colormap; colorbar; 
    set(gca, 'YLim', [2 200]);  
    set(gca, 'CLim', [-1 1]); %WAS [0 3]  
    title(['ALL DATA: Channel ',  char(eeg_channel_labels(ch_to_spect(i))) ' Time-freq'], 'FontSize', 7,'interpreter','none')
    xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');
end
figure;
for i = 1:length(ch_to_spect) 
    [S,F,T2,P] = spectrogram(eeg_without_artif_avg(ch_to_spect(i),:),WINDOW,NOOVERLAP,NFFT,samp_freq);
    h = subplot(1,2,i);
    surf(T2,F,log10(P),'edgecolor','none');  
    axis tight; shading interp; view(0,90);
    %colormap; colorbar; 
    set(gca, 'YLim', [2 200]);  
    set(gca, 'CLim', [-1 1]); %WAS [0 3]  
    title(['WITHOUT ARTIFACT: Channel ',  char(eeg_channel_labels(ch_to_spect(i))) ' Time-freq'], 'FontSize', 7,'interpreter','none')
    xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');
end

% %%%if given time point in artifact_rejection time space, then calculate what that time point is in all_data time space
% tar = 176.5; %time in seconds (in artifact reject time space)
% artif_time_sum = 0; 
% for j = 1:size(rej_time,1) 
% 	tmp_Astart =  rej_time(j,1); tmp_Astop =  rej_time(j,2); 
% 	if (tmp_Astart <=  tar && tmp_Astop >=  tar) %tar is witin artifact period
%         artif_time_sum = artif_time_sum + (tar - tmp_Astart); 
%         break; 
%     else
%         artif_time_sum = artif_time_sum + (tmp_Astop - tmp_Astart); 
% 	end    
% end
% tor = tar + artif_time_sum


error('REVIEW ARTIFACT SELECTIONS AND IF MAKING CHANGES, EXECUTE FROM LINE 121') %forces a stop
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REDO ARTIFACT REJECTIONS
%%% SELECT ARTIFACT SEGMENTS MANUALLY, THEN GO TO FIGURE > ACCEPT&CLOSE, THEN EXECUTE REMAINDER OF THIS SECTION

rej1 = TMPREJ(:,1:2);  %TMPREJ has 78 columns (not sure why), first two columns are start and stop points for the manually-chosen segments (in data points)
rej=sortrows(rej1); %artifact sections will be not be ordered temporally, if they were not chosen like that, so need to sort in ascending order
%make sure that one artifact segment does not contain another one because that causes problems with movement analysis (this is rare, but happens, not sure why)
i=1;
while (i <= size(rej,1)-1)
    if (rej(i,2) >= rej(i+1,1))
        rej(i+1,:) = [];
        disp('NOTE: Removing artifact segment')
    else
        i = i+1; 
    end    
end    
rej_time = rej./samp_freq; %convert to seconds 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVE DATA

disp('Saving preprocessed file...')
if (size(eeg_data,2) < 2*10^6)
	save([input_dir '\' input_file] , 'rej_time', 'unused_eeg_channels', '-append');
else
	disp('Large eeg matrix...using 7.3 switch')
    save([input_dir '\' input_file] , 'rej_time', 'unused_eeg_channels', '-append', '-v7.3');          
end   






