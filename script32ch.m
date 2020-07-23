%% Script for 32-channel EEG processing and analysis
% addpath(genpath('C:\Users\matto\Box\Project_FlickerBloodFlow\EEG'))
%% Directory locations:
%folder that contains recording files, (.bdf's)
raw_folderName = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\rawdata\Experiment01\allfiles';
% raw_folderName = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\rawData';
%folder where cleaned data will go
cleaned_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\CleanedData\AutoClean';
% cleaned_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\CleanedData';

%folder where analyzed data/data structures will go
analyzed_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\AnalyzedResults\';
% analyzed_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\AnalyzedResults';

recordingList_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\EEGRecordingListV2.xlsx';

%% What script will do:
cleanThenAnalyze = 0;
reanalyze = 1;
checkSpreadSheet = 0;
updateSpreadSheet = 0;
manualClean = 0; %1 if you would like to manually
% analyzedStimuliFilesOnly = 0;

% %% Run EEG Lab
% if ~exist("EEG", "var")
%     eeglab
% end

%%
if checkSpreadSheet
     T = table2array(readtable(recordingList_dir));  %open EEGrecordingSpreadsheet
     %check for 1s in runCleaning col
     %check for 1s in runAnalysis ocol
end
%%  Check files in excel sheet
% if analyzedStimuliFilesOnly
% end
%entrainedNum_per_recording = int16.empty(5,0);
entrainedNum_per_recording = [];
%% Iterate thru all bdf files in folder to clean, then analyze cleaned files
if cleanThenAnalyze
    files = dir(fullfile(raw_folderName,'*.bdf'));
    
    % entrainedNum_per_recording = int16.empty(max(size(files)),0);
    [~,idx] = sort([files.datenum]);
    for k = 1:numel(files)
%         close all  
            fileName = [files(idx(k)).folder, '\', files(idx(k)).name];
            cleaned_file = preprocess32ch(fileName,cleaned_dir);
    %         saveCleanedData()
            [analyzed_filename, num_ch_entrained] = Analysis32ch(cleaned_file, analyzed_dir);
            entrainedNum_per_recording = [entrainedNum_per_recording num_ch_entrained];
    end
    %disp(['entrainedNum_per_recording: ', entrainedNum_per_recording]);
end

%% Reanalyze cleaned data.  For data that has already been cleaned
if reanalyze
    % get all files in cleaned folder directory
    files = dir(fullfile(cleaned_dir,'*.mat')); %#ok<*UNRCH>
    for k = 1:numel(files)  %iterate thru all cleaned files
        close all
            cleaned_file = [files(k).folder, '\', files(k).name];
%             cleaned_file = files(k).name;
    %         saveCleanedData()
            [analyzed_filename, num_ch_entrained] = Analysis32ch072220(cleaned_file, analyzed_dir);
    end
end

%% End Message 
disp(['Script complete. Ran for ' , num2str(k), ' files.']);