function [] = entrainCount()
%entrainCount Collects the number of channels entrained from each recording of
%interest from RESULTS.MAT files
%   Detailed explanation goes here

%% Directory locations:
    analyzed_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\AnalyzedResults\BipolarMontage\AutoClean';
    figData_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\Fig3Data\';

    recordingList_dir = 'C:\Users\matto\Box\Project_FlickerBloodFlow\EEG\CABI Recordings\EEGRecordingListV2.xlsx';
    %  recordingList_dir = 'C:\Users\mattokaren3\Box\Project_FlickerBloodFlow\EEG\AnalysisCode\EEGRecordingListV2.xlsx';

    checkTestCol = 20;  % col in excel sheet that has 1 if to run figure analysis on that recording (column 20)
 
 %% Read Recording List Excel Sheet and get recording files of interest
    T = readtable(recordingList_dir); % imports excel as a Table in matlab
    Arr = table2array(T); % converts table to arr

    % get rownumbers of recordings with a 1 in the RunFigure3 col
    [figFiles , ~] = find(ismember(Arr(:,checkTestCol), '1'));
    %  isOne = cellfun(@(x)isequal(x,1),T);
    %  [row,col] = find(isOne);
    %  [~, col] = find(ismember(T(:,20), '1'));  %check column header
    %  [r c] = find(strcmp([T{:}], '1'))
    %  col = table2array(T(:,1))
    figStruc = struct;  %create structure called 'figStruct'
    resultsFiles = dir(fullfile(analyzed_dir,'*.mat'));  % get all RESULTs.mat files in analyzed folder

 %% Iterate thru all Results Files
    for iResult = 1:numel(figFiles) % 
    %   close all  

%         check if to be checked
        recordingFileName = char(Arr(figFiles(iResult),1));  %get EEG recording filename
    %         index = contains(files.name,fileName);
%         index = find(strcmp({resultsFiles.name}, fileName)==1);
    %         resultsFilesNames = resultsFiles.name;
        resultsFilesNames = extractfield(resultsFiles,'name')'; %get list of analyzed EEG results

        %find file that matches fileName
%         fileNames = resultsFiles.name;
%         fileNames = resultsFiles.name;
    %         idx = find(ismember(resultsFiles.name, fileName))
        [resultsFileROW, ~] = find(contains([resultsFilesNames],recordingFileName));  % find results file that corresponds to EEG recording of interest
        if isempty(resultsFileROW)
            disp(['Results for ', recordingFileName, ' does not exist in this directory.'])
         
        else
        %         [a b] = find(contains(T(:,1),fileName)); % search for file name in excel
        %         if str2double(T(a,checkTestCol))

            fullFileName = [resultsFiles(resultsFileROW).folder, '\', resultsFiles(resultsFileROW).name];

            % get info to be saved into new struct
            load(fullFileName,'num_ch_entrained'); % load variable for # ch entrained for this recording
            subjectID = char(Arr(figFiles(iResult),5));
            stimCondition = char(Arr(figFiles(iResult),8));

            % save info into the new structure
            figStruc(iResult).subjectID = subjectID;
            figStruc(iResult).stimCondition = stimCondition;
            figStruc(iResult).num_ch_entrained = num_ch_entrained;
        %         figStruc(iResul = [subjectID,stimCondition, num_ch_entrained];
        end
    end
    
%% Save and Plot Results

% save(figStruct)
[subjects,stims, fig] = plotEntrPerStim(figStruc, figData_dir);
end

