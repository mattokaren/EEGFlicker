function [subjects,stims, fig] = plotEntrPerStim(figStruc, figData_dir)
%plotEntrPerStim Plots number of channels entrained per subject per stimuli
%condition
%   
    close all
    %%  Format structure input variable as a suitable table
    T = struct2table(figStruc);
    U = unstack(T,'num_ch_entrained','stimCondition');
    U.Properties.RowNames = U.subjectID;
    U.subjectID = [];
    
    %% Get subID and stimulus conditions for graph
    subjects = U.Properties.RowNames';
    stims = U.Properties.VariableNames;
    Markers = {'+','o','*','x','v','d','^','s','>','<'};
    
    %% Plot Entrainment for each subject
    fig = figure;
    for subject = 1:length(subjects)
        plot(U{subject,:},strcat('-',Markers{subject}))
        hold on
    end
    
    %% Graph Details
    set(gca, 'XTick', 1:length(stims), 'XTickLabel',stims);
    xlabel('Stimulus Condition');
    ylabel('Number of Channels Entained');
    title('Number of Channels Entrained for each Subject per Stimulus Condition');
    legend(subjects);
    
    %% Save Figure
    savefig(fig, [figData_dir 'fig3.fig'])

end

