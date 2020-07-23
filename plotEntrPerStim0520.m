function [stims,y] = plotEntrPerStim()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all
stims = {'Baseline','7k','8k','10k','LED Only'};
% y = [0 10 5];
% z = [0 10 5; 2 3 4]';
S3 = [0 0 1 1 0];
S4 = [0 0 1 0 11];
S5 = [0 7 4 9 0];
S6 = [0 0 0 0 5];
subjects = {'S3','S4','S5','S6'};
%     hold on;
       plot(S3)
       hold on;
       plot(S4)
       plot(S5)
       plot(S6)
%        bar([neworder{:,2}])
       set(gca, 'XTick', 1:5, 'XTickLabel',stims);
       xlabel('Stimulus Condition');
       ylabel('Number of Channels Entained');
       title('Number of Channels Entrained for each Subject per Stimulus Condition');
       legend(subjects);
       %        hold off;
% plot(x,y)
end

