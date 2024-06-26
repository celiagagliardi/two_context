%% Figure 3C
% load data
data_path = '/Users/celia/Documents/two_context_data'; % change this
cal = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tet = load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_results.mat'));

best_aligned = [cal.analysisResults.BestAligned; tet.analysisResults.BestAligned];
best_aligned = reformatTbl(best_aligned); % just a helper function to make indexing session days easier
fprintf('total cell number = %d, cells with place fields in both contexts = %d\n', ...
    size(best_aligned,1), sum(~isnan(best_aligned.bestCorrelation)));

figure;
binWidth = .02;
histogram(best_aligned.bestCorrelation, 'BinWidth', binWidth, 'normalization', 'probability');

[N, edges] = histcounts(best_aligned.bestCorrelation, 'binWidth', binWidth);
[probN, edges] = histcounts(best_aligned.bestCorrelation, 'binWidth', binWidth, 'normalization', 'probability');

binCenters = edges(1:end-1) + diff(edges) / 2;
% 
% writetable(best_aligned(:,[5,20]), 'source_data.xlsx', 'Sheet', 'Fig3D');
% 
% tmpT = table(binCenters', N', probN', 'VariableNames', {'BinCenters', 'Hist_counts', 'Hist_probability'});
% writetable(tmpT, 'source_data.xlsx', 'Sheet', 'Fig3D', 'Range', 'E1');
% 
