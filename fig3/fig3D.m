%% fig3D

% load data
data_path = '/Users/celia/Documents/two_context_data'; % change this
cal = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tet = load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_results.mat'));
best_aligned = [cal.analysisResults.BestAligned; tet.analysisResults.BestAligned];
best_aligned = reformatTbl(best_aligned);

figData = nan(3,2); % each row is a day, each column is FI/FS
for d = 1:3
    tbl = best_aligned(best_aligned.dayUsed == d,:);
    total_cellNum = size(tbl,1);
    fi_cellNum = sum(tbl.isStable == 1);
    fs_cellNum = sum(tbl.isStable == 0);

    fprintf('Day %d\n', d);
    fprintf('\t %d total cells found, %d FI cells, %d FS cells\n', total_cellNum, fi_cellNum, fs_cellNum);
    fi_percent = (fi_cellNum / total_cellNum) * 100;
    fs_percent = (fs_cellNum / total_cellNum) * 100;
    figData(d,1) = fi_percent;
    figData(d,2) = fs_percent;
    fprintf('\t\t %.3f%% FI, %.3f%% FS\n', fi_percent, fs_percent);
end

figure;
barh(figData, 'stacked');

