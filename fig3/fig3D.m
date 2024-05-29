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

%% write to source data file

srcfile = '/Users/celia/Library/CloudStorage/OneDrive-UniversityofIowa/General - Two Context/celia_sourcedata/source_data.xlsx';
tmpT = best_aligned(:,[5, 22, 20, 21]);
tmpT = renamevars(tmpT, 'isStable', 'isFI');

writetable(tmpT, srcfile, 'Sheet', 'Fig3D');

tmpT = table(figData(:,1), figData(:,2), 'VariableNames', {'FI_percent', 'FS_percent'});
writetable(tmpT, srcfile, 'Sheet', 'Fig3D', 'Range', 'G1');


