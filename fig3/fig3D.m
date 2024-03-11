%% fig3D

data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));
tmptbl = reformatTbl(analysisResults.BestAlignedCounts);
load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tmptbl1 = reformatTbl(analysisResults.BestAlignedCounts);

tbl = [tmptbl; tmptbl1];
y = nan(2,3);
for d = 1:3
    sumStable = sum(tbl.numStable(tbl.dayUsed == d));
    sumUnstable = sum(tbl.numUnstable(tbl.dayUsed == d));
    y(:,d) = [100 * sumStable /(sumStable + sumUnstable); 100*sumUnstable/(sumUnstable+sumStable)];
    fprintf('Day %d: %d (%.3f) FI and %d (%.3f) FS cells\n', d, sumStable, 100*sumStable/(sumStable+sumUnstable), sumUnstable, 100*sumUnstable/(sumUnstable+sumStable));
end

figure;
barh(flipud(y'), 'stacked');

%% Supplemental figure
data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));
tbl = reformatTbl(analysisResults.BestAlignedCounts);

avg = nan(1,3);
sem = nan(1,3);
for d = 1:3
    dtbl = tbl(tbl.dayUsed == d,:);
    avg(d) = mean(dtbl.percentStable);
    sem(d) = std(dtbl.percentStable) ./ sqrt(size(dtbl,1));
end

figure;
bar(avg);
hold on;
errorbar(avg, sem, 'k', 'linestyle', 'none');
scatter(repmat(1, sum(tbl.dayUsed ==1),1 ), tbl.percentStable(tbl.dayUsed == 1), 'jitter', 'on')
scatter(repmat(2, sum(tbl.dayUsed ==2),1 ), tbl.percentStable(tbl.dayUsed == 2), 'jitter', 'on')
scatter(repmat(3, sum(tbl.dayUsed ==3),1 ), tbl.percentStable(tbl.dayUsed == 3), 'jitter', 'on')

title('ephys');

load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tmptbl1 = reformatTbl(analysisResults.BestAlignedCounts);

tbl = reformatTbl(analysisResults.BestAlignedCounts);

avg = nan(1,3);
sem = nan(1,3);
for d = 1:3
    dtbl = tbl(tbl.dayUsed == d,:);
    avg(d) = mean(dtbl.percentStable);
    sem(d) = std(dtbl.percentStable) ./ sqrt(size(dtbl,1));
end

figure;
bar(avg);
hold on;
errorbar(avg, sem, 'k', 'linestyle', 'none');
scatter(repmat(1, sum(tbl.dayUsed ==1),1 ), tbl.percentStable(tbl.dayUsed == 1), 'jitter', 'on')
scatter(repmat(2, sum(tbl.dayUsed ==2),1 ), tbl.percentStable(tbl.dayUsed == 2), 'jitter', 'on')
scatter(repmat(3, sum(tbl.dayUsed ==3),1 ), tbl.percentStable(tbl.dayUsed == 3), 'jitter', 'on')

title('calcium');