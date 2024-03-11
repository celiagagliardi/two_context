%% Predict heading and context using only CG trials 
% Rate as the predictor
% Predict one cell at a time across days.

clear all;
data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_input.mat'));
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));

output_path = '/Users/celia/Documents/two_context_data/figures/SVMpred_celltype';
if ~isfolder(output_path)
    mkdir(output_path);
end

mapsData = reformatTbl(analysisInput.MapsData);
baData = reformatTbl(analysisResults.BestAligned);

northDigs = {'Corr'};
southDigs = {'Geo'};
% add any parameters here

FIheading = cell(1,3);
FSheading = cell(1,3);
FIcontext = cell(1,3);
FScontext = cell(1,3);

%% compute
for d = 1:3
    celltype = 1; % 1 = FI, 0 = FS;
    allcells = baData.animalSessionCellName(baData.dayUsed == d & baData.isStable == celltype,:);
    heading_prediction = nan(length(allcells),1);
    context_prediction = nan(length(allcells),1);
    for c = 1:length(allcells)
        mtbl = mapsData(ismember(mapsData.animalSessionCellName, allcells{c}),:);
        validTrials = mtbl(ismember(mtbl.dig, [northDigs, southDigs]),:);
        mfr = validTrials.mfr;
        mfr(isnan(mfr)) = 0;
        digClass = categorical(validTrials.dig);
        contextClass = categorical(validTrials.contextId);

        digMdl = fitcsvm(mfr, digClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'leaveout', 'on');
        err = kfoldLoss(digMdl);
        acc = 1 - err;
        heading_prediction(c) = acc;

        conMdl = fitcsvm(mfr, contextClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'leaveout', 'on');
        err = kfoldLoss(conMdl);
        acc = 1 - err;
        context_prediction(c) = acc;
    end

    FIheading{d} = heading_prediction;
    FIcontext{d} = context_prediction;
    
    % FS cells
    celltype = 0; % 1 = FI, 0 = FS;
    allcells = baData.animalSessionCellName(baData.dayUsed == d & baData.isStable == celltype,:);
    heading_prediction = nan(length(allcells),1);
    context_prediction = nan(length(allcells),1);
    for c = 1:length(allcells)
        mtbl = mapsData(ismember(mapsData.animalSessionCellName, allcells{c}),:);
        validTrials = mtbl(ismember(mtbl.dig, [northDigs, southDigs]),:);
        mfr = validTrials.mfr;
        mfr(isnan(mfr)) = 0;
        digClass = categorical(validTrials.dig);
        contextClass = categorical(validTrials.contextId);

        digMdl = fitcsvm(mfr, digClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'leaveout', 'on');
        err = kfoldLoss(digMdl);
        acc = 1 - err;
        heading_prediction(c) = acc;

        conMdl = fitcsvm(mfr, contextClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'leaveout', 'on');
        err = kfoldLoss(conMdl);
        acc = 1 - err;
        context_prediction(c) = acc;
    end

    FSheading{d} = heading_prediction;
    FScontext{d} = context_prediction;
end

%% Stats
digavgFI = nan(1,3);
digsemFI = nan(1,3);
conavgFI = nan(1,3);
consemFI = nan(1,3);
digavgFS = nan(1,3);
digsemFS = nan(1,3);
conavgFS = nan(1,3);
consemFS = nan(1,3);

for d = 1:3
    fprintf('Day %d FI cells: \n', d);
    cdata = FIheading{d};
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t Heading Prediction avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d cells\n', ...
        avg, sem, h, p, sum(~isnan(cdata)));
    digavgFI(d) = avg; digsemFI(d) = sem;
   
    cdata = FIcontext{d};
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t Context Prediction avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d cells\n', ...
        avg, sem, h, p, sum(~isnan(cdata)));
    conavgFI(d) = avg; consemFI(d) = sem;

    fprintf('Day %d FS cells: \n',d);
    cdata = FSheading{d};
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t Heading Prediction avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d cells\n', ...
        avg, sem, h, p, sum(~isnan(cdata)));
    digavgFS(d) = avg; digsemFS(d) = sem;
   
    cdata = FScontext{d};
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t Context Prediction avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d cells\n', ...
        avg, sem, h, p, sum(~isnan(cdata)));
    conavgFS(d) = avg; consemFS(d) = sem;
end

%% plot
% Example data 
figure;
subplot(2,2,1);
bar(digavgFI); 
hold on;
errorbar(digavgFI, digsemFI, 'k--', 'linestyle', 'none');
yline(.5)
ylim([0 1]);
hold off;
title('FI Heading Prediction');

subplot(2,2,2)
bar(conavgFI); 
hold on;
errorbar(conavgFI, consemFI, 'k--',  'linestyle', 'none');
yline(.5)
ylim([0 1]);
hold off;
title('FI Context Prediction');

subplot(2,2,3);
bar(digavgFS); 
hold on;
errorbar(digavgFS, digsemFS, 'k--', 'linestyle', 'none');
yline(.5)
ylim([0 1]);
hold off;
title('FS Heading Prediction');

subplot(2,2,4)
bar(conavgFS); 
hold on;
errorbar(conavgFS, consemFS, 'k--', 'linestyle', 'none');
yline(.5)
ylim([0 1]);
hold off;
title('FS Context Prediction');



% data1 = cell2mat(cPredH(:,1));
% data2 = cell2mat(cPredH(:,2));
% data3 = cell2mat(cPredH(:,3));
% 
% scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
% scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
% scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')
% 
% hold off
% ylim([0 1])
% ylabel('Prediction Accuracy');
% xlabel('Day');
% title(sprintf('%s cell heading prediction', cell2analyze));
% 
%% FS vs FI

avg = nan(3,2);
sem = nan(3,2);
figure
for d = 1:3

    fprintf('day %d\n', d);
    fidata = cell2mat(FIheading(:,d));
    fsdata = cell2mat(FSheading(:,d));
    avg(1,d) = nanmean(fidata);
    avg(2,d) = nanmean(fsdata);
    sem(1,d) = std(fidata) ./ sqrt(sum(~isnan(fidata)));
    sem(2,d) = std(fsdata) ./ sqrt(sum(~isnan(fsdata)));
    fprintf('\tFI Heading prediction avg = %.3f, sem = %.3f\n', avg(1,d), sem(1,d));
    fprintf('\tFS Heading prediction avg = %.3f, sem = %.3f\n', avg(2,d), sem(2,d));
    [h,p, stats, ci] = ttest2(fidata, fsdata);
    fprintf('\t\t FI vs FS data: h = %d, p = %.5f\n', h, p);
end
avg = avg';
sem = sem';
subplot(1,2,1)
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
yline(.5)
hold off
ylim([0 1])
title('Heading Prediction Accuracy');
legend('FI', 'FS')
xlabel('Day');
ylabel('Prediction Accuracy');
    

avg = nan(3,2);
sem = nan(3,2);
for d = 1:3
    fprintf('day %d\n', d);
    fidata = cell2mat(FIcontext(:,d));
    fsdata = cell2mat(FScontext(:,d));
    avg(1,d) = nanmean(fidata);
    avg(2,d) = nanmean(fsdata);
    sem(1,d) = std(fidata) ./ sqrt(sum(~isnan(fidata)));
    sem(2,d) = std(fsdata) ./ sqrt(sum(~isnan(fsdata)));
    fprintf('\tFI Context prediction avg = %.3f, sem = %.3f\n', avg(1,d), sem(1,d));
    fprintf('\tFS Context prediction avg = %.3f, sem = %.3f\n', avg(2,d), sem(2,d));
    [h,p, stats, ci] = ttest2(fidata, fsdata);
    fprintf('\t\t FI vs FS data: h = %d, p = %.5f\n', h, p);
end
avg = avg';
sem = sem';

subplot(1,2,2)
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
yline(.5)
hold off
ylim([0 1])
ylabel('Prediction Accuracy');
legend('FI', 'FS')
xlabel('Day');
title('Context Prediction');