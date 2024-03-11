%% Predict heading and context using only CG trials 
% averaged per animal and two sample t test for prediction rates of FI and
% FS
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
allanimals = unique(baData.animalName);
FIheading = cell(length(allanimals),3);
FSheading = cell(length(allanimals),3);
FIcontext = cell(length(allanimals),3);
FScontext = cell(length(allanimals),3);

%% compute
for d = 1:3
    for a = 1:length(allanimals)
        celltype = 1; % 1 = FI, 0 = FS;
        allcells = baData.animalSessionCellName(ismember(baData.animalName, allanimals{a}) & baData.dayUsed == d & baData.isStable == celltype);
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

        FIheading{a,d} = nanmean(heading_prediction);
        FIcontext{a,d} = nanmean(context_prediction);

        % FS cells
        celltype = 0; % 1 = FI, 0 = FS;
         allcells = baData.animalSessionCellName(ismember(baData.animalName, allanimals{a}) & baData.dayUsed == d & baData.isStable == celltype);
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

        FSheading{a,d} = nanmean(heading_prediction);
        FScontext{a,d} = nanmean(context_prediction);
    end
end

%% FS vs FI

avg = nan(2,3);
sem = nan(2,3);
figure
for d = 1:3

    fprintf('day %d\n', d);
    fidata = cell2mat(FIheading(:,d));
    fsdata = cell2mat(FSheading(:,d));
    avg(1,d) = nanmean(fidata);
    avg(2,d) = nanmean(fsdata);
    sem(1,d) = std(fidata, 'omitnan') ./ sqrt(sum(~isnan(fidata)));
    sem(2,d) = std(fsdata, 'omitnan') ./ sqrt(sum(~isnan(fsdata)));
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
    

avg = nan(2,3);
sem = nan(2,3);
for d = 1:3
    fprintf('day %d\n', d);
    fidata = cell2mat(FIcontext(:,d));
    fsdata = cell2mat(FScontext(:,d));
    avg(1,d) = nanmean(fidata);
    avg(2,d) = nanmean(fsdata);
    sem(1,d) = std(fidata, 'omitnan') ./ sqrt(sum(~isnan(fidata)));
    sem(2,d) = std(fsdata, 'omitnan') ./ sqrt(sum(~isnan(fsdata)));
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