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
mincell = -inf;
mintrial = 4;
figure;
%% ALL cells
celltype = baData(baData.isStable == 1,:);
allCpredFi= cell(1,3); %context pred
allHpredFi = cell(1,3); % heading pred
cPredH = cell(7,3);
cPredC = cell(7,3);

truedigs = cell(7,3);
truecon = cell(7,3);
predDigs = cell(7,3);
predCons = cell(7,3);

fprintf('all cells\n\n');
for d = 1:3
    dtbl = celltype(celltype.dayUsed == d,:);
    sessions = unique([dtbl(:,1), dtbl(:,end)], 'rows');
    animals = unique(sessions.animalName);
    contextPred = cell(length(animals),1);
    headingPred = cell(length(animals),1);

    TH = [];
    TC = [];
    for a = 1:length(animals)
        cellnames = celltype.animalSessionCellName(ismember(celltype.animalName, animals{a}) & celltype.dayUsed == d);
        stbl = mapsData(ismember(mapsData.animalName, animals{a}) & mapsData.dayUsed == d,:);
        trialTbl = unique([stbl(:,6), stbl(:,7), stbl(:,9)], 'rows'); 
        tmprv = nan(max(trialTbl.trialId), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = mapsData(ismember(mapsData.animalSessionCellName, cellnames{c}),:);
            cellTrials = ctbl.trialId;
            tmprv(cellTrials, c) = ctbl.mfr;
        end
        notallnans = ~all(isnan(tmprv),1);
        tmprv = tmprv(:, notallnans);
        tmprv(isnan(tmprv)) = 0;
        validTrials = trialTbl.trialId(ismember(trialTbl.contextId, [1,2]) ...
            & ismember(trialTbl.dig, [northDigs, southDigs]));
        rv = tmprv(validTrials,:);
        if size(rv, 2) < mincell
            fprintf('animal %s day %d does not have enough cells, skipping\n', animals{a}, d);
            continue
        end

        contextClass = categorical(trialTbl.contextId(ismember(trialTbl.trialId, validTrials)));
        digClass = categorical(trialTbl.dig(ismember(trialTbl.trialId, validTrials)));

        hMdl = fitcsvm(rv, digClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'leaveout', 'on', 'verbose', 0);
        err = kfoldLoss(hMdl, 'mode', 'individual');
        predAcc = 1 - err;
        cPredH{a,d} = nanmean(predAcc);
        cMdl = fitcsvm(rv, contextClass, 'kernelfunction', 'linear', 'prior', 'empirical', 'leaveout', 'on', 'verbose', 0);
        err = kfoldLoss(cMdl, 'mode', 'individual');
        predAcc = 1 - err;
        cPredC{a,d} = nanmean(predAcc);

        iT = calc_svm_scores(hMdl, rv, digClass);
        TH = [TH; iT];
        iT = calc_svm_scores(cMdl, rv, contextClass);
        TC = [TC; iT];


    end
   
end

%%  plot the results
avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Context Prediction using ALL cell rates\n')
%fprintf(fid, 'FI cells\n');
for d = 1:3
    fprintf('Day %d\n', d)
    cdata = cell2mat(cPredC(:,d));
    avg(d) = mean(cdata, 'omitnan');
    sem1(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    sem2(d) = sqrt((mean(cdata, 'omitnan')*(1-mean(cdata, 'omitnan'))) ./ sum(~isnan(cdata)));
% 
%     [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
%     %fprintf('\tContext Prediction on day %d p val = %.3f\n', d, p);
% 
%     fprintf('\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
% 
%     %fprintf(fid, '\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));

    [p,h] = signrank(cdata, .5, 'tail', 'right');

    fprintf('\t avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg(d), sem1(d), h, p, sum(~isnan(cdata)));
   
end

sem = sem1;

% Example data 
subplot(1,2,1);
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


data1 = cell2mat(cPredC(:,1));
data2 = cell2mat(cPredC(:,2));
data3 = cell2mat(cPredC(:,3));

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('ALL cell context prediction');

% heading prediction      

avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Heading Prediction using ALL cell rates\n')
for d = 1:3
    fprintf('Day %d\n', d)
    cdata = cell2mat(cPredH(:,d));
    avg(d) = mean(cdata, 'omitnan');
    sem1(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    sem2(d) = sqrt((mean(cdata, 'omitnan')*(1-mean(cdata, 'omitnan'))) ./ sum(~isnan(cdata)));

    %[h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    %fprintf('\tContext Prediction on day %d p val = %.3f\n', d, p);
    [p,h] = signrank(cdata, .5, 'tail', 'right');

    fprintf('\t avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg(d), sem1(d), h, p, sum(~isnan(cdata)));
    %fprintf(fid, '\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
end

sem = sem1;

% Example data 
subplot(1,2,2);
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


data1 = cell2mat(cPredH(:,1));
data2 = cell2mat(cPredH(:,2));
data3 = cell2mat(cPredH(:,3));

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('ALL cell heading prediction');

sgtitle('with empirical prior');
       