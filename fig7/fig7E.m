%% contextPrediction_SVM
% The one that's used in the figure. Context Prediction using Geo axis
% trials (digs are Correct or Geo).
data_path = '/Users/celia/Documents/two_context_data'; % change this
load(fullfile(data_path, 'ratemaps-2021_12_1_12_34_49_K1_d4.mat'));
minTrial = 4;


fid = fopen('Fig7Edata.csv', 'w');
fprintf(fid, 'animal, day, acc\n');

cPred = cell(7,3);
fprintf('Context Prediction using only Corr/Geo Trials\n');
for a = 1:length(ratemaps.data)
    for s = 1:3
        if ismember(ratemaps.data(a).animal, 'HGY1_recut')
            if s == 2
                fprintf('HGY1_recut does not have day 2 data\n');
                continue
            elseif s == 3
                fprintf('HGY1_recut does not have day 3 data\n')
                continue
            end
        end

        allcells = findUniqueCells(ratemaps.data(a).session(s));
        contexts = extractfield(ratemaps.data(a).session(s).trial, 'context');
        digs = extractfield(ratemaps.data(a).session(s).trial, 'dig');
        trialNum = extractfield(ratemaps.data(a).session(s).trial, 'trialNum');
        trialsAvail = trialNum(ismember(digs, {'Corr', 'Geo'}));
        if ~any(ismember(digs, 'Corr')) || ~any(ismember(digs, 'Geo'))
            fprintf('%s does not have any Corr or Geo digs this session: %s, skipping...\n', ratemaps.data(a).animal, ratemaps.data(a).session(sIdx(s)).name);
            continue;
        end
        rv = nan(length(trialsAvail), length(allcells));
        for c = 1:length(allcells)
            mfr = calculateMfrVector(ratemaps.data(a).session(s), allcells{c});
            cellTrialsAvail = trialNum(~isnan(mfr));
            cellCounts = length(cellTrialsAvail) >= minTrial & sum(contexts(cellTrialsAvail) == 1) > 2 & sum(contexts(cellTrialsAvail) == 2) > 2;

            if ~cellCounts
                fprintf('%s, day %d, %s is disqualified from analysis\n', ratemaps.data(a).animal, s, allcells{c});
                continue
            end
            rv(:,c) = mfr(trialsAvail);
        end
        allcells = allcells(any(rv,1));
        rv = rv(:, any(rv,1));
        rv(isnan(rv)) = 0; % this step might be sketchy but you can't train a svm if theres a nan value anywhere.

        contextClass = contexts(trialsAvail)';
        contextClass = categorical(contextClass);
        mdl = fitcsvm(rv, contextClass, 'KernelFunction', 'linear', 'PredictorNames', allcells, 'Prior', 'empirical', 'Leaveout', 'on', 'Verbose',0);


        err = kfoldLoss(mdl, 'mode', 'individual');
        predAcc = 1 - err;
        cPred{a, s} = nanmean(predAcc);
        fprintf(fid, '%s, %d, %.3f\n', animals{a}, s, nanmean(predAcc));

    end
end

avg = nan(3,1);
sem = nan(3,1);

for d = 1:3
    fprintf('Day %d\n', d)
    cdata = cell2mat(cPred(:,d));
    avg(d) = mean(cdata);
    sem(d) = std(cdata) ./ sqrt(length(cdata));

    fprintf('\t avg = %.4f, sem = %.4f using %d animals\n',  avg(d), sem(d), length(cdata));

    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\t\t sign rank test: h = %d, p = %.3f using %d animals\n', h, p, length(cell2mat(cPred(:,d))));
    [h,p,ci,stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t\t one sample ttest: t(%d) = %.3f, p = %.5f\n', stats.df, stats.tstat, p);
  

end


% Example data
figure;
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
ylim([0 0.9])
title('Context Prediction using Corr/Geo Trials Only');

data1 = cell2mat(cPred(:,1));
data2 = cell2mat(cPred(:,2));
data3 = cell2mat(cPred(:,3));

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

fclose(fid);

%% 
data = readtable('Fig7Edata.csv');
fprintf('Context Prediction');
avg = nan(3,1);
sem = nan(3,1);
for d = 1:3
    fprintf('\tDay %d\n', d);
    cdata = data.acc(data.day == d);
    davg = nanmean(cdata);
    dsem = std(cdata) / sqrt(length(cdata));
    avg(d) = davg; sem(d) = dsem;
    fprintf('\t\tavg = %.3f, sem = %.3f using %d animals\n', davg, dsem, length(cdata));
    
    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\t\t\t sign rank test: h = %d, p = %.3f\n', h, p);
    [h,p,ci,stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('\t\t\t one sample ttest: t(%d) = %.3f, p = %.5f\n', stats.df, stats.tstat, p);
end

figure;
b = bar(avg, 'grouped');
hold on
[ngroups, nbars] = size(avg);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
ylim([0 1])
title('Context Prediction using Corr/Geo Trials Only');

scatter(1, data.acc(data.day == 1), 'ro', 'xjitter', 'rand');
scatter(2, data.acc(data.day == 2), 'bo', 'xjitter', 'rand');
scatter(3, data.acc(data.day == 3), 'co', 'xjitter', 'rand');

  
