%% Heading Prediction SVM
% Figure 7 E. Using firing rates to predict context in all cells.
%% load data
data_path = '/Users/celia/Documents/two_context_data'; % change this
tetin = load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_input.mat'));
maps = reformatTbl(tetin.analysisInput.MapsData);

minTrial = 2; 
validDigs = {'Corr', 'Geo'};
animals = unique(maps.animalName);
cPred = cell(length(animals),3);
fprintf('Heading prediction\n');
for a = 1:length(animals)
    for d = 1:3
        tbl = maps(ismember(maps.animalName, animals{a}) & maps.dayUsed == d & ismember(maps.dig, validDigs),:);
        if isempty(tbl)
            continue; % one animals only has data for day 1
        end
        trialId = unique(tbl.trialId);
        digId = cell(length(trialId),1);
        for t = 1:length(trialId)
            trialDig = unique(tbl.dig(tbl.trialId == trialId(t)));
            digId(t) = trialDig;
        end
        digId = categorical(digId);
        cellnames = unique(tbl.animalSessionCellName);
        rv = zeros(length(trialId), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = tbl(ismember(tbl.animalSessionCellName, cellnames{c}),:);
            cellTrials = ctbl.trialId;
            numCorrTrials = sum(ismember(ctbl.dig, 'Corr'));
            numGeoTrials = sum(ismember(ctbl.dig, 'Geo'));
            if length(cellTrials) >= minTrial && numCorrTrials >= 1 && numGeoTrials >= 1
                for t = 1:length(cellTrials)
                    trialIdx = find(cellTrials(t) == trialId);
                    rv(trialIdx, c) = ctbl.mfr(ctbl.trialId == cellTrials(t));
                end
            else
                fprintf('%s day %d cell %d excluded from analysis\n', animals{a}, d, c);
            end
        end

        cellsIncluded = any(rv);
        rv = rv(:, cellsIncluded);
      
        mdl = fitcsvm(rv, digId, 'kernelfunction', 'linear', 'leaveout', 'on', 'verbose', 0);
        err = kfoldLoss(mdl);
        acc = 1 - err;
        cPred{a,d} = acc;
    end
end

% plot

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
ylim([0 1])
title('Heading Prediction using Corr/Geo Trials Only');

data1 = cell2mat(cPred(:,1));
data2 = cell2mat(cPred(:,2));
data3 = cell2mat(cPred(:,3));

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

