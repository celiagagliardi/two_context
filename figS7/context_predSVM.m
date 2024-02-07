%% Context prediction using rates and SVMs
% Use Corr and Geo trials to predict context. Do prediction per cell type
% and see if any one cell type is driving the predictions, as per reviewer
% 2's suggestion. 
% CMG
% 1/31/24
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

fi_ba = baData(baData.isStable ==1,:);
fs_ba = baData(baData.isStable ==0,:);

fid = fopen(fullfile(output_path, 'context_pred_per_animal.txt'), 'w');
fprintf(fid, 'Context Prediction using rates, broken down by cell type, per animal\n\n');

validDigs = {'Corr', 'Geo'};
mincell = -inf;
mintrial = 4;
%% Perform Context prediction in FI cells. 
celltype = fi_ba;
cPred = cell(1,3);
for d = 1:3
    dtbl = celltype(celltype.dayUsed == d,:);
    sessions = unique([dtbl(:,1), dtbl(:,end)], 'rows');
    animals = unique(sessions.animalName);
    perAnimalPred = cell(length(animals),1);
    for a = 1:length(animals)
        cellnames = celltype.animalSessionCellName(ismember(celltype.animalName, animals{a}) & celltype.dayUsed == d);
        stbl = mapsData(ismember(mapsData.animalName, animals{a}) & mapsData.dayUsed == d,:);
        contextTbl = unique([stbl(:,6), stbl(:,7), stbl(:,9)], 'rows'); 
        tmprv = nan(length(contextTbl.trialId), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = mapsData(ismember(mapsData.animalSessionCellName, cellnames{c}),:);
            cellTrials = ctbl.trialId;
          if length(cellTrials) < mintrial
                fprintf('%s day %d cell %s does not have enough trials, skipping\n', animals{a}, d, cellnames{c});
                continue
            end
            tmprv(cellTrials, c) = ctbl.mfr;
        end
        % remove cells that do not participate.
        notallnans = ~all(isnan(tmprv),1);
        tmprv = tmprv(:, notallnans);
        tmprv(isnan(tmprv)) = 0;
        validTrials = contextTbl.trialId(ismember(contextTbl.contextId, [1,2]) ...
            & ismember(contextTbl.dig, validDigs));
        rv = tmprv(validTrials,:);
        contextClass = categorical(contextTbl.contextId(ismember(contextTbl.trialId, validTrials)));
       
        if size(rv, 2) < mincell
            fprintf('animal %s day %d does not have enough cells, skipping\n', animals{a}, d);
            continue
        end

        if sum(ismember(contextClass,'1')) < 2 || sum(ismember(contextClass,'2')) < 2
            fprintf('animal %s day %d does not have enough cellular data to predict context\n', animals{a}, d);
            continue
        end
        
        mdl = fitcsvm(rv, contextClass', 'KernelFunction', 'linear', 'Prior', 'empirical', 'Leaveout', 'on', 'Verbose',0); 
        err = kfoldLoss(mdl, 'mode', 'individual');
        predAcc = 1 - err;
        perAnimalPred{a} = nanmean(predAcc);
    end
    cPred{d} = perAnimalPred;
end

avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Context Prediction using FI cell rates\n')
fprintf(fid, 'FI cells\n');
for d = 1:3
    fprintf('Day %d\n', d)
    fprintf(fid, 'Day %d\n', d);
    cdata = cell2mat(cPred{d});
    avg(d) = mean(cdata, 'omitnan');
    sem1(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    sem2(d) = sqrt((mean(cdata, 'omitnan')*(1-mean(cdata, 'omitnan'))) ./ sum(~isnan(cdata)));

    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    %fprintf('\tContext Prediction on day %d p val = %.3f\n', d, p);

    fprintf('\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
    fprintf(fid, '\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
end

sem = sem1;

% Example data 
figure;
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


data1 = cell2mat(cPred{1});
data2 = cell2mat(cPred{2});
data3 = cell2mat(cPred{3});

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('FI cell context prediction');


%% Perform Heading prediction in FS cells. 
celltype = fs_ba;
cPred = cell(1,3);
for d = 1:3
    dtbl = celltype(celltype.dayUsed == d,:);
    sessions = unique([dtbl(:,1), dtbl(:,end)], 'rows');
    animals = unique(sessions.animalName);
    perAnimalPred = cell(length(animals),1);
    for a = 1:length(animals)
        cellnames = celltype.animalSessionCellName(ismember(celltype.animalName, animals{a}) & celltype.dayUsed == d);
        stbl = mapsData(ismember(mapsData.animalName, animals{a}) & mapsData.dayUsed == d,:);
        contextTbl = unique([stbl(:,6), stbl(:,7), stbl(:,9)], 'rows'); 
        tmprv = nan(length(contextTbl.trialId), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = mapsData(ismember(mapsData.animalSessionCellName, cellnames{c}),:);
            cellTrials = ctbl.trialId;
            if length(cellTrials) < mintrial
                fprintf('%s day %d cell %s does not have enough trials, skipping\n', animals{a}, d, cellnames{c});
                continue
            end
            tmprv(cellTrials, c) = ctbl.mfr;
        end
        % remove cells that do not participate.
        notallnans = ~all(isnan(tmprv),1);
        tmprv = tmprv(:, notallnans);
        tmprv(isnan(tmprv)) = 0;
        validTrials = contextTbl.trialId(ismember(contextTbl.contextId, [1,2]) ...
            & ismember(contextTbl.dig, validDigs));
        rv = tmprv(validTrials,:);
        contextClass = categorical(contextTbl.contextId(ismember(contextTbl.trialId, validTrials)));

        if size(rv, 2) < mincell
            fprintf('animal %s day %d does not have enough cells, skipping\n', animals{a}, d);
            continue
        end

        if sum(ismember(contextClass, '1')) < 2 || sum(ismember(contextClass, '2')) < 2
            fprintf('animal %s day %d does not have enough cellular data to predict context\n', animals{a}, d);
            continue
        end
        
        mdl = fitcsvm(rv, contextClass', 'KernelFunction', 'linear',  'Prior', 'empirical', 'Leaveout', 'on', 'Verbose',0); 
        err = kfoldLoss(mdl, 'mode', 'individual');
        predAcc = 1 - err;
        perAnimalPred{a} = nanmean(predAcc);
    end
    cPred{d} = perAnimalPred;
end

avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Context Prediction using FS cell rates\n')
fprintf(fid, 'FS Cells\n');
for d = 1:3
    fprintf('Day %d\n', d)
    fprintf(fid, 'Day %d\n', d);
    cdata = cell2mat(cPred{d});
    avg(d) = mean(cdata, 'omitnan');
    sem1(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    sem2(d) = sqrt((mean(cdata, 'omitnan')*(1-mean(cdata, 'omitnan'))) ./ sum(~isnan(cdata)));

    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    %fprintf('\tContext Prediction on day %d p val = %.3f\n', d, p);

    fprintf('\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
    fprintf(fid,'\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));

    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\tsign rank test: h = %d, p = %.3f\n', h, p)

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


data1 = cell2mat(cPred{1});
data2 = cell2mat(cPred{2});
data3 = cell2mat(cPred{3});

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('FS cell context prediction');

%%

sgtitle('Context Prediction per Animal');
fclose(fid);

saveas(gcf, fullfile(output_path, 'context_pred_per_animal.pdf'))


