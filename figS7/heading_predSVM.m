%% Heading prediction using rates and SVMs
% CMG
% 1/25/24

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

northDigs = {'Corr'};
southDigs = {'Geo'};

fid = fopen(fullfile(output_path, 'heading_pred_per_animal.txt'), 'w');
fprintf(fid, 'Heading Prediction using rates, broken down by cell type, per animal\n\n');

mincell = -inf;
mintrial = 4;
%% Perform Heading prediction in FI cells. 
celltype = fi_ba;
%sessions = unique([celltype(:,1), celltype(:,end)], 'rows');
cPred = cell(1,3);
for d = 1:3
    dtbl = celltype(celltype.dayUsed == d,:);
    sessions = unique([dtbl(:,1), dtbl(:,end)], 'rows');
    animals = unique(sessions.animalName);
    perAnimalPred = cell(length(animals),1);
    for a = 1:length(animals)
        cellnames = celltype.animalSessionCellName(ismember(celltype.animalName, animals{a}) & celltype.dayUsed == d);
        stbl = mapsData(ismember(mapsData.animalName, animals{a}) & mapsData.dayUsed == d,:);
        digs = unique([stbl(:,6), stbl(:,9)], 'rows');
        tmprv = nan(length(digs.trialId), length(cellnames));
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
        validDigTrials = digs.trialId(ismember(digs.dig, [northDigs, southDigs]));
        digClass = categorical(digs.dig(ismember(digs.trialId, validDigTrials)));

        rv = tmprv(validDigTrials,:);
        
        if size(rv,2) < mincell
            fprintf('Animal %s day %d doesnt have enough cells, skipping\n', animals{a}, d);
            continue
        end
        
        
        if sum(ismember(digClass, northDigs)) < 2 || sum(ismember(digClass, southDigs)) < 2
            fprintf('Animal %s day %d does not have enough digs\n', animals{a}, d);
            continue
        end
        
        mdl = fitcsvm(rv, digClass', 'KernelFunction', 'linear',  'Prior', 'empirical', 'Leaveout', 'on', 'Verbose',0); 
        err = kfoldLoss(mdl, 'mode', 'individual');
        predAcc = 1 - err;
        perAnimalPred{a} = nanmean(predAcc);
    end
    cPred{d} = perAnimalPred;
end

avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Heading Prediction using FI cell rates\n')
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
yline(0.5);
hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('FI cell heading prediction');


%% Perform Heading prediction in FS cells. 
celltype = fs_ba;
%sessions = unique([celltype(:,1), celltype(:,end)], 'rows');
cPred = cell(1,3);
for d = 1:3
    dtbl = celltype(celltype.dayUsed == d,:);
    sessions = unique([dtbl(:,1), dtbl(:,end)], 'rows');
    animals = unique(sessions.animalName);
    perAnimalPred = cell(length(animals),1);
    for a = 1:length(animals)
        cellnames = celltype.animalSessionCellName(ismember(celltype.animalName, animals{a}) & celltype.dayUsed == d);
        stbl = mapsData(ismember(mapsData.animalName, animals{a}) & mapsData.dayUsed == d,:);
        digs = unique([stbl(:,6), stbl(:,9)], 'rows'); 
        tmprv = nan(length(digs.trialId), length(cellnames));
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
        validDigTrials = digs.trialId(ismember(digs.dig, [northDigs, southDigs]));
        digClass = categorical(digs.dig(ismember(digs.trialId, validDigTrials)));

        rv = tmprv(validDigTrials,:);
        
        if size(rv,2) < mincell
            fprintf('Animal %s day %d doesnt have enough cells, skipping\n', animals{a}, d);
            continue
        end
        
        % 
        if sum(ismember(digClass, northDigs)) < 2 || sum(ismember(digClass, southDigs)) < 2
            fprintf('Animal %s day %d does not have enough digs\n', animals{a}, d);
            continue
        end
        
        mdl = fitcsvm(rv, digClass', 'KernelFunction', 'linear',  'Prior', 'empirical', 'Leaveout', 'on', 'Verbose',0); 
        err = kfoldLoss(mdl, 'mode', 'individual');
        predAcc = 1 - err;
        perAnimalPred{a} = nanmean(predAcc);
    end
    cPred{d} = perAnimalPred;
end

avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Heading Prediction using FS cell rates\n')
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

yline(0.5)
hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('FS cell heading prediction');

%%

sgtitle('Heading Prediction per Animal');
fclose(fid);

saveas(gcf, fullfile(output_path, 'heading_pred_per_animal.pdf'))


