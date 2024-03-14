

clear all; clc;
%% load data
data_path = '/Users/celia/Documents/two_context_data'; % change this
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tetout = load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_results.mat'));
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
tetin = load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_input.mat'));

%% combine maps from tetrodes and calcium
tetmaps = reformatTbl(tetin.analysisInput.MapsData);
calmaps = reformatTbl(calin.analysisInput.MapsData);
same_vars = intersect(tetmaps.Properties.VariableNames, calmaps.Properties.VariableNames);
maps = [tetmaps(:, ismember(tetmaps.Properties.VariableNames, same_vars)); ...
    calmaps(:, ismember(calmaps.Properties.VariableNames, same_vars))];

% combine best aligned data to identify FI and FS cells
ba = [reformatTbl(tetout.analysisResults.BestAligned); reformatTbl(calout.analysisResults.BestAligned)];
ficells_all = ba.animalSessionCellName(ba.isStable == 1);
fscells_all = ba.animalSessionCellName(ba.isStable == 0);

animals = unique(ba.animalName); % its 12.
fi_pred = cell(length(animals), 3);
fs_pred = cell(length(animals), 3);

%% compute
for a = 1:length(animals)
    for d = 1:3
        tbl = maps(ismember(maps.animalName, animals{a}) & maps.dayUsed == d,:);
        if isempty(tbl)
            continue
        end
        mapSize = size(tbl.map{1});
        ficells = unique(tbl.animalSessionCellName(ismember(tbl.animalSessionCellName, ficells_all)));
        fscells = unique(tbl.animalSessionCellName(ismember(tbl.animalSessionCellName, fscells_all)));

        fprintf('Animal %s day %d has %d FI and %d FS cells\n', animals{a}, d, length(ficells), length(fscells));
        contextTrial = unique([tbl.trialId, tbl.contextId], 'rows');
        trialId = contextTrial(:,1);
        contextId = contextTrial(:,2);

        fi_maps = cell(length(ficells), length(trialId));
        fs_maps = cell(length(fscells), length(trialId));

        %% collect maps and align them
        % fi maps
        for c = 1:length(ficells)
            ctbl = tbl(ismember(tbl.animalSessionCellName, ficells{c}),:);
            cellTrials = unique(ctbl.trialId);
            cba = ba(ismember(ba.animalSessionCellName, ficells{c}),:);
            c1Trials = ctbl.trialId(ctbl.contextId == 1);
            c2Trials = ctbl.trialId(ctbl.contextId == 2);
            rot1 = cell2mat(cba.rotationSequence1);
            rot2 = cell2mat(cba.rotationSequence2);

            for cIdx = 1:length(c1Trials)
                iTrial = c1Trials(cIdx);
                trialIdx = find(iTrial == trialId);
                cellmap = ctbl.map{ctbl.trialId == iTrial};
                if rot1(cIdx)
                    fi_maps{c,trialIdx} = rot90(cellmap,2);
                else
                    fi_maps{c,trialIdx} = cellmap;
                end
            end
            for cIdx = 1:length(c2Trials)
                iTrial = c2Trials(cIdx);
                trialIdx = find(iTrial == trialId);
                cellmap = ctbl.map{ctbl.trialId == iTrial};
                if rot2(cIdx)
                    fi_maps{c,trialIdx} = rot90(cellmap,2);
                else
                    fi_maps{c,trialIdx} = cellmap;
                end
            end
        end

        % collect and align FS cells
        for c = 1:length(fscells)
            ctbl = tbl(ismember(tbl.animalSessionCellName, fscells{c}),:);
            cellTrials = unique(ctbl.trialId);
            cba = ba(ismember(ba.animalSessionCellName, fscells{c}),:);
            c1Trials = ctbl.trialId(ctbl.contextId == 1);
            c2Trials = ctbl.trialId(ctbl.contextId == 2);
            rot1 = cell2mat(cba.rotationSequence1);
            rot2 = cell2mat(cba.rotationSequence2);

            for cIdx = 1:length(c1Trials)
                iTrial = c1Trials(cIdx);
                trialIdx = find(iTrial == trialId);
                cellmap = ctbl.map{ctbl.trialId == iTrial};
                if rot1(cIdx)
                    fs_maps{c,trialIdx} = rot90(cellmap,2);
                else
                    fs_maps{c,trialIdx} = cellmap;
                end
            end
            for cIdx = 1:length(c2Trials)
                iTrial = c2Trials(cIdx);
                trialIdx = find(iTrial == trialId);
                cellmap = ctbl.map{ctbl.trialId == iTrial};
                if rot2(cIdx)
                    fs_maps{c,trialIdx} = rot90(cellmap,2);
                else
                    fs_maps{c,trialIdx} = cellmap;
                end
            end
        end

%  t
        % put nans in empty trial maps
        fi_maps(cellfun(@isempty, fi_maps)) = {nan(mapSize)};
        fs_maps(cellfun(@isempty, fs_maps)) = {nan(mapSize)};

        %% averages and predict one cell at a time
        fi_predict = nan(length(trialId),1);
        fs_predict = nan(length(trialId),1);
        if ~isempty(fi_maps)
            for w = 1:size(fi_maps,2)
                wMaps = cat(3,fi_maps{:,w});
                c1Trials = find(contextId == 1);
                c2Trials = find(contextId == 2);
                wC1 = c1Trials(c1Trials ~= w);
                wC2 = c2Trials(c2Trials ~= w);
                c1maps = fi_maps(:, wC1);
                c2maps = fi_maps(:,wC2);
                % build average maps
               
                avgC1 = nan([mapSize, length(ficells)]);
                avgC2 = nan([mapSize, length(ficells)]);
                for c = 1:length(ficells)
                    cellmaps = c1maps(c,:);
                    catmaps = cat(3, cellmaps{:});
                    if ~isempty(catmaps)
                        avgC1(:,:,c) = nanmean(catmaps,3);
                    end

                    cellmaps = c2maps(c,:);
                    catmaps = cat(3, cellmaps{:});
                    if ~isempty(catmaps)
                        avgC2(:,:,c) = nanmean(catmaps,3);
                    end
                end

                wValid = ~squeeze(any(isnan(wMaps), [1 2]));
                c1Valid = ~squeeze(any(isnan(avgC1), [1,2]));
                c2Valid = ~squeeze(any(isnan(avgC2), [1,2]));

                if ~any(wValid) || ~any(c1Valid) || ~any(c2Valid)
                    fprintf('%s day %d:\5either withheld, avg c1, or avg c2 maps are all nans, skipping\n', animals{a}, d);
                    continue;
                end

                % only compute with non nans
                validCells = all([wValid, c1Valid],2);
                [c1Dp, c1DpAvg] = ml_alg_popvectors_compute_dotproducts(wMaps(:,:,validCells), avgC1(:,:,validCells));
                validCells = all([wValid, c2Valid],2);
                [c2Dp, c2DpAvg] = ml_alg_popvectors_compute_dotproducts(wMaps(:,:,validCells), avgC2(:,:,validCells));

                if c1DpAvg > c2DpAvg
                    predTrial = 1;
                elseif c2DpAvg > c1DpAvg
                    predTrial = 2;
                else
                    % this indicates only one cell had enough maps to
                    % predict using dot product, use correlation for these
                    % cases
                    validCells = all([wValid, c1Valid],2);
                    c1corr = corr2(wMaps(:,:,validCells), avgC1(:,:,validCells));
                    validCells = all([wValid, c2Valid],2);
                    c2corr = corr2(wMaps(:,:,validCells), avgC2(:,:,validCells));

                    if c1corr > c2corr
                        predTrial = 1;
                    elseif c2corr > c1corr
                        predTrial = 2;
                    else
                        fprintf('%s day %d still showing nans\n', animals{a}, d);
                        predTrial = nan;
                    end

                end

                fi_predict(w) = ismember(contextId(w), predTrial);
            end
        end

        fi_pred{a,d} = nanmean(fi_predict);

        % now FS maps
        if ~isempty(fs_maps)
            for w = 1:size(fs_maps,2)
                wMaps = cat(3,fs_maps{:,w});
                c1Trials = find(contextId == 1);
                c2Trials = find(contextId == 2);
                wC1 = c1Trials(c1Trials ~= w);
                wC2 = c2Trials(c2Trials ~= w);
                c1maps = fs_maps(:,wC1);
                c2maps = fs_maps(:,wC2);

                % build average maps
                avgC1 = nan([mapSize, length(fscells)]);
                avgC2 = nan([mapSize, length(fscells)]);
                for c = 1:length(fscells)
                    cellmaps = c1maps(c,:);
                    catmaps = cat(3, cellmaps{:});
                    if ~isempty(catmaps)
                        avgC1(:,:,c) = nanmean(catmaps,3);
                    end

                    cellmaps = c2maps(c,:);
                    catmaps = cat(3, cellmaps{:});
                    if ~isempty(catmaps)
                        avgC2(:,:,c) = nanmean(catmaps,3);
                    end
                end

                wValid = ~squeeze(any(isnan(wMaps), [1 2]));
                c1Valid = ~squeeze(any(isnan(avgC1), [1,2]));
                c2Valid = ~squeeze(any(isnan(avgC2), [1,2]));

                if ~any(wValid) || ~any(c1Valid) || ~any(c2Valid)
                    fprintf('%s day %d:\5either withheld, avg c1, or avg c2 maps are all nans, skipping\n', animals{a}, d);
                    continue;
                end

                % only compute with non nans
                validCells = all([wValid, c1Valid],2);
                [c1Dp, c1DpAvg] = ml_alg_popvectors_compute_dotproducts(wMaps(:,:,validCells), avgC1(:,:,validCells));
                validCells = all([wValid, c2Valid],2);
                [c2Dp, c2DpAvg] = ml_alg_popvectors_compute_dotproducts(wMaps(:,:,validCells), avgC2(:,:,validCells));

                if c1DpAvg > c2DpAvg
                    predTrial = 1;
                elseif c2DpAvg > c1DpAvg
                    predTrial = 2;
                else
                    % this indicates only one cell had enough maps to
                    % predict using dot product, use correlation for these
                    % cases

                    validCells = all([wValid, c1Valid],2);
                    c1corr = corr2(wMaps(:,:,validCells), avgC1(:,:,validCells));
                    validCells = all([wValid, c2Valid],2);
                    c2corr = corr2(wMaps(:,:,validCells), avgC2(:,:,validCells));
                    if c1corr > c2corr
                        predTrial = 1;
                    elseif c2corr > c1corr
                        predTrial = 2;
                    else
                        fprintf('%s day %d still showing nans\n', animals{a}, d);
                        predTrial = nan;
                    end

                end

                fs_predict(w) = ismember(contextId(w), predTrial);
            end
        end

        fs_pred{a,d} = nanmean(fs_predict);
    end
end

%% Graph results
avg = nan(3,1);
sem = nan(3,1);
fprintf('Context Prediction using average dot product of population vectors OR CORRELATION \n')
for d = 1:3
    fprintf('Day %d FI cells\n', d)
    cdata = cell2mat(fi_pred(:,d));
    avg(d) = mean(cdata, 'omitnan');
    sem(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));

    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');

    %fprintf('\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f,  using %s cells and %d %ss\n',  avg(d), sem(d), stats.df, stats.tstat, p, 'FI', sum(~isnan(cdata)), 'animals');
    fprintf('\t avg = %.4f, sem = %.4f\n',  avg(d), sem(d));
    fprintf('\t\t one sample ttest: t(%d) = %.3f, p = %.5f, using %s cells and %d %s\n', stats.df, stats.tstat, p, 'FI', sum(~isnan(cdata)), 'animals');
    p = signrank(cdata, .5, 'tail', 'right');
    fprintf('\t\t signrank test: p = %.3f\n', p);

end

% Example data
%figure;
figure
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
yline(0.5, 'k', 'linewidth', 2);
scatter(1, cell2mat(fi_pred(:,1)), 'ko', 'jitter', 'on');
scatter(2, cell2mat(fi_pred(:,2)), 'ro', 'jitter', 'on');
scatter(3, cell2mat(fi_pred(:,3)), 'bo', 'jitter', 'on');

hold off
ylim([0 1])
title('FI cells')
%% Graph FS results

avg = nan(3,1);
sem = nan(3,1);
fprintf('Context Prediction using average dot product of population vectors OR CORRELATION \n')
for d = 1:3
    fprintf('Day %d FS cells\n', d)
    cdata = cell2mat(fs_pred(:,d));
    avg(d) = mean(cdata, 'omitnan');
    sem(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));

    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');

    %fprintf('\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f,  using %s cells and %d %ss\n',  avg(d), sem(d), stats.df, stats.tstat, p, 'FI', sum(~isnan(cdata)), 'animals');
    fprintf('\t avg = %.4f, sem = %.4f, using %s cells and %d %s\n',  avg(d), sem(d), 'FS', sum(~isnan(cdata)), 'animals');
    fprintf('\t\t one tail t test: t(%d) = %.3f, p = %.5f\n',  stats.df, stats.tstat, p);
    p = signrank(cdata, .5, 'tail', 'right');
    fprintf('\t\t sign rank test: p = %.3f\n', p);
end


% Example data
%figure;
figure
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
yline(0.5, 'k', 'linewidth', 2);
scatter(1, cell2mat(fs_pred(:,1)), 'ko', 'jitter', 'on');
scatter(2, cell2mat(fs_pred(:,2)), 'ro', 'jitter', 'on');
scatter(3, cell2mat(fs_pred(:,3)), 'bo', 'jitter', 'on');
title('FS context prediction');

hold off
ylim([0 1])
