%% Confirm heading prediction using data from same day and previous days

data_path = '/Users/celia/Documents/two_context_data';
tetin = load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_input.mat'));
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));

tetmaps = reformatTbl(tetin.analysisInput.MapsData);
calmaps = reformatTbl(calin.analysisInput.MapsData);

tetvars = tetmaps.Properties.VariableNames;
calvars = calmaps.Properties.VariableNames;
intvar = intersect(tetvars, calvars);

maps = [tetmaps(:, ismember(tetvars, intvar)); calmaps(:, ismember(calvars, intvar))];
northDigs = {'Corr'};
southDigs = {'Geo'};
%% Figure 2E
allanimals = unique(maps.animalName);
cPred = cell(length(allanimals),3);
for d = 1:3
    animals = unique(maps.animalName(maps.dayUsed == d));
    for a = 1:length(animals)
        atbl = maps(ismember(maps.animalName, animals{a}) & maps.dayUsed == d,:);
        vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]),:);
        if isempty(vtbl)
            continue
        end
        validTrials = unique(vtbl.trialId);
        if length(validTrials) < 2
            continue
        end
        allcells = unique(vtbl.animalSessionCellName);
       
        digClass = cell(length(validTrials), 1);
        for t = 1:length(validTrials)
            iTrial = validTrials(t);
            digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
        end
        digClass = categorical(digClass);
        
        rv = cell(length(validTrials), length(allcells));
        for c = 1:length(allcells)
            ctbl = vtbl(ismember(vtbl.animalSessionCellName, allcells{c}),:);
            celltrials = ctbl.trialId;
            for t = 1:length(celltrials)
                iTrial = celltrials(t);
                tIdx = find(validTrials == iTrial);
                rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
            end
        end

        co_angles = cellfun(@centerOutRadians, rv);
        nancols = any(isnan(co_angles));
        co_angles = co_angles(:,~nancols);
        mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf', 'Leaveout', 'on');
        err = kfoldLoss(mdl);
        cPred{a,d} = 1-err;
    end
end
avg = nan(1,3);
sem = nan(1,3);
for d = 1:3
    cdata = cell2mat(cPred(:,d));
    davg = nanmean(cdata);
    dsem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    avg(d) = davg;
    sem(d) = dsem;
    %p = signrank(cdata, .5, 'tail', 'right');
    [h,p] = ttest(cdata, .5, 'tail', 'right');
    fprintf('Day %d\n', d);
    fprintf('avg = %.3f, sem = %.3f, sign rank p = %.3f\n', davg, dsem, p);
end

figure;
bar(avg); hold on;
errorbar(avg, sem, 'k', 'linestyle', 'none');
title('Figure2E')
xlabel('day');
ylabel('heading pred acc')

%% Figure 6C
fprintf('Registered Predictions\n');
% left panel (train day 1, predict day 2)
trainDay = 1; predictDay = 2;
cellnames = intersect(calmaps.registeredCellName(calmaps.dayUsed == trainDay), ...
    calmaps.registeredCellName(calmaps.dayUsed == predictDay));
tbl = calmaps(ismember(calmaps.registeredCellName, cellnames),:);
animals = unique(tbl.animalName);

cPred = nan(length(animals),1);
for a = 1:length(animals)
    %train
    atbl = tbl(ismember(tbl.animalName, animals{a}),:);
    vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == trainDay,:);
    validTrials = unique(vtbl.trialId);
    allcells = unique(vtbl.registeredCellName);

    digClass = cell(length(validTrials), 1);
    for t = 1:length(validTrials)
        iTrial = validTrials(t);
        digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
    end
    digClass = categorical(digClass);
    rv = cell(length(validTrials), length(allcells));
    for c = 1:length(allcells)
        ctbl = vtbl(ismember(vtbl.registeredCellName, allcells{c}),:);
        celltrials = ctbl.trialId;
        for t = 1:length(celltrials)
            iTrial = celltrials(t);
            tIdx = find(validTrials == iTrial);
            rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
        end
    end

    co_angles = cellfun(@centerOutRadians, rv);
    nancols = any(isnan(co_angles));
    co_angles = co_angles(:,~nancols);
    %mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf', 'Leaveout', 'on');
   mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf');

   % predidct
    vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == predictDay,:);
    validTrials = unique(vtbl.trialId);
     digClass = cell(length(validTrials), 1);
    for t = 1:length(validTrials)
        iTrial = validTrials(t);
        digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
    end
    digClass = categorical(digClass);
    rv = cell(length(validTrials), length(allcells));
    for c = 1:length(allcells)
        ctbl = vtbl(ismember(vtbl.registeredCellName, allcells{c}),:);
        celltrials = ctbl.trialId;
        for t = 1:length(celltrials)
            iTrial = celltrials(t);
            tIdx = find(validTrials == iTrial);
            rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
        end
    end
     co_angles = cellfun(@centerOutRadians, rv);
     co_angles = co_angles(:,~nancols);

     cmdl = compact(mdl);
     predictedDigs = predict(cmdl, co_angles);
     cPred(a) = mean(predictedDigs == digClass);
end
[sortedPred, sortIdx] = sort(cPred, 'ascend');
sAnimals = animals(sortIdx);
overall = mean(cPred);
%[h,p] = ttest(cPred, .5, 'tail', 'right');
[p,h] = signrank(cPred, .5, 'tail', 'right');
fprintf('Train %d Predict %d results\n', trainDay, predictDay)
fprintf('\t sign rank test h = %d, p = %.3f\n', h,p);
acc = [overall; sortedPred];
figure
barh(acc); hold on;
xlabel('Prediction Accuracy');
yticklabels(['overall'; sAnimals]);
xlim([0 1]);
xline(.5);
title(sprintf('Train %d, Predict %d', trainDay, predictDay));
hold off;

% middle panel (train day 2, predict day 3)
trainDay = 2; predictDay = 3;
cellnames = intersect(calmaps.registeredCellName(calmaps.dayUsed == trainDay), ...
    calmaps.registeredCellName(calmaps.dayUsed == predictDay));
tbl = calmaps(ismember(calmaps.registeredCellName, cellnames),:);
animals = unique(tbl.animalName);

cPred = nan(length(animals),1);
for a = 1:length(animals)
    %train
    atbl = tbl(ismember(tbl.animalName, animals{a}),:);
    vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == trainDay,:);
    validTrials = unique(vtbl.trialId);
    allcells = unique(vtbl.registeredCellName);

    digClass = cell(length(validTrials), 1);
    for t = 1:length(validTrials)
        iTrial = validTrials(t);
        digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
    end
    digClass = categorical(digClass);
    rv = cell(length(validTrials), length(allcells));
    for c = 1:length(allcells)
        ctbl = vtbl(ismember(vtbl.registeredCellName, allcells{c}),:);
        celltrials = ctbl.trialId;
        for t = 1:length(celltrials)
            iTrial = celltrials(t);
            tIdx = find(validTrials == iTrial);
            rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
        end
    end

    co_angles = cellfun(@centerOutRadians, rv);
    nancols = any(isnan(co_angles));
    co_angles = co_angles(:,~nancols);
    %mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf', 'Leaveout', 'on');
   mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf');

   % predidct
    vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == predictDay,:);
    validTrials = unique(vtbl.trialId);
     digClass = cell(length(validTrials), 1);
    for t = 1:length(validTrials)
        iTrial = validTrials(t);
        digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
    end
    digClass = categorical(digClass);
    rv = cell(length(validTrials), length(allcells));
    for c = 1:length(allcells)
        ctbl = vtbl(ismember(vtbl.registeredCellName, allcells{c}),:);
        celltrials = ctbl.trialId;
        for t = 1:length(celltrials)
            iTrial = celltrials(t);
            tIdx = find(validTrials == iTrial);
            rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
        end
    end
     co_angles = cellfun(@centerOutRadians, rv);
     co_angles = co_angles(:,~nancols);

     cmdl = compact(mdl);
     predictedDigs = predict(cmdl, co_angles);
     cPred(a) = mean(predictedDigs == digClass);
end
[sortedPred, sortIdx] = sort(cPred, 'ascend');
sAnimals = animals(sortIdx);
overall = mean(cPred);
%[h,p] = ttest(cPred, .5, 'tail', 'right');
[p,h] = signrank(cPred, .5, 'tail', 'right');
fprintf('Train %d Predict %d results\n', trainDay, predictDay)
fprintf('\t sign rank test h = %d, p = %.3f\n', h,p);
acc = [overall; sortedPred];

figure
barh(acc); hold on;
xlabel('Prediction Accuracy');
yticklabels(['overall'; sAnimals]);
xlim([0 1]);
xline(.5);
title(sprintf('Train %d, Predict %d', trainDay, predictDay));
hold off;

% right panel (train day 1, predict day 3)
trainDay = 1; predictDay = 3;
cellnames = intersect(calmaps.registeredCellName(calmaps.dayUsed == trainDay), ...
    calmaps.registeredCellName(calmaps.dayUsed == predictDay));
tbl = calmaps(ismember(calmaps.registeredCellName, cellnames),:);
animals = unique(tbl.animalName);

cPred = nan(length(animals),1);
for a = 1:length(animals)
    %train
    atbl = tbl(ismember(tbl.animalName, animals{a}),:);
    vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == trainDay,:);
    validTrials = unique(vtbl.trialId);
    allcells = unique(vtbl.registeredCellName);

    digClass = cell(length(validTrials), 1);
    for t = 1:length(validTrials)
        iTrial = validTrials(t);
        digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
    end
    digClass = categorical(digClass);
    rv = cell(length(validTrials), length(allcells));
    for c = 1:length(allcells)
        ctbl = vtbl(ismember(vtbl.registeredCellName, allcells{c}),:);
        celltrials = ctbl.trialId;
        for t = 1:length(celltrials)
            iTrial = celltrials(t);
            tIdx = find(validTrials == iTrial);
            rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
        end
    end

    co_angles = cellfun(@centerOutRadians, rv);
    nancols = any(isnan(co_angles));
    co_angles = co_angles(:,~nancols);
    %mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf', 'Leaveout', 'on');
   mdl = fitcsvm(co_angles, digClass, 'kernelfunction', 'rbf');

   % predidct
    vtbl = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == predictDay,:);
    validTrials = unique(vtbl.trialId);
     digClass = cell(length(validTrials), 1);
    for t = 1:length(validTrials)
        iTrial = validTrials(t);
        digClass(t) = unique(vtbl.dig(vtbl.trialId == iTrial));
    end
    digClass = categorical(digClass);
    rv = cell(length(validTrials), length(allcells));
    for c = 1:length(allcells)
        ctbl = vtbl(ismember(vtbl.registeredCellName, allcells{c}),:);
        celltrials = ctbl.trialId;
        for t = 1:length(celltrials)
            iTrial = celltrials(t);
            tIdx = find(validTrials == iTrial);
            rv(tIdx, c) = ctbl.map(ctbl.trialId == iTrial);
        end
    end
     co_angles = cellfun(@centerOutRadians, rv);
     co_angles = co_angles(:,~nancols);

     cmdl = compact(mdl);
     predictedDigs = predict(cmdl, co_angles);
     cPred(a) = mean(predictedDigs == digClass);
end
[sortedPred, sortIdx] = sort(cPred, 'ascend');
sAnimals = animals(sortIdx);
overall = mean(cPred);
%[h,p] = ttest(cPred, .5, 'tail', 'right');
[p,h] = signrank(cPred, .5, 'tail', 'right');
fprintf('Train %d Predict %d results\n', trainDay, predictDay)
fprintf('\t signrank test: h = %d, p = %.3f\n', h,p);
acc = [overall; sortedPred];
figure
barh(acc); hold on;
xlabel('Prediction Accuracy');
yticklabels(['overall'; sAnimals]);
xlim([0 1]);
xline(.5);
title(sprintf('Train %d, Predict %d', trainDay, predictDay));
hold off;


