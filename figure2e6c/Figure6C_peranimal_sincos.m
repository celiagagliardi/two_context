%% per animal prediction across time
% use just the center out angle, if this doesn't work use sin and cosin as
% well

data_path = '/Users/celia/Documents/two_context_data';
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
maps = reformatTbl(calin.analysisInput.MapsData);
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
ba = reformatTbl(calout.analysisResults.BestAligned);
ficells = ba.animalSessionCellName(ba.isStable == 1);
tbl = maps(ismember(maps.animalSessionCellName, ficells),:);

panels = [1,2;1,3; 2, 3]; % col 1 is train day and col2 is predict day
northDigs = {'Corr'}; southDigs = {'Geo'};
animals = unique(tbl.animalName);

%%
for t = 1:size(panels,1)
    fprintf('Train Day %d, predict day %d\n', panels(t,1), panels(t,2));
    cPred = nan(length(animals),1);
    for a = 1:length(animals)
        atbl = tbl(ismember(tbl.animalName, animals{a}),:);
        fprintf('\tAnalysing %s\n', animals{a});
        cellnames = intersect(atbl.registeredCellName(atbl.dayUsed == panels(t,1)), atbl.registeredCellName(atbl.dayUsed == panels(t,2)));
        % train
        validData = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == panels(t,1),:);
        validTrials = unique(validData.trialId);
        digClass = cell(length(validTrials),1);
        for j = 1:length(validTrials)
            iTrial = validTrials(j);
            digClass(j) = unique(validData.dig(validData.trialId == iTrial));
        end
        digClass = categorical(digClass);
        rv = cell(length(digClass), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = validData(ismember(validData.registeredCellName, cellnames{c}),:);
            celltrials = ctbl.trialId;
            for j = 1:length(celltrials)
                iTrial = celltrials(j);
                tIdx = find(validTrials == iTrial);
                rv(tIdx,c) = ctbl.map(ctbl.trialId == iTrial);
            end
        end

        co_angles = cellfun(@centerOutRadians, rv);
        sinData = sin(co_angles);
        cosData = cos(co_angles);
        trainData = [sinData, cosData];

        nancols = any(isnan(trainData));
        trainData = trainData(:, ~nancols);
        mdl = fitcsvm(trainData, digClass, 'kernelfunction', 'linear', 'prior', 'empirical');

        % predict
        validData = atbl(ismember(atbl.dig, [northDigs, southDigs]) & atbl.dayUsed == panels(t,2),:);
        validTrials = unique(validData.trialId);
        digClass = cell(length(validTrials),1);
        for j = 1:length(validTrials)
            iTrial = validTrials(j);
            digClass(j) = unique(validData.dig(validData.trialId == iTrial));
        end
        testDig = categorical(digClass);
        rv = cell(length(testDig), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = validData(ismember(validData.registeredCellName, cellnames{c}),:);
            celltrials = ctbl.trialId;
            for j = 1:length(celltrials)
                iTrial = celltrials(j);
                tIdx = find(validTrials == iTrial);
                rv(tIdx,c) = ctbl.map(ctbl.trialId == iTrial);
            end
        end

        coangles = cellfun(@centerOutRadians, rv);
        sinData = sin(coangles);
        cosData = cos(coangles);
        testData = [sinData, cosData];
        testData = testData(:, ~nancols);
        label = predict(mdl, testData);
        pred_acc = mean(label == testDig);
        fprintf('\t\t %.3f accuracy\n', pred_acc);
        cPred(a) = pred_acc;
    end

   fprintf('overall prediction\n')
    avg = mean(cPred);
    sem = std(cPred) ./ sqrt(length(cPred));
    fprintf('\t\t avg = %.3f, sem = %.3f\n', avg, sem);

     h = kstest(cPred);
    if h
        fprintf('\t\tdata not normal, performing sign rank test\n');
        p = signrank(cPred, .5, 'tail', 'right');
        fprintf('\t\t p =  %.3f\n', p);
    else
        [~,p,ci, stats] = ttest(cPred, .5, 'tail', 'right');
        fprintf('One sample t test: t(%d) = %.3f, p = %.3f\n', stats.df, stats.tstat, p);
    end
end








        
    