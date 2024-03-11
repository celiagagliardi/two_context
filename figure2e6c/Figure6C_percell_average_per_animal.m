clear all; clc;

data_path = '/Users/celia/Documents/two_context_data';
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
maps = reformatTbl(calin.analysisInput.MapsData);
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
ba = reformatTbl(calout.analysisResults.BestAligned);
ficells = ba.animalSessionCellName(ba.isStable == 1);
%tbl = maps(ismember(maps.animalSessionCellName, ficells),:);
tbl = maps;
panels = [1,2;1,3; 2, 3]; % col 1 is train day and col2 is predict day
northDigs = {'Corr'}; southDigs = {'Geo'};

%%
figure
for t = 1:size(panels,1)
    fprintf('Train Day %d, predict day %d\n', panels(t,1), panels(t,2));
    cellnames = intersect(tbl.registeredCellName(tbl.dayUsed == panels(t,1)), tbl.registeredCellName(tbl.dayUsed == panels(t,2)));
    animals = unique(tbl.animalName(ismember(tbl.registeredCellName, cellnames)));
    cPred = nan(length(animals),1);
   
    for a = 1:length(animals)
        atbl = tbl(ismember(tbl.animalName, animals{a}) & ismember(tbl.registeredCellName, cellnames),:);
        acells = unique(atbl.registeredCellName);
        cell_pred = nan(length(acells),1);
        for c = 1:length(acells)
            % train
            trainData = atbl(ismember(atbl.registeredCellName, acells{c}) & atbl.dayUsed == panels(t,1),:);
            validData = trainData(ismember(trainData.dig, [northDigs, southDigs]), :);
            if isempty(validData)
                fprintf('not enough training data, skipping\n');
                continue
            end
            digClass = categorical(validData.dig);
            digClass(ismember(digClass, northDigs)) = 'North';
            digClass(ismember(digClass, southDigs)) = 'South';
            rv = cellfun(@centerOutRadians, validData.map);
            
            if length(unique(digClass(~isnan(rv))))  < 2
                continue
            end
            
           trainVector = [sin(rv), cos(rv)];
            %trainVector = rv;
            mdl = fitcsvm(trainVector, digClass); %,'standardize', 'on', 'kernelfunction', 'linear'); %% NO CROSS VALIDATION
            %pmdl = fitPosterior(mdl, trainVector, digClass);
            cmdl = crossval(mdl, 'leaveout', 'on');

            % predict
            testData = atbl(ismember(atbl.registeredCellName, acells{c}) & atbl.dayUsed == panels(t,2),:);
            validData = testData(ismember(testData.dig, [northDigs, southDigs]), :);
            testDig = categorical(validData.dig);
            testDig(ismember(testDig, northDigs)) = 'North';
            testDig(ismember(testDig, southDigs)) = 'South';
            rv = cellfun(@centerOutRadians, validData.map);
            testVector = [sin(rv), cos(rv)];
          % testVector = rv;

            %label = predict(mdl, testVector);
            %label = predict(pmdl, testVector);
            pred_acc = mean(testDig == label);
            cell_pred(c) = pred_acc;
        end
        cPred(a) = nanmean(cell_pred);
    end

    avg = nanmean(cPred);
    sem = std(cPred, 'omitnan') ./ sqrt(sum(~isnan(cPred)));
    fprintf('\t avg = %.3f, sem = %.3f\n', avg, sem);

     p = signrank(cPred, .5, 'tail', 'right');
     fprintf('\t\t sign rank test: p =  %.3f\n', p);
      [~,p,ci, stats] = ttest(cPred, .5, 'tail', 'right');
        fprintf('\t\tOne sample t test: t(%d) = %.3f, p = %.3f\n', stats.df, stats.tstat, p);

    h = kstest(cPred);
    if h
        fprintf('\n\t\tdata not normal, recommend performing sign rank test\n\n');  
    end

     subplot(1,3,t)
    ynames = [animals; {'overall'}];
    xvals = [cPred; avg];
    barh(xvals); hold on;
    yticklabels(ynames);
    xline(.5);
    xlim([0 1]);
    title(sprintf('Train Day %d, Predict Day %d', panels(t,1), panels(t,2)));
end

sgtitle('Per cell averaged per animal')


