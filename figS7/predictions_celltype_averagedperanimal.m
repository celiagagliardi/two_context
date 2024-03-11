%% Predict heading and context using only CG trials 
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

FIheading = cell(7,3);
FSheading = cell(7,3);
FIcontext = cell(7,3);
FScontext = cell(7,3);

%% compute

for d = 1:3
    animals = unique(baData.animalName(baData.dayUsed == d));
    for a = 1:length(animals)
        celltype = 1; % 1 = FI, 0 = FS;
        allcells = baData.animalSessionCellName(baData.dayUsed == d & ismember(baData.animalName, animals{a}) &  baData.isStable == celltype,:);
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

        FIheading{a,d} = heading_prediction;
        FIcontext{a,d} = context_prediction;

        % FS cells
        celltype = 0; % 1 = FI, 0 = FS;
        allcells = baData.animalSessionCellName(baData.dayUsed == d & ismember(baData.animalName, animals{a}) &  baData.isStable == celltype,:);
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

        FSheading{a,d} = heading_prediction;
        FScontext{a,d} = context_prediction;
    end
end

%% Stats
fiheadingavg = cellfun(@nanmean, FIheading);
fsheadingavg = cellfun(@nanmean, FSheading);
ficontextavg = cellfun(@nanmean, FIcontext);
fscontextavg = cellfun(@nanmean, FScontext);

digavgFI = nan(1,3);
digsemFI = nan(1,3);
conavgFI = nan(1,3);
consemFI = nan(1,3);
digavgFS = nan(1,3);
digsemFS = nan(1,3);
conavgFS = nan(1,3);
consemFS = nan(1,3);
% sign rank test
for d = 1:3
    fprintf('Day %d FI cells: \n', d);
    cdata = fiheadingavg(:,d);
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\t heading avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg, sem, h, p, sum(~isnan(cdata)));
    digavgFI(d) = avg; digsemFI(d) = sem;

    cdata = ficontextavg(:,d);
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\tcontext avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg, sem, h, p, sum(~isnan(cdata)));
    conavgFI(d) = avg; consemFI(d) = sem;

    fprintf('Day %d FS cells:\n', d);
    cdata = fsheadingavg(:,d);
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\theading avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg, sem, h, p, sum(~isnan(cdata)));
    digavgFS(d) = avg; digsemFS(d) = sem;

    cdata = fscontextavg(:,d);
    avg = nanmean(cdata);
    sem = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    [p,h] = signrank(cdata, .5, 'tail', 'right');
    fprintf('\tcontext avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg, sem, h, p, sum(~isnan(cdata)));
    conavgFS(d) = avg; consemFS(d) = sem;
end

 figure;
    subplot(2,2,1);
    bar(digavgFI);
    hold on;
    errorbar(digavgFI, digsemFI, 'k--', 'linestyle', 'none');
    yline(.5)
    ylim([0 1]);
    hold off;
    title('FI Heading Prediction');

    subplot(2,2,2)
    bar(conavgFI);
    hold on;
    errorbar(conavgFI, consemFI, 'k--',  'linestyle', 'none');
    yline(.5)
    ylim([0 1]);
    hold off;
    title('FI Context Prediction');

    subplot(2,2,3);
    bar(digavgFS);
    hold on;
    errorbar(digavgFS, digsemFS, 'k--', 'linestyle', 'none');
    yline(.5)
    ylim([0 1]);
    hold off;
    title('FS Heading Prediction');

    subplot(2,2,4)
    bar(conavgFS);
    hold on;
    errorbar(conavgFS, consemFS, 'k--', 'linestyle', 'none');
    yline(.5)
    ylim([0 1]);
    hold off;
    title('FS Context Prediction');

%% FI vs FS
