%% Fig3E
% Examples of FI and FS cells

data_path = '/Users/celia/Documents/two_context_data'; % change this
cal = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
tet = load(fullfile(data_path, 'tetrodes_20220316_100902/analysis_results.mat'));

best_aligned = [cal.analysisResults.BestAligned; tet.analysisResults.BestAligned];

%%
% ephys, FI cells
excells = {'AK42_CA1_d7_TT2_04.t', 'AK42_CA1_d7_TT5_01.t','AK74_CA1_d3_TT4_01.t'}
for x = 1:length(excells)
    ind = ismember(best_aligned.animalSessionCellName, excells{x});
    figure
    subplot(1,2,1)
    pcolor(best_aligned.averageMapContext1{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;
    subplot(1,2,2);
    pcolor(best_aligned.averageMapContext2{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;

    sgtitle(sprintf('r = %.3f', best_aligned.bestCorrelation(ind)));
end
% calcium, FI cells
excells = {'CMG162_CA1_s1_115.t', 'CMG154_CA1_s1_151.t', 'CMG154_CA1_s1_110.t'};
for x = 1:length(excells)
    ind = ismember(best_aligned.cellName, excells{x});
    figure
    subplot(1,2,1)
    pcolor(best_aligned.averageMapContext1{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;
    subplot(1,2,2);
    pcolor(best_aligned.averageMapContext2{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;

    sgtitle(sprintf('r = %.3f', best_aligned.bestCorrelation(ind)));
end

%%
% ephys, FS cells
excells = {'AK74_CA1_d2_TT4_03.t', 'AK74_CA1_d2_TT4_01.t','JJ9_CA1_d3_TT5_05.t'};
for x = 1:length(excells)
    ind = ismember(best_aligned.animalSessionCellName, excells{x});
    figure
    subplot(1,2,1)
    pcolor(best_aligned.averageMapContext1{ind});
    daspect([1 1 1]); colormap jet; shading interp;  axis ij;
    subplot(1,2,2);
    pcolor(best_aligned.averageMapContext2{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;

    sgtitle(sprintf('r = %.3f', best_aligned.bestCorrelation(ind)));
end

% calcium, FS cells
excells = {'CMG154_CA1_s1_304.t', 'CMG162_CA1_s1_174.t', 'CMG162_CA1_s5_24.t'};
for x = 1:length(excells)
    ind = ismember(best_aligned.cellName, excells{x});
    figure
    subplot(1,2,1)
    pcolor(best_aligned.averageMapContext1{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;
    subplot(1,2,2);
    pcolor(best_aligned.averageMapContext2{ind});
    daspect([1 1 1]); colormap jet; shading interp; axis ij;

    sgtitle(sprintf('r = %.3f', best_aligned.bestCorrelation(ind)));
end
