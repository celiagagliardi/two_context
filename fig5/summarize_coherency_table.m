function out = summarize_coherency_table(fn)
if nargin < 1
    fn = uigetfile('output/*_table.csv');
end

fprintf('Performing summary of %s output\n', fn);
outputFolder = 'output';
if ~isfolder(outputFolder)
    mkdir(outputFolder)
end

filename = 'coherency_summary.csv';

fid = fopen(fullfile(outputFolder, filename), 'w');
fprintf(fid, [repmat('%s,', 1, 16), '%s\n'], 'animal', 'day', 'trialA', 'trialB', ...
    'contextA', 'contextB', 'comparison', 'FIisCoherent', 'FI_rot1', 'FI_rot2', ...
    'FI_rot3', 'FI_rot4', 'FSisCoherent', 'FS_rot1', 'FS_rot2', 'FS_rot3', 'FS_rot4');
formatSpec = '%s, %d, %d, %d, %d, %d, %s, %d, %.3f, %.3f, %.3f, %.3f, %d, %.3f, %.3f, %.3f, %.3f\n';

minCell = 0;

data = readtable(fullfile(outputFolder, fn));
animals = unique(data.animal);

fi_within_all = cell(length(animals),3);
fi_across_all = cell(length(animals),3);
fs_within_all = cell(length(animals),3);
fs_across_all = cell(length(animals),3);

fi_within_cellpair = cell(length(animals),3);
fi_across_cellpair = cell(length(animals),3);
fs_within_cellpair = cell(length(animals),3);
fs_across_cellpair = cell(length(animals),3);

for a = 1:length(animals)
    for d = 1:3
        dtbl = data(ismember(data.animal, animals{a}) & data.day == d,:);
        pwTrials = unique([dtbl.trialA, dtbl.trialB], 'rows');
        contextChanges = [mod(pwTrials(:,1),2), mod(pwTrials(:,2),2)];
        withinContext = contextChanges(:,1)== contextChanges(:,2);

        FIpwCoh = nan(size(pwTrials,1),1);
        FSpwCoh = nan(size(pwTrials,1),1);
        fi_rotdist = zeros(size(pwTrials,1),4);
        fs_rotdist = zeros(size(pwTrials,1),4);
    
        for pw = 1:size(pwTrials,1)
            pwTbl = dtbl(dtbl.trialA == pwTrials(pw,1) & dtbl.trialB == pwTrials(pw,2),:);
            
            
            contextA = unique(pwTbl.contextA);
            contextB = unique(pwTbl.contextB);
            
            if contextA == contextB
                comparison = 'within';
            else
                comparison = 'across';
            end

            % FI CELLS

            bmr = pwTbl.bmr(~isnan(pwTbl.bmr) & pwTbl.isStable == 1);
            numCells = sum(~isnan(pwTbl.bmr) & pwTbl.isStable == 1);
            
            [pval] = calc_chi2(bmr, 1:5);
            if pval < 0.05
                Fih = 1;
            elseif isnan(pval)
                Fih = nan;
            else
                Fih = 0;
            end
            if numCells < minCell
                FIpwCoh(pw) = nan;
            else
                FIpwCoh(pw) = Fih;
            end

            n = histcounts(bmr, 1:5) / numCells;
            sort_n = sort(n, 'descend') ;
            fi_rotdist(pw,:) = sort_n;
            fi_rot1 = sort_n(1); fi_rot2 = sort_n(2); fi_rot3 = sort_n(3); fi_rot4 = sort_n(4);

            % FS CELLS
            bmr = pwTbl.bmr(~isnan(pwTbl.bmr) & pwTbl.isStable == 0);
            numCells = sum(~isnan(pwTbl.bmr) & pwTbl.isStable == 0);

            pval = calc_chi2(bmr, 1:5);
            if pval < 0.05
                Fsh = 1;
            elseif isnan(pval)
                Fsh = nan;
            else
                Fsh = 0;
            end

            if numCells < minCell
                FSpwCoh(pw) = nan;
            else
                FSpwCoh(pw) = Fsh;
            end

            n = histcounts(bmr, 1:5)/ numCells;
            sort_n = sort(n, 'descend');
            fs_rotdist(pw,:) = sort_n;
            fs_rot1 = sort_n(1); fs_rot2 = sort_n(2); fs_rot3 = sort_n(3); fs_rot4 = sort_n(4);



           fprintf(fid, formatSpec, animals{a}, d, pwTrials(pw,1), pwTrials(pw,2), ...
               contextA, contextB, comparison, Fih, fi_rot1, fi_rot2, fi_rot3, fi_rot4, ...
               Fsh, fs_rot1, fs_rot2, fs_rot3, fs_rot4);


        end
        fi_within_all{a,d} = FIpwCoh(withinContext);
        fi_across_all{a,d} = FIpwCoh(~withinContext);
        fi_within_cellpair{a,d} = fi_rotdist(withinContext);
        fi_across_cellpair{a,d} = fi_rotdist(~withinContext);

        fs_within_all{a,d} = FSpwCoh(withinContext);
        fs_across_all{a,d} = FSpwCoh(~withinContext);
        fs_within_cellpair{a,d} = fs_rotdist(withinContext);
        fs_across_cellpair{a,d} = fs_rotdist(~withinContext);

    end
end


fclose(fid);
