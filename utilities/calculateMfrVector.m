function [rateVector] = calculateMfrVector(s, label)

nTrials = length(s.trial);
rateVector = nan(1, nTrials);

for tr = 1:nTrials
    idx = strcmp(s.trial(tr).label, label);
    if ~any(idx)
        rateVector(tr) = nan;
    else
        rateVector(tr) = s.trial(tr).mfr{idx};
    end
end