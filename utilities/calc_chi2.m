%% Chi2stat and pval sub-function
function [pval, chi2stat_out] = calc_chi2(metric_use, bins_use)
n = histcounts(metric_use, bins_use);
e_unif = sum(n)/length(n);  % Expected count in each bin for uniform distribution
chi2stat_out = sum((n-e_unif).^2/e_unif);
pval = chi2pdf(chi2stat_out,length(n)-1);
end
