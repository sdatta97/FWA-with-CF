function [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell,mean_rate_dl_FWA)
FWA_util = min(1,mean_rate_dl_FWA./params.r_max_FWA); %max(min(1,(exp(10^(-8)*mean_rate_dl_FWA) - exp(10^(-8)*params.r_min_FWA))/(exp(10^(-8)*params.r_max_FWA) - exp(10^(-8)*params.r_min_FWA))),0);
cell_util = mean_rate_dl_cell >= params.r_min_cell; %1./(1+exp(-10^(-8)*(mean_rate_dl_cell-params.r_min_cell)));
end