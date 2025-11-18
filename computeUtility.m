function [cell_util, FWA_util] = computeUtility(params,mean_rate_dl)
K_FWA = params.numCPE;
K = params.numCPE + params.numUE;
% obj = params.price_FWA*sum(mean_rate_dl(1:K_FWA) > params.r_min_FWA) + params.price_cell*sum(mean_rate_dl((1+K_FWA):K) > params.r_min_cell);
% FWA_util = params.price_FWA*sum(mean_rate_dl(1:K_FWA) > params.r_min_FWA);
% cell_util = params.price_cell*sum(mean_rate_dl((1+K_FWA):K) > params.r_min_cell);
FWA_util = max(min(1,(exp(10^(-8)*mean_rate_dl(1:K_FWA)) - exp(10^(-8)*params.r_min_FWA))/(exp(10^(-8)*params.r_max_FWA) - exp(10^(-8)*params.r_min_FWA))),0);
cell_util = 1./(1+exp(-10^(-8)*(mean_rate_dl((1+K_FWA):K)-params.r_min_cell)));
end