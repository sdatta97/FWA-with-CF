function obj = computeUtility(params,mean_rate_dl)
K_FWA = params.numCPE;
K = params.numCPE + params.numUE;
obj = params.price_FWA*sum(mean_rate_dl(1:K_FWA) > params.r_min_FWA) + params.price_cell*sum(mean_rate_dl((1+K_FWA):K) > params.r_min_cell);
end