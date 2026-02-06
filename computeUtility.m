function [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell,mean_rate_dl_FWA)
FWA_util = max(min(1,(mean_rate_dl_FWA - params.r_min_FWA)./(params.r_max_FWA - params.r_min_FWA)),0); 
cell_util = mean_rate_dl_cell >= params.r_min_cell;
end