function [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell,mean_rate_dl_FWA)
%FWA utility per CPE against its zone's ABSOLUTE demand rate for the
%active demand configuration (params.fwa_demand_rates, from the LOW/HIGH
%plan-anchored configs in the main script); falls back to the scalar
rmin = params.fwa_demand_rates(:); %per-CPE zone demands (mandatory)
FWA_util = max(min(1,(mean_rate_dl_FWA - rmin + 1)./(params.r_max_FWA - rmin)),0); 
cell_util = mean_rate_dl_cell >= params.r_min_cell;
end
