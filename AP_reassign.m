% function D_after_handoff =  AP_reassign(D, chgains, K_mmW, k_idx, p_fac, ue_idxs)
function [D_after_handoff, ue_idxs] =  AP_reassign(params, ue_idx)
D = params.D;
chgains = params.BETA;
D_after_handoff = D;
loss_pc_thresh = params.loss_pc_thresh;
% M = size(chgains,1);
% K = size(chgains,2);
M = params.numGNB;
K_mmW = params.numUE;
K = params.numUE + params.numUE_sub6;
ap_idxs = find(D(:,ue_idx));
ue_idxs = [];
for i = 1:length(ap_idxs)
    ap = ap_idxs(i);
    ue_idxs = union(ue_idxs, find(D(ap,:)));
end
% for k = 1:(K-K_mmW)
ue_idxs = ue_idxs ((1+K_mmW):end);
loss_pc_arr = zeros(size(ue_idxs));
for kk_idx = 1:length(ue_idxs)
    k = ue_idxs(kk_idx);
    ap_idxs_k = find(D(:,k));
    ap_idxs_affected = intersect(ap_idxs,ap_idxs_k);
    loss_pc_arr(kk_idx) = 100*(sum(chgains(ap_idxs_affected,k))/sum(chgains(ap_idxs_k,k)));
end
ue_idxs = ue_idxs(loss_pc_arr > loss_pc_thresh);
for kk_idx = 1:length(ue_idxs)
    k = ue_idxs(kk_idx);
    ap_idxs_k = find(D(:,k));
    other_ap_idxs = setdiff(1:M,union(ap_idxs,ap_idxs_k));
    [~,other_ap_idxs_idxs] = sort(chgains(other_ap_idxs,k),'descend'); 
    other_ap_idxs = other_ap_idxs(other_ap_idxs_idxs);
    ap_idxs_affected = intersect(ap_idxs,ap_idxs_k);
    % for m = 1:length(ap_idxs_affected)
    for m = 1:min(length(ap_idxs_affected),length(other_ap_idxs))
        if (numel(other_ap_idxs) > 0)
%             other_ue_idxs = setdiff(find(D(other_ap_idxs(m),:)),k);
%             [~,k_idx_idx] = min(chgains(other_ap_idxs(m),other_ue_idxs)); 
%             k_idx = other_ue_idxs(k_idx_idx);
%             if (numel(k_idx)>0)
%                 D_after_handoff(other_ap_idxs(m),k) = 1;
%                 D_after_handoff(other_ap_idxs(m),k_idx) = 0;
%             end
            D_after_handoff(other_ap_idxs(m),k) = 1;
            D_after_handoff(ap_idxs_affected(m),k) = 0;
        end
    end
end
end