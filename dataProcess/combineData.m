close all;
clear;
%combineData: DATA GATHERING ONLY (plots live in plotData.m). Aggregates
%the per-seed (SLURM array) outputs of SimulationMain.m into stored
%files:
%  1) per-seed sweep tables (one row per sweep-grid point) in the
%     resultData/FWA_const_* folders (sweepNaming.m convention) are
%     combined into per-folder summary.csv/.txt with mean/std/median
%     across seeds;
%  2) raw packing matrices (per-user x per-(snapshot x realization)
%     rates) in resultData/FWA_packing_analysis_gIsite are reduced to the
%     spectrum-utilization curves f(eps) = q_eps(rate)/mean(rate) and
%     stored in packing_f_curves.csv.
%CSV-combining pattern adapted from
%https://in.mathworks.com/matlabcentral/answers/538119
repoRoot = fullfile(fileparts(mfilename('fullpath')),'..');
packDir = fullfile(repoRoot,'resultData','FWA_packing_analysis_gIsite');

%% 1) Combine the sweep CSVs across seeds and summarize.
%Result folders follow the sweepNaming.m convention: each
%resultData/FWA_const_<constant names>/ folder holds one schema (the
%swept parameters are its table columns), with one file per seed named
%by the constant values. Each folder is summarized independently; group
%keys are ALL non-outcome columns, so the aggregation adapts to
%whatever axes were swept.
outcomeVars = {'numCPE','numUE','init_FWA','Band_FWA','cell_se_ue','FWA_se'};
constDirs = dir(fullfile(repoRoot,'resultData','FWA_const_*'));
constDirs = constDirs([constDirs.isdir]);
for cd_i = 1:numel(constDirs)
    thisDir = fullfile(constDirs(cd_i).folder, constDirs(cd_i).name);
    dinfo = dir(fullfile(thisDir,'results_*.csv'));
    if isempty(dinfo), continue, end
    tables = cell(numel(dinfo),1);
    for f = 1:numel(dinfo)
        try
            tables{f} = readtable(fullfile(thisDir,dinfo(f).name));
        catch Exception
            disp(Exception);
        end
    end
    combinedTable = vertcat(tables{:});
    groupVars = setdiff(combinedTable.Properties.VariableNames, outcomeVars, 'stable');
    ovars = intersect(outcomeVars, combinedTable.Properties.VariableNames, 'stable');
    summaryTable = groupsummary(combinedTable, groupVars, {'mean','std','median'}, ovars);
    writetable(summaryTable, fullfile(thisDir,'summary.csv'));
    writetable(summaryTable, fullfile(thisDir,'summary.txt'));
    fprintf('sweep summary [%s]: %d rows from %d seed files -> summary.csv\n', ...
        constDirs(cd_i).name, height(summaryTable), numel(dinfo));
end
if isempty(constDirs)
    fprintf('no FWA_const_* result folders found under %s\n', fullfile(repoRoot,'resultData'));
end

%% 2) Tabulate the full-band per-terminal SE comparison across seeds
%(the gross physical-layer contrast written by SimulationMain at
%num_rep = 0: all UEs vs all CPEs, each on the whole band, no split, no
%feasibility filtering)
seDir = fullfile(repoRoot,'resultData','FWA_SE_comparison_gIsite');
seFiles = dir(fullfile(seDir,'se_comp_*.csv'));
if ~isempty(seFiles)
    seT = cell(numel(seFiles),1);
    for f = 1:numel(seFiles)
        seT{f} = readtable(fullfile(seDir,seFiles(f).name));
    end
    seAll = vertcat(seT{:});
    m = seAll.multiple;
    seSum = table(numel(m), mean(seAll.se_cell_ue), mean(seAll.se_fwa_cpe), ...
        mean(m), median(m), std(m)/sqrt(numel(m)), ...
        'VariableNames',{'n_seeds','mean_se_cell_ue','mean_se_fwa_cpe', ...
                         'mean_multiple','median_multiple','se_multiple'});
    writetable(seSum, fullfile(seDir,'se_comparison_summary.csv'));
    fprintf('full-band SE comparison (%d seeds): cell %.2f vs FWA %.2f b/s/Hz per terminal -> %.1fx (median %.1fx)\n', ...
        numel(m), seSum.mean_se_cell_ue, seSum.mean_se_fwa_cpe, seSum.mean_multiple, seSum.median_multiple);
end

%% 3) Reduce the packing matrices to spectrum-utilization curves
%MEAN-vs-PERCENTILE RATE AND SPECTRUM UTILIZATION. A user's delivered
%rate r fluctuates over the variability (fading, mobility). Let
%  r_bar   = E[r]        (mean deliverable rate), and
%  q_eps(r)              (the eps-quantile: rate exceeded w.p. 1 - eps).
%To guarantee a committed rate at outage eps the operator must plan for
%the WORST CASE, i.e. size the spectrum so the commitment is met even at
%q_eps(r). The link delivers r_bar on average, so the fraction of that
%provisioned capacity actually committed (monetized) is
%  f(eps) = q_eps(r) / r_bar   in (0, 1],
%with f -> 1 only as the rate variance -> 0. Equivalently, to deliver a
%target rate R per user the operator provisions B(eps) = R/q_eps(SE),
%versus R/mean(SE) if rates were constant, so 1/f(eps) is the spectrum
%over-provisioning factor and (1 - f(eps)) is the idle-spectrum margin.
%f(eps) is thus the fraction of spectrum utilized under outage-eps
%planning; f_FWA/f_cell is FWA's spectrum-utilization (revenue-per-Hz)
%multiple. Below, per user we sum q_eps over users / sum r_bar over users.
%Planning conventions (the asymmetry IS the point of the metric):
%cellular provisions for a RANDOM population of mobile users, so rates
%pool across drops (seeds), mobility snapshots, and fading realizations
%before taking the quantile - mirroring busy-hour worst-case
%dimensioning practice, where expansion triggers at ~70% PRB utilization
%(Survey on 5G Network Expansion Methods, IUCC 2021) and busy-hour
%utilization stays below ~80-85% (telecomHall capacity guidelines; RF
%Wireless World). FWA links are measured at installation and service can
%be declined on poor links (FCC 24-27), so the usable fraction is
%computed per drop from the per-CPE quantiles, then averaged over drops.
%The ratio f_FWA/f_cell is the packing multiple of revenue potential.
fwaFiles = dir(fullfile(packDir,'packing_FWA_*CPE_*act_*rep_*.csv'));
if ~isempty(fwaFiles)
    eps_grid = (0.01:0.01:0.20)';
    %one consistent (CPE count, pool) configuration - the smallest tags -
    %reduced separately at EVERY activity operating point present
    tok = regexp({fwaFiles.name},'packing_FWA_(\d+)CPE_([\d.]+)act_(\d+)rep_','tokens','once');
    cpeTags = cellfun(@(x) str2double(x{1}), tok);
    actTags = cellfun(@(x) str2double(x{2}), tok);
    repTags = cellfun(@(x) str2double(x{3}), tok);
    ncpe = min(cpeTags);
    nrep = min(repTags(cpeTags == ncpe));
    acts = unique(actTags(cpeTags == ncpe & repTags == nrep));
    allRows = [];
    for ai = 1:numel(acts)
        act = acts(ai);
        tag = sprintf('%dCPE_%gact_%drep', ncpe, act, nrep);
        sel = fwaFiles(cpeTags == ncpe & actTags == act & repTags == nrep);
        pool_cell_rep = [];
        pool_cell_plain = [];
        f_FWA_per_seed = zeros(numel(sel),numel(eps_grid));
        for f = 1:numel(sel)
            seed_id = erase(erase(sel(f).name,sprintf('packing_FWA_%s_',tag)),'.csv');
            rate_FWA = readmatrix(fullfile(packDir,sel(f).name));
            for ie = 1:numel(eps_grid)
                f_FWA_per_seed(f,ie) = sum(quantile(rate_FWA,eps_grid(ie),2))/sum(mean(rate_FWA,2));
            end
            fcr = fullfile(packDir,sprintf('packing_cell_rep_%s_%s.csv',tag,seed_id));
            if isfile(fcr)
                m = readmatrix(fcr); pool_cell_rep = [pool_cell_rep; m(:)]; %#ok<AGROW>
            end
            %plain-cellular baseline: an explicit CELL_REPEAT = 0 file if
            %present, else the num_rep = 0 cellular matrix (statistically
            %equivalent: CPEs present but inert, verified to 0.07%)
            fcp = fullfile(packDir,sprintf('packing_cell_plain_%s_%s.csv',tag,seed_id));
            if ~isfile(fcp)
                fcp = fullfile(packDir,sprintf('packing_cell_rep_%dCPE_%gact_0rep_%s.csv',ncpe,act,seed_id));
            end
            if isfile(fcp)
                m = readmatrix(fcp); pool_cell_plain = [pool_cell_plain; m(:)]; %#ok<AGROW>
            end
        end
        f_FWA = mean(f_FWA_per_seed,1)';
        f_cell_rep = nan(numel(eps_grid),1);
        f_cell_plain = nan(numel(eps_grid),1);
        if ~isempty(pool_cell_rep)
            f_cell_rep = arrayfun(@(e) quantile(pool_cell_rep,e)/mean(pool_cell_rep), eps_grid);
        end
        if ~isempty(pool_cell_plain)
            f_cell_plain = arrayfun(@(e) quantile(pool_cell_plain,e)/mean(pool_cell_plain), eps_grid);
        end
        allRows = [allRows; table(eps_grid, act*ones(size(eps_grid)), f_FWA, f_cell_rep, f_cell_plain, ...
            'VariableNames',{'eps','activity','f_FWA','f_cell_rep','f_cell_plain'})]; %#ok<AGROW>
        [~,i5] = min(abs(eps_grid - 0.05));
        if ~isnan(f_cell_rep(i5))
            fprintf('packing act %g (%d seeds, %s): eps 5%%: f_FWA %.3f, f_cell(NCR) %.3f, multiple %.2fx\n', ...
                act, numel(sel), tag, f_FWA(i5), f_cell_rep(i5), f_FWA(i5)/f_cell_rep(i5));
        end
    end
    writetable(allRows, fullfile(packDir,'packing_f_curves.csv'));
    fprintf('packing curves -> %s\n', fullfile(packDir,'packing_f_curves.csv'));
else
    fprintf('no packing matrices found in %s\n', packDir);
end
