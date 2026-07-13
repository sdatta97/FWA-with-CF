close all;
clear;
%combineData: DATA GATHERING ONLY (plots live in plotData.m). Aggregates
%the per-seed (SLURM array) outputs of SimulationMain.m into stored
%files:
%  1) sweep CSVs (one row per parameter combination per seed) in
%     resultData/FWA_multi_cell_repeater_fix_comp_alloc are combined into
%     summary.csv/.txt with mean/std/median across seeds;
%  2) raw packing matrices (per-user x per-(snapshot x realization)
%     rates) in resultData/FWA_packing_analysis are reduced to the
%     safe-load-fraction curves f(eps) = q_eps(rate)/mean(rate) and
%     stored in packing_f_curves.csv.
%CSV-combining pattern adapted from
%https://in.mathworks.com/matlabcentral/answers/538119
repoRoot = fullfile(fileparts(mfilename('fullpath')),'..');
sweepDir = fullfile(repoRoot,'resultData','FWA_multi_cell_repeater_fix_comp_alloc');
packDir = fullfile(repoRoot,'resultData','FWA_packing_analysis');

%% 1) Combine the sweep CSVs across seeds and summarize
dinfo = dir(fullfile(sweepDir,'*.csv'));
dinfo = dinfo(~startsWith({dinfo.name},'summary'));
%only the current output schema (older runs used different columns and
%cannot be concatenated); the activity tag marks the schema that carries
%the activity, demand_profile and rep_frac columns
dinfo = dinfo(contains({dinfo.name},'activity_'));
if ~isempty(dinfo)
    filenames = fullfile({dinfo.folder},{dinfo.name});
    tables = cell(numel(filenames),1);
    for f = 1:numel(filenames)
        try
            tables{f} = readtable(filenames{f});
        catch Exception
            disp(Exception);
        end
    end
    combinedTable = vertcat(tables{:});
    %group by the configuration columns; summarize the outcome columns
    groupVars = {'numCPE','activity','deployRange','r_min_cell', ...
        'demand_profile','num_rep','rep_frac','Band'};
    outcomeVars = {'init_FWA','max_FWA','Band_FWA','cell_se','FWA_se'};
    groupVars = intersect(groupVars, combinedTable.Properties.VariableNames,'stable');
    outcomeVars = intersect(outcomeVars, combinedTable.Properties.VariableNames,'stable');
    summaryTable = groupsummary(combinedTable,groupVars,{'mean','std','median'},outcomeVars);
    writetable(summaryTable, fullfile(sweepDir,'summary.csv'));
    writetable(summaryTable, fullfile(sweepDir,'summary.txt'));
    fprintf('sweep summary: %d rows from %d files -> %s\n', ...
        height(summaryTable), numel(filenames), fullfile(sweepDir,'summary.csv'));
else
    fprintf('no sweep CSVs found in %s\n', sweepDir);
end

%% 2) Reduce the packing matrices to safe-load-fraction curves
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
fwaFiles = dir(fullfile(packDir,'packing_FWA_*CPE_*rep_*.csv'));
if ~isempty(fwaFiles)
    eps_grid = (0.01:0.01:0.20)';
    %use the smallest CPE-count and num_rep tags present, so the curves
    %come from one consistent configuration
    tok = regexp({fwaFiles.name},'packing_FWA_(\d+)CPE_(\d+)rep_','tokens','once');
    cpeTags = cellfun(@(c) str2double(c{1}), tok);
    repTags = cellfun(@(c) str2double(c{2}), tok);
    ncpe = min(cpeTags);
    nrep = min(repTags(cpeTags == ncpe));
    tag = sprintf('%dCPE_%drep', ncpe, nrep);
    fwaFiles = fwaFiles(cpeTags == ncpe & repTags == nrep);
    pool_cell_rep = [];
    pool_cell_plain = [];
    f_FWA_per_seed = zeros(numel(fwaFiles),numel(eps_grid));
    for f = 1:numel(fwaFiles)
        seed_id = erase(erase(fwaFiles(f).name,sprintf('packing_FWA_%s_',tag)),'.csv');
        rate_FWA = readmatrix(fullfile(packDir,fwaFiles(f).name));
        for ie = 1:numel(eps_grid)
            f_FWA_per_seed(f,ie) = sum(quantile(rate_FWA,eps_grid(ie),2))/sum(mean(rate_FWA,2));
        end
        fcr = fullfile(packDir,sprintf('packing_cell_rep_%s_%s.csv',tag,seed_id));
        if isfile(fcr)
            m = readmatrix(fcr); pool_cell_rep = [pool_cell_rep; m(:)]; %#ok<AGROW>
        end
        fcp = fullfile(packDir,sprintf('packing_cell_plain_%s_%s.csv',tag,seed_id));
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
    packingTable = table(eps_grid, f_FWA, f_cell_rep, f_cell_plain, ...
        'VariableNames',{'eps','f_FWA','f_cell_rep','f_cell_plain'});
    writetable(packingTable, fullfile(packDir,'packing_f_curves.csv'));
    [~,i5] = min(abs(eps_grid - 0.05));
    fprintf('packing curves (%d seeds, %s pool) -> %s\n', numel(fwaFiles), tag, ...
        fullfile(packDir,'packing_f_curves.csv'));
    if ~isnan(f_cell_rep(i5))
        fprintf('eps = 5%%: f_FWA = %.3f, f_cell(NCR) = %.3f, packing multiple = %.2fx\n', ...
            f_FWA(i5), f_cell_rep(i5), f_FWA(i5)/f_cell_rep(i5));
    end
else
    fprintf('no packing matrices found in %s\n', packDir);
end
