close all;
clear;
%plot_packing: aggregates the per-seed outputs of PackingAnalysis.m and
%plots the safe load fraction f(eps) = q_eps(rate)/mean(rate) versus the
%planning outage target eps, for cellular UEs and FWA CPEs.
%
%Planning conventions (the asymmetry IS the point of the figure):
% - Cellular: the carrier provisions for a RANDOM population of mobile
%   users, so rates are pooled across deployment drops (seeds) AND
%   channel realizations before taking the quantile. This mirrors
%   busy-hour worst-case dimensioning practice, where expansion triggers
%   at ~70% PRB utilization (Survey on 5G Network Expansion Methods,
%   IUCC 2021) and busy-hour utilization is kept below ~80-85%
%   (telecomHall capacity guidelines; RF Wireless World).
% - FWA: each CPE's link is measured at installation and service can be
%   declined on poor links (FCC 24-27), so the usable fraction is
%   computed per drop from the per-CPE temporal quantiles, then averaged
%   across drops.
%The ratio f_FWA/f_cell is the packing multiple of revenue potential.
projectdir = fullfile(fileparts(mfilename('fullpath')),'..','resultData','FWA_packing_analysis');
eps_grid = 0.01:0.01:0.20;

pool_cell_rep = [];
pool_cell_plain = [];
f_FWA_per_seed = []; %seeds x numel(eps_grid)
files = dir(fullfile(projectdir,'packing_FWA_*.csv'));
assert(~isempty(files),'no PackingAnalysis outputs found in %s',projectdir);
for ff = 1:numel(files)
    seed_id = erase(erase(files(ff).name,'packing_FWA_'),'.csv');
    rate_FWA = readmatrix(fullfile(projectdir,files(ff).name));
    rate_cr = readmatrix(fullfile(projectdir,strcat('packing_cell_rep_',seed_id,'.csv')));
    rate_cp = readmatrix(fullfile(projectdir,strcat('packing_cell_plain_',seed_id,'.csv')));
    pool_cell_rep = [pool_cell_rep; rate_cr(:)]; %#ok<AGROW>
    pool_cell_plain = [pool_cell_plain; rate_cp(:)]; %#ok<AGROW>
    f_seed = zeros(1,numel(eps_grid));
    for ie = 1:numel(eps_grid)
        f_seed(ie) = sum(quantile(rate_FWA,eps_grid(ie),2))/sum(mean(rate_FWA,2));
    end
    f_FWA_per_seed = [f_FWA_per_seed; f_seed]; %#ok<AGROW>
end
fprintf('aggregated %d seeds\n', numel(files));

f_cell_rep = arrayfun(@(e) quantile(pool_cell_rep,e)/mean(pool_cell_rep), eps_grid);
f_cell_plain = arrayfun(@(e) quantile(pool_cell_plain,e)/mean(pool_cell_plain), eps_grid);
f_FWA = mean(f_FWA_per_seed,1);

figure; hold on; grid on;
plot(100*eps_grid, f_FWA, 'b-o', 'LineWidth', 1.5, 'MarkerIndices', 1:2:numel(eps_grid));
plot(100*eps_grid, f_cell_rep, 'r-s', 'LineWidth', 1.5, 'MarkerIndices', 1:2:numel(eps_grid));
plot(100*eps_grid, f_cell_plain, 'r--d', 'LineWidth', 1.5, 'MarkerIndices', 1:2:numel(eps_grid));
xlabel('Planning outage target, \epsilon (in %)');
ylabel('Safe load fraction of mean capacity, f(\epsilon)');
legend('FWA CPE (link known at install)', ...
       'Cellular UE, NCR-assisted (random population)', ...
       'Cellular UE, no NCRs (random population)', 'Location','southeast');
%annotate the packing multiple at the 5% planning target
[~,i5] = min(abs(eps_grid - 0.05));
mult5 = f_FWA(i5)/f_cell_rep(i5);
title(sprintf('Packing multiple at \\epsilon = 5%%: %.2fx', mult5));
fprintf('eps = 5%%: f_FWA = %.3f, f_cell(NCR) = %.3f, f_cell(plain) = %.3f, multiple = %.2fx\n', ...
    f_FWA(i5), f_cell_rep(i5), f_cell_plain(i5), mult5);

plotdir = fullfile(fileparts(mfilename('fullpath')),'..','plots');
savefig(fullfile(plotdir,'packing_analysis.fig'));
saveas(gcf, fullfile(plotdir,'packing_analysis.png'));
saveas(gcf, fullfile(plotdir,'packing_analysis.eps'), 'epsc');
