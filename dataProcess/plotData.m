close all;
clear;
%plotData: PLOTTING ONLY. Reads the files stored by combineData.m
%(summary.csv for the sweep, packing_f_curves.csv for the packing
%analysis) and produces the journal figures, styled for the IEEE
%template: single-column width (3.5 in), 8 pt Times, saved as
%.eps/.png/.fig in plots/. Run combineData.m first.
repoRoot = fullfile(fileparts(mfilename('fullpath')),'..');
sweepDir = fullfile(repoRoot,'resultData','FWA_multi_cell_repeater_fix_comp_alloc');
packDir = fullfile(repoRoot,'resultData','FWA_packing_analysis');
plotDir = fullfile(repoRoot,'plots');

%% 1) Sweep figures
summaryFile = fullfile(sweepDir,'summary.csv');
if isfile(summaryFile)
    summaryTable = readtable(summaryFile);
    markers = {'-o','-s','-d','-^','-v','->'};
    %(a) FWA subscribers served vs the subscriber-plan rate tier, one
    %curve per number of enabled NCRs (error bars: std across seeds)
    fig = figure; hold on; grid on;
    reps = unique(summaryTable.num_rep);
    legends = cell(numel(reps),1);
    for i = 1:numel(reps)
        rows = summaryTable.num_rep == reps(i);
        [~,x] = ismember(lower(string(summaryTable.demand_profile(rows))), ["low" "medium" "high"]);
        [x,ord] = sort(x); %profile index: 1=low, 2=medium, 3=high
        y = summaryTable.mean_max_FWA(rows); y = y(ord);
        e = summaryTable.std_max_FWA(rows); e = e(ord);
        e(isnan(e)) = 0; %single-seed data has no spread yet
        errorbar(x, y, e, markers{1+mod(i-1,numel(markers))}, 'LineWidth', 1);
        legends{i} = sprintf('%d NCRs', reps(i));
    end
    xlabel('Minimum FWA rate, $r^{\mathrm{F}}_{\min}$ (in Mbps)','Interpreter','latex');
    ylabel('Number of FWA subscribers served');
    legend(legends,'Location','northeast');
    styleIEEE(fig);
    saveIEEE(fig, plotDir, 'K_FWA_vs_rmin');
    %(b) bandwidth freed for FWA vs the number of enabled NCRs, at the
    %entry-tier plan rate (Fig. 5 of the paper)
    fig = figure; hold on; grid on;
    rows0 = lower(string(summaryTable.demand_profile)) == "low"; %lowest profile
    [x,ord] = sort(summaryTable.num_rep(rows0));
    y = summaryTable.mean_Band_FWA(rows0)/1e6; y = y(ord);
    e = summaryTable.std_Band_FWA(rows0)/1e6; e = e(ord); e(isnan(e)) = 0;
    errorbar(x, y, e, '-o', 'LineWidth', 1);
    xlabel('Number of NCRs enabled, $N_{\mathrm{rep}}$','Interpreter','latex');
    ylabel('Bandwidth available for FWA (in MHz)');
    styleIEEE(fig);
    saveIEEE(fig, plotDir, 'Band_FWA_vs_num_NCR');
else
    fprintf('%s not found - run combineData.m first\n', summaryFile);
end

%% 2) Packing figure: safe load fraction vs planning outage target
packingFile = fullfile(packDir,'packing_f_curves.csv');
if isfile(packingFile)
    packingTable = readtable(packingFile);
    fig = figure; hold on; grid on;
    legends = {};
    plot(100*packingTable.eps, packingTable.f_FWA, '-o', 'LineWidth', 1, ...
        'MarkerIndices', 1:3:height(packingTable));
    legends{end+1} = 'FWA CPE (link known at install)';
    if any(~isnan(packingTable.f_cell_rep))
        plot(100*packingTable.eps, packingTable.f_cell_rep, '-s', 'LineWidth', 1, ...
            'MarkerIndices', 1:3:height(packingTable));
        legends{end+1} = 'Cellular UE, NCR-assisted';
    end
    if any(~isnan(packingTable.f_cell_plain))
        plot(100*packingTable.eps, packingTable.f_cell_plain, '--d', 'LineWidth', 1, ...
            'MarkerIndices', 1:3:height(packingTable));
        legends{end+1} = 'Cellular UE, no NCRs';
    end
    xlabel('Planning outage target, $\epsilon$ (in \%)','Interpreter','latex');
    ylabel('Safe load fraction of capacity, $f(\epsilon)$','Interpreter','latex');
    legend(legends,'Location','east');
    styleIEEE(fig);
    saveIEEE(fig, plotDir, 'packing_analysis');
else
    fprintf('%s not found - run combineData.m first\n', packingFile);
end

%% IEEE journal template styling and export
function styleIEEE(fig)
%single-column IEEE figure: 3.5 in wide, 8 pt Times throughout
set(fig,'Units','inches','Position',[1 1 3.5 2.4],'Color','w');
ax = findall(fig,'Type','axes');
set(ax,'FontName','Times New Roman','FontSize',8,'Box','on','LineWidth',0.5);
lg = findall(fig,'Type','legend');
set(lg,'FontName','Times New Roman','FontSize',7);
end

function saveIEEE(fig, plotDir, name)
if not(isfolder(plotDir))
    mkdir(plotDir)
end
savefig(fig, fullfile(plotDir,[name '.fig']));
saveas(fig, fullfile(plotDir,[name '.png']));
saveas(fig, fullfile(plotDir,[name '.eps']), 'epsc');
end
