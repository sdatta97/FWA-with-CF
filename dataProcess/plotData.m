close all;
clear;
%plotData: PLOTTING ONLY. Reads the files stored by combineData.m
%(summary.csv for the sweep, packing_f_curves.csv for the packing
%analysis) and produces the journal figures, styled for the IEEE
%template: single-column width (3.5 in), 8 pt Times, saved as
%.eps/.png/.fig in plots/. Run combineData.m first.
%Figures are drawn at a REFERENCE operating point (the largest take
%rate, the activity closest to the 5%% busy-hour reference, and a
%mid-size NCR pool). Metrics that depend on the FWA demand tier
%(subscribers served, FWA SE) get one curve per tier; the freed
%bandwidth is tier-INVARIANT (set by the cellular phase before the FWA
%demand loop), so it is split by activity instead; packing is
%tier-invariant and drawn once.
repoRoot = fullfile(fileparts(mfilename('fullpath')),'..');
%newest sweepNaming.m result folder whose summary carries the plot axes
cdirs = dir(fullfile(repoRoot,'resultData','FWA_const_*'));
cdirs = cdirs([cdirs.isdir]);
[~,ord] = sort([cdirs.datenum],'descend');
sweepDir = '';
for ci = ord
    cand = fullfile(cdirs(ci).folder, cdirs(ci).name);
    if isfile(fullfile(cand,'summary.csv'))
        vars = readtable(fullfile(cand,'summary.csv')).Properties.VariableNames;
        if all(ismember({'num_rep','demand_profile','activity'}, vars))
            sweepDir = cand; break
        end
    end
end
if isempty(sweepDir)
    error('no FWA_const_* folder with a plottable summary.csv (run combineData first)');
end
packDir = fullfile(repoRoot,'resultData','FWA_packing_analysis');
plotDir = fullfile(repoRoot,'plots');

T = readtable(fullfile(sweepDir,'summary.csv'));
takeRef = max(T.take_rate);                       %largest take rate = most CPEs
[~,ia] = min(abs(unique(T.activity)-0.05)); ua=unique(T.activity); actRef=ua(ia);
reps = unique(T.num_rep); repRef = reps(min(1+find(reps>=6,1,'first')-1, numel(reps)));
if isempty(repRef), repRef = reps(round(end/2)); end
tiers = {'low','medium','high'}; tierLbl = {'Low demand','Medium demand','High demand'};
tierMark = {'-o','-s','-d'};

%% Fig 1: bandwidth freed for FWA vs #NCRs, one curve per activity
%(Band_FWA is tier-invariant, so activity is the informative legend axis)
fig=figure; hold on; grid on; L={};
for a = unique(T.activity)'
    rows = T.take_rate==takeRef & T.activity==a & strcmp(T.demand_profile,'low');
    if ~any(rows), continue, end
    [x,o]=sort(T.num_rep(rows)); y=T.mean_Band_FWA(rows)/1e6; y=y(o);
    e=T.std_Band_FWA(rows)/1e6; e=e(o); e(isnan(e))=0;
    errorbar(x,y,e,'-o','LineWidth',1); L{end+1}=sprintf('%.1f%% activity',100*a); %#ok<SAGROW>
end
xlabel('Number of NCRs, $N_{\mathrm{rep}}$','Interpreter','latex');
ylabel('Bandwidth freed for FWA (MHz)');
legend(L,'Location','northeast'); styleIEEE(fig); saveIEEE(fig,plotDir,'Band_FWA_vs_num_NCR');

%% Fig 2: FWA subscribers served vs #NCRs, one curve per demand tier
fig=figure; hold on; grid on;
for ti=1:numel(tiers)
    rows = T.take_rate==takeRef & T.activity==actRef & strcmp(T.demand_profile,tiers{ti});
    [x,o]=sort(T.num_rep(rows)); y=T.mean_init_FWA(rows); y=y(o);
    e=T.std_init_FWA(rows); e=e(o); e(isnan(e))=0;
    errorbar(x,y,e,tierMark{ti},'LineWidth',1);
end
xlabel('Number of NCRs, $N_{\mathrm{rep}}$','Interpreter','latex');
ylabel('FWA subscribers served');
legend(tierLbl,'Location','southeast'); styleIEEE(fig); saveIEEE(fig,plotDir,'K_FWA_vs_num_NCR');

%% Fig 3: FWA subscribers served vs busy-hour activity, one curve per tier
fig=figure; hold on; grid on;
for ti=1:numel(tiers)
    rows = T.take_rate==takeRef & T.num_rep==repRef & strcmp(T.demand_profile,tiers{ti});
    [x,o]=sort(T.activity(rows)); y=T.mean_init_FWA(rows); y=y(o);
    e=T.std_init_FWA(rows); e=e(o); e(isnan(e))=0;
    errorbar(100*x,y,e,tierMark{ti},'LineWidth',1);
end
xlabel('Busy-hour device activity (\%)','Interpreter','latex');
ylabel('FWA subscribers served');
legend(tierLbl,'Location','southwest'); styleIEEE(fig); saveIEEE(fig,plotDir,'K_FWA_vs_activity');

%% Fig 4: packing safe-load fraction vs outage target (tier-invariant), 5% activity
packingFile = fullfile(packDir,'packing_f_curves.csv');
if isfile(packingFile)
    P = readtable(packingFile);
    if ismember('activity', P.Properties.VariableNames)
        uap=unique(P.activity); [~,iap]=min(abs(uap-0.05)); P=P(P.activity==uap(iap),:);
    end
    fig=figure; hold on; grid on; L={};
    plot(100*P.eps,P.f_FWA,'-o','LineWidth',1,'MarkerIndices',1:3:height(P)); L{end+1}='FWA CPE (link known at install)';
    if any(~isnan(P.f_cell_rep))
        plot(100*P.eps,P.f_cell_rep,'-s','LineWidth',1,'MarkerIndices',1:3:height(P)); L{end+1}='Cellular UE';
    end
    xlabel('Planning outage target, $\epsilon$ (\%)','Interpreter','latex');
    ylabel('Safe load fraction, $f(\epsilon)$','Interpreter','latex');
    legend(L,'Location','east'); styleIEEE(fig); saveIEEE(fig,plotDir,'packing_analysis');
else
    fprintf('%s not found - run combineData first\n', packingFile);
end

%% IEEE journal template styling and export
function styleIEEE(fig)
set(fig,'Units','inches','Position',[1 1 3.5 2.5],'Color','w');
ax=findall(fig,'Type','axes'); set(ax,'FontName','Times New Roman','FontSize',8,'Box','on','LineWidth',0.5);
set(findall(fig,'Type','legend'),'FontName','Times New Roman','FontSize',7);
end
function saveIEEE(fig, plotDir, name)
if ~isfolder(plotDir), mkdir(plotDir), end
savefig(fig, fullfile(plotDir,[name '.fig']));
saveas(fig, fullfile(plotDir,[name '.png']));
saveas(fig, fullfile(plotDir,[name '.eps']), 'epsc');
end
