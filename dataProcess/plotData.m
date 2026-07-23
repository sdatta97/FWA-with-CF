close all;
clear;
set(0,'DefaultFigureVisible','off');
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
takeRef = max(T.take_rate);       %largest take rate = most CPEs
%FIXED busy-hour device activity: 5%% - the most defensible US-suburban
%value (top of the 3-5%% suburban range, and it reproduces the ITU-R
%M.2412 ~10-users-per-TRxP evaluation density; see SimulationMain.m).
actRef = 0.05;
ua = unique(T.activity); if ~any(ua==actRef), [~,ia]=min(abs(ua-0.05)); actRef=ua(ia); end
tiers = {'low','medium','high'}; tierLbl = {'Low demand','Medium demand','High demand'};
tcol = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19]; tmk = {'-s','-d','-^'};

%% Dual-axis figure at the fixed activity: bandwidth freed for FWA (left,
%tier-invariant) and FWA subscribers served (right, one curve per demand
%tier), vs the number of NCRs. Error bars = +/-1 s.e. over seeds.
fig=figure('Visible','off'); hold on; grid on;
yyaxis left
rB = T.take_rate==takeRef & T.activity==actRef & strcmp(T.demand_profile,'low'); %tier-invariant
[xb,ob]=sort(T.num_rep(rB)); yb=T.mean_Band_FWA(rB)/1e6; yb=yb(ob);
eb=(T.std_Band_FWA(rB)/1e6)./sqrt(T.GroupCount(rB)); eb=eb(ob); eb(isnan(eb))=0;
hB=errorbar(xb,yb,eb,'-o','LineWidth',1.4,'Color','k','MarkerFaceColor','k','MarkerSize',4);
ylabel('Bandwidth freed for FWA (MHz)'); set(gca,'YColor','k');
yyaxis right
%only the high-demand tier is informative here: low/medium are already
%served in full (flat), so NCRs add subscribers only at the high tier
r = T.take_rate==takeRef & T.activity==actRef & strcmp(T.demand_profile,'high');
[x,o]=sort(T.num_rep(r)); y=T.mean_init_FWA(r); y=y(o);
e=(T.std_init_FWA(r))./sqrt(T.GroupCount(r)); e=e(o); e(isnan(e))=0;
hS=errorbar(x,y,e,'-s','LineWidth',1,'Color',tcol(1,:), ...
    'MarkerFaceColor',tcol(1,:),'MarkerSize',4);
ylabel('FWA subscribers served'); set(gca,'YColor',tcol(1,:));
xlabel('Number of NCRs, $N_{\mathrm{rep}}$','Interpreter','latex');
legend([hB;hS],{'Bandwidth freed','Subscribers served'},'Location','east','FontSize',7);
styleIEEE(fig); saveIEEE(fig,plotDir,'Band_and_K_vs_num_NCR');

%% Fig 4: packing safe-load fraction vs outage target (tier-invariant), 5% activity
packingFile = fullfile(packDir,'packing_f_curves.csv');
if isfile(packingFile)
    P = readtable(packingFile);
    if ismember('activity', P.Properties.VariableNames)
        uap=unique(P.activity); [~,iap]=min(abs(uap-0.05)); P=P(P.activity==uap(iap),:);
    end
    P = P(100*P.eps <= 10,:);   %truncate the outage axis to 1-10%
    fig=figure('Visible','off'); hold on; grid on; L={};
    plot(100*P.eps,P.f_FWA,'-o','LineWidth',1,'MarkerIndices',1:height(P)); L{end+1}='Fixed Wireless Access';
    if any(~isnan(P.f_cell_rep))
        plot(100*P.eps,P.f_cell_rep,'-s','LineWidth',1,'MarkerIndices',1:height(P)); L{end+1}='Mobile Cellular';
    end
    xlabel('Planning outage target, $\epsilon$ (\%)','Interpreter','latex');
    ylabel('Fraction of spectrum utilized, $f(\epsilon)$','Interpreter','latex');
    xlim([1 10]); xticks(1:1:10);
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
drawnow;
savefig(fig, fullfile(plotDir,[name '.fig']));
saveas(fig, fullfile(plotDir,[name '.png']));
saveas(fig, fullfile(plotDir,[name '.eps']), 'epsc');
end
