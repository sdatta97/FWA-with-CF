close all;
clear;
set(0,'DefaultFigureVisible','off');
%plotData: PLOTTING ONLY. Reads the files stored by combineData.m
%(summary.csv for the sweep, packing_f_curves.csv for the packing
%analysis) and produces the journal figures, styled for the IEEE
%template: single-column width (3.5 in), 8 pt Times, saved as
%.eps/.png/.fig in plots/. Run combineData.m first.
%Figures are drawn at the model's DEFAULT operating point (20% take,
%5% activity, low demand tier); sensitivity campaigns that sweep those
%axes are filtered back to the defaults via refMask. The available FWA
%bandwidth is tier-invariant (set by the cellular phase before the FWA demand
%loop); packing is tier-invariant and drawn once.
repoRoot = fullfile(fileparts(mfilename('fullpath')),'..');
%newest sweepNaming.m result folder whose summary carries the plot axes
cdirs = dir(fullfile(repoRoot,'resultData','FWA_const_*'));
cdirs = cdirs([cdirs.isdir]);
[~,ord] = sort([cdirs.datenum],'descend');
sweepDir = '';
%default-config summaries carry only num_rep as a swept axis (take rate,
%activity, and demand tier are constants baked into the folder name), so
%require only the NCR axis plus the plotted outcomes
for ci = ord
    cand = fullfile(cdirs(ci).folder, cdirs(ci).name);
    if isfile(fullfile(cand,'summary.csv'))
        vars = readtable(fullfile(cand,'summary.csv')).Properties.VariableNames;
        if all(ismember({'num_rep','mean_Band_FWA','mean_init_FWA'}, vars))
            sweepDir = cand; break
        end
    end
end
if isempty(sweepDir)
    error('no FWA_const_* folder with a plottable summary.csv (run combineData first)');
end
packDir = fullfile(repoRoot,'resultData','FWA_packing_analysis_gIsite');
plotDir = fullfile(repoRoot,'plots');

T = readtable(fullfile(sweepDir,'summary.csv'));
%REFERENCE operating point = the model defaults (20%% take, 5%% activity,
%low tier - see SimulationMain.m). Any of these axes may be a CONSTANT of
%the campaign (column absent from the summary); refMask filters only on
%the columns that exist.
actRef = 0.05;
if ismember('activity',T.Properties.VariableNames)
    ua = unique(T.activity); [~,ia]=min(abs(ua-actRef)); actRef=ua(ia);
end
refMask = true(height(T),1);
if ismember('take_rate',T.Properties.VariableNames), refMask = refMask & T.take_rate==max(T.take_rate); end
if ismember('activity',T.Properties.VariableNames),  refMask = refMask & T.activity==actRef; end
cSubs = [0 0.45 0.74]; %subscriber-series color

%% Dual-axis figure at the fixed activity: TOTAL bandwidth available for
%FWA use (left, tier-invariant) and FWA subscribers served (right), vs
%the number of NCRs. Error bars = +/-1 s.e. over seeds. The N = 0 value
%is the fallow band the cellular phase leaves without any NCR; the rise
%above it is the NCR-attributable relief.
fig=figure('Visible','off'); hold on; grid on;
yyaxis left
rB = refMask; %tier-invariant left axis
if ismember('demand_profile',T.Properties.VariableNames)
    rB = rB & strcmp(T.demand_profile,T.demand_profile{find(refMask,1)});
end
[xb,ob]=sort(T.num_rep(rB)); yb=T.mean_Band_FWA(rB)/1e6; yb=yb(ob);
eb=(T.std_Band_FWA(rB)/1e6)./sqrt(T.GroupCount(rB)); eb=eb(ob); eb(isnan(eb))=0;
hB=errorbar(xb,yb,eb,'-o','LineWidth',1.4,'Color','k','MarkerFaceColor','k','MarkerSize',4);
ylabel('Bandwidth available for FWA use (MHz)'); set(gca,'YColor','k');
yyaxis right
%subscribers served at the LOW tier - the model default, where service is
%band-limited and subscriber growth tracks the freed bandwidth ~1:1
r = refMask;
if ismember('demand_profile',T.Properties.VariableNames)
    r = r & strcmp(T.demand_profile,'low');
end
[x,o]=sort(T.num_rep(r)); y=T.mean_init_FWA(r); y=y(o);
e=(T.std_init_FWA(r))./sqrt(T.GroupCount(r)); e=e(o); e(isnan(e))=0;
hS=errorbar(x,y,e,'-s','LineWidth',1,'Color',cSubs, ...
    'MarkerFaceColor',cSubs,'MarkerSize',4);
ylabel('FWA subscribers served'); set(gca,'YColor',cSubs);
xlabel('Number of NCRs, $N_{\mathrm{rep}}$','Interpreter','latex');
legend([hB;hS],{'Bandwidth available','Subscribers served'},'Location','east','FontSize',7);
styleIEEE(fig); saveIEEE(fig,plotDir,'Band_and_K_vs_num_NCR');

%% Fig 4: spectrum utilization vs outage target (tier-invariant), at actRef
packingFile = fullfile(packDir,'packing_f_curves.csv');
if isfile(packingFile)
    P = readtable(packingFile);
    if ismember('activity', P.Properties.VariableNames)
        uap=unique(P.activity); [~,iap]=min(abs(uap-actRef)); P=P(P.activity==uap(iap),:);
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
