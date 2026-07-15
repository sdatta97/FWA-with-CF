%MECHANISM SCHEMATICS for the FWA-cellular paper: side-elevation vignettes
%reusing the deployment figure's glyph/color vocabulary.
% (A) NCR-aided support: relayed cellular path over a blocker.
% (B) Interference nulling: static CPE array nulls the interfering gNB.
outdir = '/Users/sdatta/FWA-with-CF/plots';
if ~isfolder(outdir), mkdir(outdir), end
close all; set(0,'DefaultFigureVisible','off');

col.mast   = [0.15 0.15 0.18];
col.apt    = [0.35 0.55 0.80];
col.office = [0.90 0.55 0.30];
col.home   = [0.85 0.80 0.62];
col.cpe    = [0.05 0.15 0.50];
col.tree   = [0.23 0.49 0.25];
col.ue     = [0.72 0.10 0.16];
col.beam   = [0.00 0.45 0.74];
col.relay  = [0.20 0.55 0.25];
col.intf   = [0.72 0.10 0.16];
col.blk    = [0.55 0.55 0.55];

%% ============ Figure A: NCR-aided cellular support ============
fig = figure('Position',[100 100 520 300],'Color','w'); hold on
groundLine(0, 100);
drawMast(8, 33, col);                              % serving gNB
drawBuilding(50, 11, 30, col.apt, col);            % tall apartment = blocker
drawBuilding(29, 9, 15, col.office, col);          % office w/ NCR-CPE
drawCPE(29, 15, col);
for xt = [64 71 78 88 94], drawTree(xt, 7+2*mod(xt,3), col); end
drawUE(80, col);                                   % cell-edge UE in NLOS

% direct path (blocked by the apartment)
arc([8 33],[80 3],[30 44], col.blk, 1.1, ':');
plot(50, 30, 'x', 'Color',col.blk,'MarkerSize',11,'LineWidth',2);
text(67,22,'weak direct path (NLOS)','Color',col.blk,'FontSize',8,'FontName','Times New Roman','HorizontalAlignment','center');

% relayed path: gNB -> NCR, then NCR -> UE (with a service beam cone)
arc([8 33],[29 16],[18 33], col.relay, 1.6, '-');
beamCone([29 16],[80 3], 3, col.relay);
arc([29 16],[80 3],[55 20], col.relay, 1.6, '-');
text(17,29,'gNB\rightarrowNCR','Color',col.relay,'FontSize',8,'FontName','Times New Roman');
text(52,15,'NCR\rightarrowUE (cellular band)','Color',col.relay,'FontSize',8,'FontName','Times New Roman');
text(29,20,'NCR','Color',col.cpe,'FontSize',8,'FontWeight','bold','FontName','Times New Roman','HorizontalAlignment','center');

text(50,-9,'NCR boosts cell-edge SE \rightarrow frees bandwidth for FWA', ...
    'FontSize',8.5,'FontAngle','italic','FontName','Times New Roman','HorizontalAlignment','center');
finishAxes([-2 100],[-13 50]);
drawnow; exportgraphics(fig, fullfile(outdir,'mech_ncr_support.png'),'Resolution',300);
drawnow; exportgraphics(fig, fullfile(outdir,'mech_ncr_support.eps'));
savefig(fig, fullfile(outdir,'mech_ncr_support.fig'));

%% ============ Figure B: interference nulling at the CPE ============
fig = figure('Position',[100 100 520 300],'Color','w'); hold on
groundLine(0, 100);
drawMast(8, 33, col);                              % serving gNB
drawMast(92, 33, col);                             % interfering gNB
drawBuilding(47, 12, 18, col.office, col);         % FWA CPE building
drawArray(47, 18);                                 % 8-element rooftop array
for xt = [24 70 78], drawTree(xt, 7, col); end

% CPE receive pattern first (so paths overlay it): lobe to serving gNB, null to interferer
rxPattern(47, 25, 11, col.cpe);
text(40,10,'CPE','Color',col.cpe,'FontSize',8,'FontWeight','bold','FontName','Times New Roman','HorizontalAlignment','right');

% desired beam from serving gNB
beamCone([8 33],[45 22], 2.4, col.beam);
arc([8 33],[45 22],[26 31], col.beam, 1.7, '-');
text(19,32,'desired FWA beam','Color',col.beam,'FontSize',8,'FontName','Times New Roman');

% interference from the other gNB, nulled at the CPE
arc([92 33],[54 24],[72 31], col.intf, 1.3, '--');
plot(57, 24, 'x','Color',col.intf,'MarkerSize',12,'LineWidth',2.2);
text(74,31,'interfering gNB','Color',col.intf,'FontSize',8,'FontName','Times New Roman');
text(60,17,'receive null','Color',col.intf,'FontSize',8,'FontName','Times New Roman');

text(50,-9,'Static CPE + 8-antenna array \rightarrow interference null \rightarrow higher FWA SE', ...
    'FontSize',8.5,'FontAngle','italic','FontName','Times New Roman','HorizontalAlignment','center');
finishAxes([-2 102],[-13 52]);
drawnow; exportgraphics(fig, fullfile(outdir,'mech_interf_null.png'),'Resolution',300);
drawnow; exportgraphics(fig, fullfile(outdir,'mech_interf_null.eps'));
savefig(fig, fullfile(outdir,'mech_interf_null.fig'));
fprintf('mechanism figures saved to %s\n', outdir);

%% ---------- local drawing helpers ----------
function groundLine(x0,x1)
    plot([x0-2 x1+2],[0 0],'-','Color',[0.4 0.4 0.4],'LineWidth',1);
end
function drawMast(x,h,col)
    fill(x+[-0.6 0.6 0.6 -0.6], [0 0 h h], col.mast,'EdgeColor','none');
    fill(x+[-3 3 0],[h h h+5], [0 0 0],'EdgeColor','none'); % panel/triangle head
    fill(x+[-3.5 3.5 3.5 -3.5],[h-2 h-2 h h], [0.1 0.1 0.12],'EdgeColor','none');
end
function drawBuilding(x,w,h,c,col)
    fill(x+[-w/2 w/2 w/2 -w/2],[0 0 h h], c,'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.4);
    fill(x+[-w/2 w/2 w/2 -w/2],[h-1.5 h-1.5 h h],0.8*c,'EdgeColor','none'); % roof band
    for wy = 4:6:h-3   % window rows
        for wx = x-w/2+2 : 3.5 : x+w/2-2
            fill(wx+[-0.7 0.7 0.7 -0.7],[wy wy wy+2 wy+2],[0.95 0.97 1],'EdgeColor','none');
        end
    end
end
function drawCPE(x,z,col)
    fill(x+[-1.6 1.6 1.6 -1.6],[z z z+2.4 z+2.4], col.cpe,'EdgeColor','none');
    plot([x x],[z+2.4 z+4.2],'-','Color',col.cpe,'LineWidth',1); % small stub antenna
end
function drawArray(x,z)
    for e = -2.5:1:2.5
        plot(x+e+[0 0],[z z+3.2],'-','Color',[0.05 0.15 0.5],'LineWidth',1.4);
    end
end
function drawTree(x,h,col)
    fill(x+[-0.5 0.5 0.5 -0.5],[0 0 h*0.35 h*0.35],[0.4 0.28 0.15],'EdgeColor','none');
    fill(x+[-2.6 2.6 0],[h*0.3 h*0.3 h], col.tree,'EdgeColor','none');
    fill(x+[-2.1 2.1 0],[h*0.55 h*0.55 h*1.2], col.tree,'EdgeColor','none');
end
function drawUE(x,col)
    fill(x+[-2.5 2.5 2.5 -2.5],[0 0 2.6 2.6], col.ue,'EdgeColor','none');       % car body
    fill(x+[-1.4 1.4 1.0 -1.0],[2.6 2.6 4.4 4.4], col.ue,'EdgeColor','none');   % cabin
end
function arc(p0,p2,ctrl,c,lw,ls)
    t=linspace(0,1,60)';
    b=(1-t).^2.*p0 + 2*(1-t).*t.*ctrl + t.^2.*p2;
    plot(b(:,1),b(:,2),ls,'Color',c,'LineWidth',lw);
    d=b(end,:)-b(end-1,:); d=d/norm(d); n=[-d(2) d(1)];
    hd=b(end,:); s=2.2;
    fill([hd(1) hd(1)-s*d(1)+0.5*s*n(1) hd(1)-s*d(1)-0.5*s*n(1)], ...
         [hd(2) hd(2)-s*d(2)+0.5*s*n(2) hd(2)-s*d(2)-0.5*s*n(2)], c,'EdgeColor','none');
end
function beamCone(tx,rx,halfw,c)
    d=rx-tx; L=norm(d); u=d/L; n=[-u(2) u(1)];
    apex=tx+u*4;
    fill([apex(1) rx(1)+halfw*n(1) rx(1)-halfw*n(1)], ...
         [apex(2) rx(2)+halfw*n(2) rx(2)-halfw*n(2)], c,'EdgeColor','none','FaceAlpha',0.12);
end
function rxPattern(x,z,sc,c)
    th=linspace(0,2*pi,220);
    r=sc*(0.08+0.92*(0.5*(1-cos(th))).^1.7); % deep null toward +x (interferer), sharp lobe toward -x
    px=x+r.*cos(th); pz=z+0.62*r.*sin(th);
    fill(px,pz,c,'EdgeColor','none','FaceAlpha',0.18);
    plot(px,pz,'-','Color',c,'LineWidth',0.8);
end
function finishAxes(xl,yl)
    axis equal; axis off; xlim(xl); ylim(yl);
    set(gcf,'InvertHardcopy','off');
end
