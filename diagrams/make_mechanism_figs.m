%MECHANISM SCHEMATICS for the FWA-cellular paper: side-elevation vignettes
%reusing the deployment figure's glyph/color vocabulary. Sized ~1.3:1 for
%one column of a two-column template; no baked-in caption (use \caption).
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
fig = figure('Position',[100 100 430 325],'Color','w'); hold on
groundLine(0, 70);
drawMast(7, 32, col);                               % serving gNB
drawBuilding(43, 11, 30, col.apt, col);             % tall apartment = blocker
drawBuilding(24, 9, 15, col.office, col);           % office w/ NCR-CPE
drawCPE(24, 15, col);
for xt = [53 59 68], drawTree(xt, 8, col); end
drawUE(64, col);                                    % cell-edge UE in NLOS

% direct path (blocked) and relayed path gNB -> NCR -> UE
arc([7 32],[64 3],[25 39], col.blk, 1.1, ':');
plot(43, 30, 'x','Color',col.blk,'MarkerSize',10,'LineWidth',1.8);
arc([7 32],[24 16],[14 32], col.relay, 1.6, '-');
beamCone([24 16],[64 3], 2.6, col.relay);
arc([24 16],[64 3],[44 20], col.relay, 1.6, '-');

nodeLabel(7,-4,'gNB',[0 0 0]);
nodeLabel(24,21,'NCR',col.cpe);
nodeLabel(64,-4,'cellular UE',col.ue);
legendBox(4,47,{{':',col.blk,'weak direct path'},{'-',col.relay,'NCR-relayed path'}});
finishAxes([-3 71],[-7 49]);
drawnow; exportgraphics(fig, fullfile(outdir,'mech_ncr_support.png'),'Resolution',300);
drawnow; exportgraphics(fig, fullfile(outdir,'mech_ncr_support.eps'));
savefig(fig, fullfile(outdir,'mech_ncr_support.fig'));

%% ============ Figure B: interference nulling at the CPE ============
fig = figure('Position',[100 100 430 337],'Color','w'); hold on
groundLine(0, 70);
drawMast(7, 32, col);                               % serving gNB
drawMast(63, 32, col);                              % interfering gNB
drawBuilding(35, 12, 18, col.office, col);          % FWA CPE building
drawArray(35, 18);                                  % 8-element rooftop array
for xt = [17 50 55], drawTree(xt, 7, col); end

% CPE receive pattern (lobe to serving gNB, null to interferer), then the paths
rxPattern(35, 25, 10, col.cpe);
beamCone([7 32],[34 21], 2.2, col.beam);
arc([7 32],[34 21],[19 30], col.beam, 1.7, '-');
arc([63 32],[41 24],[52 30], col.intf, 1.3, '--');
plot(44, 24, 'x','Color',col.intf,'MarkerSize',11,'LineWidth',2);
text(47,15,'receive null','Color',col.intf,'FontSize',8,'FontName','Times New Roman');

nodeLabel(7,-4,'serving gNB',[0 0 0]);
nodeLabel(63,-4,'interfering gNB',col.intf);
nodeLabel(35,-4,'CPE',col.cpe);
legendBox(16,49,{{'-',col.beam,'desired FWA beam'},{'--',col.intf,'interference (nulled)'}});
finishAxes([-3 71],[-7 51]);
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
    hd=b(end,:); s=2.0;
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
function nodeLabel(x,y,txt,c)
    text(x,y,txt,'Color',c,'FontSize',8,'FontWeight','bold', ...
        'FontName','Times New Roman','HorizontalAlignment','center');
end
function legendBox(x0,y0,entries)
    dy=4.4; lw=7;
    for i=1:numel(entries)
        yy=y0-(i-1)*dy; e=entries{i};
        plot([x0 x0+lw],[yy yy], e{1},'Color',e{2},'LineWidth',1.6);
        text(x0+lw+2, yy, e{3},'Color',[0.15 0.15 0.15],'FontSize',8, ...
            'FontName','Times New Roman','VerticalAlignment','middle');
    end
end
function finishAxes(xl,yl)
    axis equal; axis off; xlim(xl); ylim(yl);
    set(gcf,'InvertHardcopy','off');
end
