%COMBINED SYSTEM-MODEL FIGURE: one single scene - the two-cell icon
%deployment - with BOTH paper mechanisms drawn in place:
%  - NCR-aided cellular support (left cell): the serving gNB's weak direct
%    path to a cell-edge in-car UE is blocked (grey dotted, x), and the
%    NCR hosted on a cell-edge office-park rooftop relays it (green
%    radiation lobes);
%  - CPE interference nulling (right cell): a rooftop CPE receives the
%    desired beam of its serving gNB (blue transmit lobe + CPE receive
%    pattern) while the co-channel lobe from the OTHER site's gNB (red)
%    is nulled by the CPE array (dashed remnant, x).
%Beams are drawn as directional radiation lobes, not arrows. Small path
%legends sit in the sky corners; the icon legend is inherited.
%REUSES make_deployment_3d.m rather than duplicating it: sceneThin = 1
%asks that script for a sparser display (2 apartment blocks per cell,
%half the homes/trees/cars) and suppresses its standalone outputs; the
%overlays are then drawn onto its axes using the actual scene object
%positions left in the workspace.
%Output: plots/deployment_combined.{png,eps,fig}.
sceneThin = 1;
run(fullfile(fileparts(mfilename('fullpath')),'make_deployment_3d.m'));
%workspace now holds the scene: fig, axMap, D (drawables with display
%positions/heights), prj/kx/ky (projection), dispC/Dx (cell centres), ICO,
%plotDir, params, ...

%% mechanism colors (same vocabulary as make_mechanism_figs.m)
colRelay = [0.20 0.55 0.25];
colBlk   = [0.50 0.50 0.50];
colBeam  = [0.00 0.45 0.74];
colIntf  = [0.72 0.10 0.16];
colNcr   = [0.05 0.15 0.50];

keys  = {D.key}; xs = [D.x]; roofs = [D.roof];
topR = towerTop(D, dispC(1,1), prj);   %right-cell antenna head
topL = towerTop(D, dispC(2,1), prj);   %left-cell antenna head

%% --- NCR-aided support, LEFT cell ---
%NCR host: the outermost office-park building (per the model, NCRs are
%CPE-hosted at the cell edge near the noise-limited road users); prefer
%one drawn without a rooftop dish so the NCR unit stands alone.
io = find(strcmp(keys,'office') & xs < 0 & ~roofs);
if isempty(io), io = find(strcmp(keys,'office') & xs < 0); end
[~,ii] = max(abs(xs(io))); o = D(io(ii));
[ox,oy] = prj(o.x,o.y); roofY = oy + o.H + 3;
%assisted UE: an in-car user a readable distance from the NCR office
ic = find(strcmp(keys,'car') & xs < 0);
dw = hypot(xs(ic)-o.x, [D(ic).y]-o.y); dw(dw < 80) = inf; %not right under it
[~,ii] = min(dw); c = D(ic(ii));
[ux,uy] = prj(c.x,c.y); carPt = [ux, uy+22];
%the NCR unit on the office rooftop (label on the open, outward side)
patch(axMap, ox+[-14 14 14 -14], roofY+[0 0 11 11], colNcr, 'EdgeColor','none');
plot(axMap, [ox ox], roofY+[11 22], '-', 'Color',colNcr, 'LineWidth',1.2);
text(axMap, ox-22, roofY+16, 'NCR', 'Color',colNcr, 'FontName','Times New Roman', ...
    'FontSize',11, 'FontWeight','bold', 'HorizontalAlignment','right');
%blocked direct path (grey dotted, x) and the relayed radiation lobes
arcPath(axMap, topL, carPt, 150, colBlk, 1.3, ':');
bmid = bezPoint(topL, carPt, 150, 0.5);
plot(axMap, bmid(1), bmid(2), 'x', 'Color',colBlk, 'MarkerSize',13, 'LineWidth',2.4);
beamLobe(axMap, topL, [ox, roofY+16], 1.0, colRelay, 80);
beamLobe(axMap, [ox, roofY+16], carPt, 1.0, colRelay, 60);

%% --- CPE interference nulling, RIGHT cell ---
%victim CPE: a dish-bearing building on the join-facing half of the right
%cell, on the UPPER (far) side so both beams approach through open sky
ib = find(roofs & xs > 280 & xs < Dx - 100);
if isempty(ib), ib = find(roofs & xs > 100 & xs < Dx); end %thinned fallback
[~,ii] = min(xs(ib)); b = D(ib(ii)); %nearest the join: longest, clearest links
[bx,by] = prj(b.x,b.y); dishPt = [bx, by + b.H*1.04 + 19];
%CPE receive pattern: main lobe toward the serving gNB (right), deep null
%toward the interfering site (left) - the nulling mechanism itself
thp = linspace(0,2*pi,220);
rp = 66*(0.08 + 0.92*(0.5*(1+cos(thp))).^1.7); %null at theta = pi (toward -x)
patch(axMap, dishPt(1)+rp.*cos(thp), dishPt(2)+0.62*rp.*sin(thp), colBeam, ...
      'EdgeColor',colBeam, 'FaceAlpha',0.16, 'LineWidth',0.7);
%desired transmit lobe from the serving (right) gNB
beamLobe(axMap, topR, dishPt, 1.0, colBeam, 60);
%interfering lobe from the LEFT site's gNB: radiates toward the victim,
%continues as a dashed remnant and is nulled short of the CPE (x)
nulStop = dishPt + [-75 24];
dIL = norm(dishPt - topL); dirIL = (dishPt - topL)/dIL;
beamLobe(axMap, topL, dishPt, 0.22, colIntf, 80);
arcPath(axMap, topL + 0.20*dIL*dirIL, nulStop, 130, colIntf, 1.7, '--');
nulX = nulStop + 0.4*(dishPt - nulStop);
plot(axMap, nulX(1), nulX(2), 'x', 'Color',colIntf, 'MarkerSize',13, 'LineWidth',2.4);

%% --- path legends in the sky corners ---
xL = xlim(axMap); yL = ylim(axMap); rx = diff(xL); ry = diff(yL);
segL = 0.040*rx; dy = 0.085*ry;
lx = xL(1) + 0.012*rx; ly = yL(2) - 0.055*ry;      %top-left: NCR mechanism
miniLine(axMap, lx, ly,    segL, ':', colBlk, 1.3, 'weak direct path');
miniLobe(axMap, lx, ly-dy, segL, colRelay, 'NCR-relayed beam');
lx = xL(2) - 0.255*rx; ly = yL(2) - 0.055*ry;      %top-right: nulling mechanism
miniLobe(axMap, lx, ly,    segL, colBeam, 'desired FWA beam');
miniLobe(axMap, lx, ly-dy, segL, colIntf, 'interference (nulled)');

drawnow;
savefig(fig, fullfile(plotDir,'deployment_combined.fig'));
exportgraphics(fig, fullfile(plotDir,'deployment_combined.png'), 'Resolution', 260);
exportgraphics(fig, fullfile(plotDir,'deployment_combined.eps'));
fprintf('combined system-model figure saved to %s\n', fullfile(plotDir,'deployment_combined.png'));

%% ---------- local helpers ----------
function t = towerTop(D, cx, prj)
    %display point of the antenna head of the tower nearest world x = cx
    it = find(strcmp({D.key},'tower'));
    [~,ii] = min(abs([D(it).x] - cx)); d = D(it(ii));
    [sx,sy] = prj(d.x,d.y);
    t = [sx, sy + 0.88*d.H];
end
function beamLobe(ax, S, T, frac, c, p)
    %directional radiation lobe from S toward T: petal r = L*cos-shaped
    %pattern of sharpness p, spanning frac of the link length
    L = frac*norm(T - S); th0 = atan2(T(2)-S(2), T(1)-S(1));
    th = linspace(-pi, pi, 240);
    r = L*(0.5*(1+cos(th))).^p;
    patch(ax, S(1)+r.*cos(th+th0), S(2)+r.*sin(th+th0), c, ...
          'EdgeColor',c, 'FaceAlpha',0.26, 'LineWidth',1.5);
end
function p = bezPoint(p0, p2, lift, t)
    %point at parameter t on the quadratic arc used by arcPath
    ctrl = (p0 + p2)/2 + [0 lift];
    p = (1-t)^2*p0 + 2*(1-t)*t*ctrl + t^2*p2;
end
function arcPath(ax, p0, p2, lift, c, lw, ls)
    %quadratic bezier arc p0 -> p2 (control point lifted above the chord)
    ctrl = (p0 + p2)/2 + [0 lift];
    t = linspace(0,1,80)';
    b = (1-t).^2.*p0 + 2*(1-t).*t.*ctrl + t.^2.*p2;
    plot(ax, b(:,1), b(:,2), ls, 'Color',c, 'LineWidth',lw);
end
function miniLine(ax, x0, yc, segL, ls, c, lw, txt)
    plot(ax, [x0 x0+segL], [yc yc], ls, 'Color',c, 'LineWidth',lw);
    text(ax, x0+segL+0.012*diff(xlim(ax)), yc, txt, 'FontName','Times New Roman', ...
        'FontSize',10, 'VerticalAlignment','middle', 'Color',[0.15 0.15 0.15]);
end
function miniLobe(ax, x0, yc, segL, c, txt)
    %small horizontal petal glyph + text
    th = linspace(-pi, pi, 120);
    r = 1.15*segL*(0.5*(1+cos(th))).^6;
    patch(ax, x0+r.*cos(th), yc+0.5*r.*sin(th), c, 'EdgeColor',c, ...
          'FaceAlpha',0.3, 'LineWidth',1.1);
    text(ax, x0+1.15*segL+0.012*diff(xlim(ax)), yc, txt, 'FontName','Times New Roman', ...
        'FontSize',10, 'VerticalAlignment','middle', 'Color',[0.15 0.15 0.15]);
end
