// System-model deployment figure for the FWA-vs-cellular paper, as an editable
// PowerPoint slide sized for a two-column (full-width) figure.
//
// Replaces the pseudo-3D scene of make_deployment_combined.m. The geometry is a
// 2.5D map: coverage drawn as ground-plane ellipses, objects as upright icon
// billboards (same trick, and the same Apache-2.0 MDI icon set in
// diagrams/icons, as make_deployment_3d.m). The ground plane is what buys the
// clarity - the NCR host gets a real lateral offset from the direct path, which
// in a side elevation collapses onto it.
//
// Shows the two paper mechanisms plus the orthogonal band split, qualitatively:
//   - NCR-relayed cellular support to a cell-edge in-car UE whose direct beam
//     is stale at 40 km/h;
//   - other-site interference nulled by the cell-edge CPE array;
//   - cellular and FWA on disjoint sub-bands, so neither interferes with the
//     other (inset; the measured split lives in the results figures).
//
//   npm install pptxgenjs
//   node diagrams/make_deployment_pptx.js plots/deployment_model.pptx
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.defineLayout({ name: "FIG", width: 13.33, height: 5.55 });
pres.layout = "FIG";

const F = "Times New Roman";
const C = {
  ink: "1A1A1A", mute: "5A5A5A", faint: "9A9A9A",
  beam: "0073BA",   // FWA sub-band: desired beam, CPE pattern
  intf: "B81A29",   // other-site interference (nulled)
  relay: "338C40",  // NCR-relayed path
  ncr: "0D2680",    // NCR unit
  blk: "8A8A8A",    // weak / stale direct path
  cell: "5B6B7A",   // cellular sub-band
  ground: "F2F5F1", groundEdge: "BFCBBF",
  panel: "D8D8D8",
};

/* ---------- icons: the same Apache-2.0 MDI set the MATLAB figures use ---------- */
const ICODIR = path.join(__dirname, "icons");
const ASPECT = {}, DATA = {};
fs.readFileSync(path.join(ICODIR, "manifest.csv"), "utf8")
  .trim().split("\n").slice(1).forEach(function (ln) {
    const c = ln.split(",");
    ASPECT[c[0]] = parseFloat(c[2]);
    DATA[c[0]] = "image/png;base64," +
      fs.readFileSync(path.join(ICODIR, c[0] + ".png")).toString("base64");
  });

const sl = pres.addSlide();
sl.background = { color: "FFFFFF" };

/* ---------- primitives ---------- */
function txt(t, o) { sl.addText(t, Object.assign({ fontFace: F, color: C.ink, margin: 0 }, o)); }
// icon billboard: base centred on the ground point (gx, gy), drawn upright
function icon(key, gx, gy, H) {
  const W = H * ASPECT[key];
  sl.addImage({ data: DATA[key], x: gx - W / 2, y: gy - H, w: W, h: H });
  return W;
}
// directional beam: filled isoceles triangle, apex at S, widening toward T
function beam(S, T, frac, color, halfW, transp) {
  const dx = T[0] - S[0], dy = T[1] - S[1];
  const D = Math.hypot(dx, dy), L = D * frac;
  sl.addShape(pres.ShapeType.triangle, {
    x: S[0] + 0.5 * L * dx / D - halfW, y: S[1] + 0.5 * L * dy / D - L / 2,
    w: 2 * halfW, h: L, rotate: Math.atan2(-dx, dy) * 180 / Math.PI,
    fill: { color: color, transparency: transp === undefined ? 55 : transp },
    line: { color: color, width: 0.75, transparency: 25 },
  });
}
function seg(A, B, color, w, dash) {
  const o = { x: Math.min(A[0], B[0]), y: Math.min(A[1], B[1]),
              w: Math.abs(B[0] - A[0]), h: Math.abs(B[1] - A[1]),
              line: { color: color, width: w } };
  if (dash) o.line.dashType = dash;
  if ((B[0] - A[0]) * (B[1] - A[1]) < 0) o.flipV = true;
  sl.addShape(pres.ShapeType.line, o);
}
function xMark(p, r, color) {
  seg([p[0] - r, p[1] - r], [p[0] + r, p[1] + r], color, 2.5);
  seg([p[0] - r, p[1] + r], [p[0] + r, p[1] - r], color, 2.5);
}

/* ---------- scene geometry ---------- */
const G1 = [3.40, 3.05], G2 = [9.30, 3.05];   // gNB ground points
const TOP1 = [3.40, 1.90], TOP2 = [9.30, 1.90]; // antenna heads (beam origins)
const UE = [1.30, 3.40];                       // in-car UE, cell edge
const NCR = [2.30, 1.84];                      // NCR unit on the office roof
const CPE = [6.35, 2.23];                      // rooftop CPE dish at the boundary

// coverage: two ground-plane ellipses meeting exactly at the cell boundary
[G1, G2].forEach(function (g) {
  sl.addShape(pres.ShapeType.ellipse, { x: g[0] - 2.95, y: g[1] - 0.95, w: 5.90, h: 1.90,
    fill: { color: C.ground }, line: { color: C.groundEdge, width: 1, dashType: "dash" } });
});
seg([6.35, 2.05], [6.35, 4.02], C.faint, 1, "sysDot");
txt("cell boundary", { x: 5.85, y: 4.04, w: 1.00, h: 0.18, fontSize: 9, italic: true,
  align: "center", color: C.mute });

// background scenery, sparse (far objects first)
icon("tree", 0.75, 3.15, 0.38);
icon("home", 2.60, 3.80, 0.40);
icon("th", 4.55, 3.60, 0.44);
icon("person", 4.20, 3.28, 0.30);
icon("th", 8.10, 3.62, 0.44);
icon("person", 9.75, 3.85, 0.30);
icon("home", 10.75, 3.50, 0.40);
icon("tree", 11.85, 3.25, 0.38);
// a second FWA household, to show the CPE fleet is not a single terminal
icon("th", 10.30, 2.72, 0.42);
icon("cpe", 10.30, 2.30, 0.26);

/* ---------- mechanism 1: NCR-relayed support at the cell edge ---------- */
beam(TOP1, UE, 0.55, C.blk, 0.115, 62);          // stale direct beam, dies short
xMark([2.10, 2.82], 0.085, C.blk);
beam(TOP1, NCR, 0.82, C.relay, 0.100, 48);       // gNB -> NCR
beam(NCR, UE, 0.80, C.relay, 0.100, 48);         // NCR -> UE

/* ---------- mechanism 2: interference nulling at the cell-edge CPE ---------- */
beam(TOP2, CPE, 0.92, C.beam, 0.130, 55);        // desired FWA beam
beam(TOP1, CPE, 0.62, C.intf, 0.130, 58);        // other-site interference
xMark([5.47, 2.13], 0.095, C.intf);

/* ---------- foreground actors, drawn over the beams ---------- */
icon("gnb", G1[0], G1[1], 1.30);
icon("gnb", G2[0], G2[1], 1.30);
icon("office", 2.30, 2.45, 0.55);                                     // NCR host
sl.addShape(pres.ShapeType.rect, { x: 2.15, y: 1.75, w: 0.30, h: 0.15,
  fill: { color: C.ncr }, line: { type: "none" } });
icon("office", 6.35, 2.95, 0.55);                                     // CPE host
icon("cpe", 6.35, 2.40, 0.34);
icon("car", UE[0], 3.55, 0.30);

/* ---------- labels ---------- */
txt("gNB site 1", { x: 2.45, y: 3.10, w: 1.90, h: 0.20, fontSize: 11, bold: true,
  align: "center" });
txt("gNB site 2", { x: 8.35, y: 3.10, w: 1.90, h: 0.20, fontSize: 11, bold: true,
  align: "center" });
txt("NCR", { x: 1.85, y: 1.53, w: 0.90, h: 0.19, fontSize: 11, bold: true,
  align: "center", color: C.ncr });
txt("in-car UE", { x: 0.65, y: 3.60, w: 1.30, h: 0.19, fontSize: 10, align: "center" });
txt("FWA CPE (cell edge)", { x: 5.50, y: 1.76, w: 1.70, h: 0.19, fontSize: 10, bold: true,
  align: "center" });
txt("stale CSI at 40 km/h —\nthe relay hop bypasses it",
  { x: 0.50, y: 1.95, w: 1.25, h: 0.42, fontSize: 9, italic: true, color: C.mute,
    lineSpacingMultiple: 0.95 });
txt("other-site interference\nnulled by the CPE array",
  { x: 4.30, y: 2.32, w: 1.65, h: 0.40, fontSize: 9, italic: true, align: "right",
    color: C.intf, lineSpacingMultiple: 0.95 });

/* ---------- inset: the two services sit on orthogonal sub-bands ---------- */
sl.addShape(pres.ShapeType.rect, { x: 9.62, y: 0.34, w: 3.36, h: 1.12,
  fill: { color: "FFFFFF" }, line: { color: C.panel, width: 1 } });
txt("Orthogonal sub-bands", { x: 9.76, y: 0.43, w: 3.08, h: 0.20, fontSize: 11, bold: true });
sl.addShape(pres.ShapeType.rect, { x: 9.80, y: 0.70, w: 1.42, h: 0.28,
  fill: { color: C.cell }, line: { color: "FFFFFF", width: 1 } });
sl.addShape(pres.ShapeType.rect, { x: 11.22, y: 0.70, w: 1.56, h: 0.28,
  fill: { color: C.beam }, line: { color: "FFFFFF", width: 1 } });
txt("cellular", { x: 9.80, y: 0.76, w: 1.42, h: 0.18, fontSize: 10, bold: true,
  align: "center", color: "FFFFFF" });
txt("FWA", { x: 11.22, y: 0.76, w: 1.56, h: 0.18, fontSize: 10, bold: true,
  align: "center", color: "FFFFFF" });
txt("frequency  →", { x: 11.28, y: 1.00, w: 1.50, h: 0.16, fontSize: 8.5, italic: true,
  align: "right", color: C.mute });
txt("the two services never share a resource — no cross-interference",
  { x: 9.76, y: 1.20, w: 3.08, h: 0.20, fontSize: 8.5, italic: true, color: C.mute });

/* ---------- path legend ---------- */
[[C.blk, "weak direct path"], [C.relay, "NCR-relayed path"],
 [C.beam, "desired FWA beam"], [C.intf, "other-site interference (nulled)"]]
.forEach(function (e, i) {
  const x = 0.50 + i * 3.05;
  sl.addShape(pres.ShapeType.triangle, { x: x, y: 4.30, w: 0.115, h: 0.30, rotate: 90,
    fill: { color: e[0], transparency: 45 }, line: { color: e[0], width: 0.75 } });
  txt(e[1], { x: x + 0.38, y: 4.34, w: 2.55, h: 0.22, fontSize: 10 });
});

/* ---------- icon legend ---------- */
[["gnb", "gNB site"], ["th", "Townhouse"], ["home", "Single-family home"],
 ["office", "Office park"], ["cpe", "Rooftop FWA CPE"], ["person", "Pedestrian UE"],
 ["car", "In-car UE"]]
.forEach(function (e, i) {
  const x = 0.50 + i * 1.78;
  const w = icon(e[0], x + 0.17 * ASPECT[e[0]], 5.15, 0.34);
  txt(e[1], { x: x + w + 0.10, y: 4.90, w: 1.78 - w - 0.14, h: 0.20, fontSize: 9.5 });
});

pres.writeFile({ fileName: process.argv[2] || "deployment_model.pptx" })
  .then(function (f) { console.log("wrote " + f); });
