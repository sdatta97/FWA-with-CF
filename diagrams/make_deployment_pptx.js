// System-model deployment figure for the FWA-vs-cellular paper, as an editable
// PowerPoint slide sized for a two-column (full-width) figure.
//
// 2.5D map: coverage as ground-plane ellipses (tall, near-top-down
// perspective so the scene fills the slide), objects as upright icon
// billboards from the Apache-2.0 MDI set in diagrams/icons. Shows the two
// paper mechanisms with inline path labels:
//   - NCR-relayed cellular support to a cell-edge in-car UE whose direct
//     beam is weak/stale (grey, dies at the x);
//   - the desired FWA beam to the cell-edge CPE and the other-site
//     interference the CPE array nulls (red, dies at the x).
//
//   npm install pptxgenjs
//   node diagrams/make_deployment_pptx.js diagrams/deployment_model.pptx
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.defineLayout({ name: "FIG", width: 13.33, height: 5.55 });
pres.layout = "FIG";

const F = "Times New Roman";
const C = {
  ink: "1A1A1A",
  beam: "0073BA",   // desired FWA beam
  intf: "B81A29",   // other-site interference (nulled)
  relay: "338C40",  // NCR-relayed path
  ncr: "0D2680",    // NCR unit
  blk: "8A8A8A",    // weak direct path
  ground: "F2F5F1", groundEdge: "BFCBBF",
};

/* ---------- icons ---------- */
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
function icon(key, gx, gy, H) { // base centred on the ground point (gx, gy)
  const W = H * ASPECT[key];
  sl.addImage({ data: DATA[key], x: gx - W / 2, y: gy - H, w: W, h: H });
  return W;
}
function beam(S, T, frac, color, halfW, transp) { // apex at S, widening toward T
  const dx = T[0] - S[0], dy = T[1] - S[1];
  const D = Math.hypot(dx, dy), L = D * frac;
  sl.addShape(pres.ShapeType.triangle, {
    x: S[0] + 0.5 * L * dx / D - halfW, y: S[1] + 0.5 * L * dy / D - L / 2,
    w: 2 * halfW, h: L, rotate: Math.atan2(-dx, dy) * 180 / Math.PI,
    fill: { color: color, transparency: transp === undefined ? 55 : transp },
    line: { color: color, width: 0.75, transparency: 25 },
  });
}
function seg(A, B, color, w) {
  const o = { x: Math.min(A[0], B[0]), y: Math.min(A[1], B[1]),
              w: Math.abs(B[0] - A[0]), h: Math.abs(B[1] - A[1]),
              line: { color: color, width: w } };
  if ((B[0] - A[0]) * (B[1] - A[1]) < 0) o.flipV = true;
  sl.addShape(pres.ShapeType.line, o);
}
function xMark(p, r, color) {
  seg([p[0] - r, p[1] - r], [p[0] + r, p[1] + r], color, 3);
  seg([p[0] - r, p[1] + r], [p[0] + r, p[1] - r], color, 3);
}

/* ---------- scene geometry (tall ellipses, near-top-down) ---------- */
const G1 = [3.43, 3.00], G2 = [9.90, 3.00];   // cell/ellipse centres
const RX = 3.23, RY = 1.55;                    // ground ellipse radii
const TOP1 = [3.43, 1.66], TOP2 = [9.90, 1.66]; // antenna heads (beam origins)
const UE = [1.30, 3.95];                       // in-car UE, cell edge
const NCR = [2.18, 1.72];                      // NCR unit on the office roof
const CPE = [6.66, 1.98];                      // rooftop CPE dish at the boundary

[G1, G2].forEach(function (g) {
  sl.addShape(pres.ShapeType.ellipse, { x: g[0] - RX, y: g[1] - RY, w: 2 * RX, h: 2 * RY,
    fill: { color: C.ground }, line: { color: C.groundEdge, width: 1, dashType: "dash" } });
});

// background scenery (far objects first), spread over the taller ground
icon("tree", 0.72, 2.95, 0.50);
icon("home", 2.45, 4.40, 0.52);
icon("th", 4.60, 3.55, 0.58);
icon("person", 4.18, 3.15, 0.40);
icon("th", 5.45, 4.12, 0.58);
icon("person", 7.60, 4.10, 0.40);
icon("th", 8.15, 3.60, 0.58);
icon("home", 11.05, 3.42, 0.52);
icon("tree", 12.55, 3.30, 0.50);
icon("home", 9.05, 4.42, 0.52);
// a second FWA household, to show the CPE fleet is not a single terminal
icon("th", 11.60, 2.55, 0.55);
icon("cpe", 11.60, 2.02, 0.33);

/* ---------- mechanism 1: NCR-relayed support at the cell edge ---------- */
beam(TOP1, UE, 0.52, C.blk, 0.15, 62);           // weak direct beam, dies short
xMark([2.32, 2.90], 0.10, C.blk);
beam(TOP1, NCR, 0.80, C.relay, 0.13, 48);        // gNB -> NCR
beam(NCR, UE, 0.82, C.relay, 0.13, 48);          // NCR -> UE

/* ---------- mechanism 2: interference nulling at the cell-edge CPE ---------- */
beam(TOP2, CPE, 0.90, C.beam, 0.17, 55);         // desired FWA beam
beam(TOP1, CPE, 0.60, C.intf, 0.17, 58);         // other-site interference
xMark([5.62, 1.82], 0.11, C.intf);

/* ---------- foreground actors, drawn over the beams ---------- */
icon("gnb", G1[0], G1[1] + 0.35, 1.75);
icon("gnb", G2[0], G2[1] + 0.35, 1.75);
icon("office", 2.18, 2.42, 0.72);                                     // NCR host
sl.addShape(pres.ShapeType.rect, { x: 2.00, y: 1.62, w: 0.36, h: 0.18,
  fill: { color: C.ncr }, line: { type: "none" } });
icon("office", 6.66, 2.75, 0.72);                                     // CPE host
icon("cpe", 6.66, 2.06, 0.42);
icon("car", UE[0], 4.13, 0.38);

/* ---------- labels (the user's edited set) ---------- */
txt("gNB 1", { x: 2.48, y: 3.48, w: 1.90, h: 0.24, fontSize: 13, bold: true,
  align: "center" });
txt("gNB 2", { x: 8.95, y: 3.48, w: 1.90, h: 0.24, fontSize: 13, bold: true,
  align: "center" });
txt("NCR", { x: 1.73, y: 1.30, w: 0.90, h: 0.22, fontSize: 12.5, bold: true,
  align: "center", color: C.ncr });
txt("FWA CPE (cell edge)", { x: 5.71, y: 1.32, w: 1.90, h: 0.22, fontSize: 12,
  bold: true, align: "center" });
txt("NCR-relayed path", { x: 0.28, y: 2.10, w: 1.45, h: 0.42, fontSize: 12,
  align: "center" });
txt("weak direct path", { x: 1.90, y: 3.06, w: 0.95, h: 0.44, fontSize: 12,
  align: "center" });
txt("other-site interference (nulled)", { x: 2.95, y: 1.10, w: 2.60, h: 0.24,
  fontSize: 12, align: "center" });
txt("desired FWA beam", { x: 7.90, y: 1.14, w: 1.75, h: 0.24, fontSize: 12,
  align: "center" });

/* ---------- icon legend ---------- */
[["gnb", "gNB site"], ["th", "Townhouse"], ["home", "Single-family home"],
 ["office", "Office park"], ["cpe", "Rooftop FWA CPE"], ["person", "Pedestrian UE"],
 ["car", "In-car UE"]]
.forEach(function (e, i) {
  const x = 0.52 + i * 1.87;
  const w = icon(e[0], x + 0.19 * ASPECT[e[0]], 5.38, 0.38);
  txt(e[1], { x: x + w + 0.11, y: 5.10, w: 1.87 - w - 0.15, h: 0.22, fontSize: 10.5 });
});

pres.writeFile({ fileName: process.argv[2] || "diagrams/deployment_model.pptx" })
  .then(function (f) { console.log("wrote " + f); });
