// 3Dmol.js renderer for the gromos.education viewer: spheres with ambient occlusion and outlines.
// 3Dmol.js arrives as the dmol_src string trait and runs from a script tag that the shell's
// deps.loadScript points at a blob: (or data:) URL, so nothing is fetched.

let loading = null;
async function load3Dmol(src, loadScript) {
  if (!window.$3Dmol) await (loading ??= loadScript(src).catch((e) => { loading = null; throw e; }));
  const $3Dmol = window.$3Dmol || window["3Dmol"];
  if (!$3Dmol) throw new Error("3Dmol.js ran but defined no $3Dmol");
  return $3Dmol;
}

function drawBox(viewer, [a, b, c]) {
  const P = (v) => ({ x: v[0] * a, y: v[1] * b, z: v[2] * c });
  for (let d = 0; d < 3; d++) {
    for (const [i, j] of [[0, 0], [0, 1], [1, 0], [1, 1]]) {
      const s = [i, j];
      s.splice(d, 0, 0);
      const e = [...s];
      e[d] = 1;
      viewer.addCylinder({ start: P(s), end: P(e), radius: 0.04, color: "#5c7cfa", fromCap: 1, toCap: 1 });
    }
  }
}

export async function create(stage, model, { loadScript }) {
  const $3Dmol = await load3Dmol(model.get("dmol_src"), loadScript);
  const viewer = $3Dmol.createViewer(stage, { backgroundColor: model.get("background") });
  const mol = viewer.addModel(model.get("pdb"), "pdb", { keepH: true });
  for (const [sel, style] of model.get("styles")) viewer.setStyle(sel, style);
  viewer.setViewStyle(model.get("view_style"));
  const box = model.get("box_angstrom");
  if (box.length === 3 && box.every((x) => x > 0)) drawBox(viewer, box);
  viewer.zoomTo();
  viewer.render();

  const n = mol.selectedAtoms({}).length;
  const r0 = model.get("radii"), c0 = model.get("colors");
  let radii = new Float32Array(r0.buffer, r0.byteOffset, r0.byteLength / 4);
  let rgb = new Uint8Array(c0.buffer, c0.byteOffset, c0.byteLength).slice();

  // One style per distinct radius; the colour of each atom is looked up by its index.
  function restyle() {
    const groups = new Map();
    for (let i = 0; i < n; i++) {
      const r = Math.round(radii[i] * 1e4) / 1e4;
      if (!groups.has(r)) groups.set(r, []);
      groups.get(r).push(i);
    }
    const colorfunc = (atom) => {
      const j = 3 * atom.index;
      return (rgb[j] << 16) | (rgb[j + 1] << 8) | rgb[j + 2];
    };
    for (const [radius, index] of groups) viewer.setStyle({ index }, { sphere: { radius, colorfunc } });
    viewer.render();
  }

  return {
    setPositions(p) {
      const xyz = new Array(n);
      for (let i = 0; i < n; i++) xyz[i] = [p[3 * i], p[3 * i + 1], p[3 * i + 2]];
      mol.setCoordinates([xyz], "array");
      viewer.setFrame(0);
      viewer.render();
    },
    setColors(c) { rgb = c.slice(); restyle(); },
    setRadii(r) { radii = r; restyle(); },
    resize() { viewer.resize(); viewer.render(); },
    dispose() { viewer.clear(); },
  };
}
