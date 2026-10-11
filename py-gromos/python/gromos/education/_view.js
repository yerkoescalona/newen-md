// gromos.education.viewer: spheres beside energy and temperature strip charts.
// The renderer module (three.js or 3Dmol.js) and three.js itself arrive as string traits and are
// imported from blob: URLs (data: URLs if the page refuses blobs), so nothing is fetched from the
// network. A renderer module exports create(stage, model, deps) -> { setPositions(f32), setColors(u8),
// setRadii(f32), resize(), dispose() }. Coordinates arrive as float32 bytes in Angstrom, colours as
// uint8 RGB; the energy trait is a list of rows [t, kinetic, potential, total, T].

// window.__gromos_view_load records which route (blob or data) loaded each source, for diagnosis.
const modules = new Map(), routes = (window.__gromos_view_load ??= {});
const base64 = (bytes, s = "") => {
  for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(s);
};
function viaUrl(what, src, load) {
  const done = (route) => (out) => { routes[what] = route; return out; };
  return load(URL.createObjectURL(new Blob([src], { type: "text/javascript" }))).then(done("blob"), (blobError) => {
    console.warn(`gromos view: ${what} from a blob: URL failed, retrying as a data: URL`, blobError);
    return load(`data:text/javascript;base64,${base64(new TextEncoder().encode(src))}`).then(done("data"), (e) => {
      throw Object.assign(new Error(`${e.message} (data: URL; blob: URL: ${blobError.message})`), { name: e.name });
    });
  });
}
function importSource(key, what, src) {
  if (!modules.has(key)) modules.set(key, viaUrl(what, src, (url) => import(url)).catch((e) => { modules.delete(key); throw e; }));
  return modules.get(key);
}
const loadScript = (src) => viaUrl("dmol", src, (url) => new Promise((ok, bad) => {
  const s = Object.assign(document.createElement("script"), { src: url, onload: ok });
  s.onerror = () => { s.remove(); bad(new Error(`the page refused the script tag (${url.slice(0, 5)} URL)`)); };
  document.head.appendChild(s);
}));

const MONO = "var(--jp-code-font-family, ui-monospace, SFMono-Regular, Menlo, monospace)";
const CSS = `
.gv { --gv-ink: var(--jp-ui-font-color1, #141a24); --gv-ink2: var(--jp-ui-font-color2, #546072);
  --gv-rule: var(--jp-border-color2, #d4dae3); --gv-panel: var(--jp-layout-color1, #ffffff);
  --gv-accent: var(--jp-brand-color1, #3b5bdb); --gv-ke: #e8590c; --gv-pe: #1c7ed6; --gv-temp: #ae3ec9;
  --gv-mono: ${MONO};
  display: flex; flex-wrap: wrap; gap: 12px; align-items: flex-start; color: var(--gv-ink);
  font: 13px var(--jp-ui-font-family, system-ui, sans-serif); }
.gv.gv-dark { --gv-ink: var(--jp-ui-font-color1, #e4e9f1); --gv-ink2: var(--jp-ui-font-color2, #98a5b8);
  --gv-rule: var(--jp-border-color2, #2a3442); --gv-panel: var(--jp-layout-color1, #161d28);
  --gv-accent: var(--jp-brand-color1, #8ea6ff); --gv-ke: #ff922b; --gv-pe: #4dabf7; --gv-temp: #da77f2; }
.gv-stage { position: relative; border: 1px solid var(--gv-rule); border-radius: 8px; overflow: hidden; touch-action: none; }
.gv-legend { position: absolute; right: 10px; top: 10px; z-index: 1; display: grid; gap: 4px; justify-items: end;
  padding: 6px 8px; border-radius: 6px; background: rgba(12, 18, 28, 0.55); color: #c9d3e2;
  font: 11px ${MONO}; pointer-events: none; }
.gv-legend:empty { display: none; }
.gv-legend div { width: 120px; height: 8px; border-radius: 4px; }
.gv-legend p { display: flex; justify-content: space-between; width: 120px; margin: 0; }
.gv-transport { display: flex; gap: 8px; align-items: center; margin-top: 6px; }
.gv-transport input { flex: 1; min-width: 0; accent-color: var(--gv-accent); }
.gv button { font: 600 13px var(--jp-ui-font-family, system-ui, sans-serif); color: var(--gv-ink); background: var(--gv-panel); min-width: 5em;
  border: 1px solid var(--gv-rule); border-radius: 6px; padding: 4px 12px; cursor: pointer; }
.gv button:hover { border-color: var(--gv-accent); }
.gv-time { font: 12px ${MONO}; color: var(--gv-ink2); min-width: 9em; text-align: right; }
.gv-side { display: grid; gap: 8px; width: 280px; }
.gv-readouts { display: grid; grid-template-columns: repeat(3, 1fr); gap: 6px; }
.gv-box { background: var(--gv-panel); border: 1px solid var(--gv-rule); border-radius: 6px; padding: 6px 8px; }
.gv-box span, .gv-keys { font-size: 11px; color: var(--gv-ink2); }
.gv-box b { display: block; font: 600 14px ${MONO}; font-variant-numeric: tabular-nums; }
.gv-box h4 { margin: 0 0 2px; font-size: 12px; }
.gv-box canvas { display: block; width: 100%; height: 130px; }
.gv-keys i { display: inline-block; width: 12px; height: 3px; margin: 0 4px 2px 0; vertical-align: middle; }
.gv-keys span + span { margin-left: 10px; }
.gv-error { position: absolute; inset: 0; margin: 0; padding: 12px; overflow: auto; white-space: pre-wrap;
  font: 12px ${MONO}; color: var(--gv-ink); background: var(--gv-panel); }`;

function el(tag, cls, html) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html) e.innerHTML = html;
  return e;
}

function isDark() {
  const lab = document.body.dataset.jpThemeLight;
  return lab ? lab === "false" : matchMedia("(prefers-color-scheme: dark)").matches;
}

const fmt = (v) => (Math.abs(v) >= 100 ? v.toFixed(1) : v.toPrecision(4));

function niceTicks(lo, hi, n) {
  const raw = (hi - lo) / n, mag = 10 ** Math.floor(Math.log10(raw)), err = raw / mag;
  const step = (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1) * mag;
  const dec = Math.max(0, -Math.floor(Math.log10(step))), out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) out.push(v.toFixed(dec));
  return out;
}

// A strip chart of rows[*][cols]: the whole run faint, the part up to row k solid, a dot at row k.
function chart(canvas, rows, cols, k, css, ref) {
  const dpr = Math.min(devicePixelRatio, 2), w = canvas.clientWidth, h = canvas.clientHeight;
  if (!w) return;
  if (canvas.width !== Math.round(w * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const g = canvas.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  if (!rows.length) return;
  const t0 = rows[0][0], t1 = Math.max(rows[rows.length - 1][0], t0 + 1e-9), p = { l: 48, r: 8, t: 6, b: 18 };
  let lo = Infinity, hi = -Infinity;
  for (const [c] of cols) for (const r of rows) { lo = Math.min(lo, r[c]); hi = Math.max(hi, r[c]); }
  if (ref != null) { lo = Math.min(lo, ref); hi = Math.max(hi, ref); }
  const pad = (hi - lo || Math.abs(hi) || 1) * 0.08;
  lo -= pad; hi += pad;
  const X = (t) => p.l + ((t - t0) / (t1 - t0)) * (w - p.l - p.r);
  const Y = (v) => p.t + (1 - (v - lo) / (hi - lo)) * (h - p.t - p.b);
  g.font = `10px ${css("--gv-mono")}`;
  g.fillStyle = css("--gv-ink2");
  g.strokeStyle = css("--gv-rule");
  g.lineWidth = 1;
  g.textAlign = "right";
  g.textBaseline = "middle";
  for (const v of niceTicks(lo, hi, 4)) {
    g.beginPath(); g.moveTo(p.l, Y(+v)); g.lineTo(w - p.r, Y(+v)); g.stroke();
    g.fillText(v, p.l - 5, Y(+v));
  }
  g.textAlign = "center";
  g.textBaseline = "alphabetic";
  for (const v of niceTicks(t0, t1, 4)) g.fillText(v, X(+v), h - 4);
  if (ref != null) {
    g.setLineDash([4, 4]); g.strokeStyle = css("--gv-ink2");
    g.beginPath(); g.moveTo(p.l, Y(ref)); g.lineTo(w - p.r, Y(ref)); g.stroke(); g.setLineDash([]);
  }
  const line = (c, upto) => {
    g.beginPath();
    for (let i = 0; i <= upto; i++) (i ? g.lineTo : g.moveTo).call(g, X(rows[i][0]), Y(rows[i][c]));
    g.stroke();
  };
  for (const [c, colour] of cols) {
    g.strokeStyle = g.fillStyle = css(colour);
    g.globalAlpha = 0.25; g.lineWidth = 1.2; line(c, rows.length - 1);
    g.globalAlpha = 1; g.lineWidth = 1.8; line(c, k);
    g.beginPath(); g.arc(X(rows[k][0]), Y(rows[k][c]), 3, 0, 2 * Math.PI); g.fill();
  }
}

function legend(box, ramp, vmax) {
  if (!(vmax > 0)) { box.innerHTML = ""; return; }
  const stops = ramp.map((c, i) => `${c} ${(100 * i) / (ramp.length - 1)}%`).join(",");
  box.innerHTML = `<span>speed (nm/ps)</span><div style="background:linear-gradient(90deg,${stops})"></div>
    <p><span>0</span><span>${vmax.toFixed(2)}</span></p>`;
}

function layout(root, width, height) {
  const stageWidth = Math.max(240, width - 292);
  root.innerHTML = `
    <div><div class="gv-stage" style="width:${stageWidth}px;height:${height}px"><div class="gv-legend"></div></div>
      <div class="gv-transport"><button type="button">Play</button>
        <input type="range" min="0" max="0" value="0" disabled><span class="gv-time">t = –</span></div></div>
    <div class="gv-side"><div class="gv-readouts">
        <div class="gv-box"><span>Temperature</span><b>–</b></div>
        <div class="gv-box"><span>Kinetic</span><b>–</b></div>
        <div class="gv-box"><span>Potential</span><b>–</b></div></div>
      <div class="gv-box"><h4>Energy (kJ/mol)</h4><div class="gv-keys">
        <span><i style="background:var(--gv-ke)"></i>kinetic</span><span><i style="background:var(--gv-pe)"></i>potential</span>
        <span><i style="background:var(--gv-ink)"></i>total</span></div><canvas></canvas></div>
      <div class="gv-box"><h4>Temperature (K)</h4><canvas></canvas></div></div>`;
  const [energyCanvas, tempCanvas] = root.querySelectorAll("canvas");
  return {
    stage: root.querySelector(".gv-stage"), legend: root.querySelector(".gv-legend"), button: root.querySelector("button"),
    slider: root.querySelector("input"), clock: root.querySelector(".gv-time"),
    readouts: root.querySelectorAll(".gv-readouts b"), energyCanvas, tempCanvas,
  };
}

async function render({ model, el: host }) {
  if (!document.getElementById("gv-style")) document.head.appendChild(Object.assign(el("style"), { id: "gv-style", textContent: CSS }));
  const root = el("div", "gv");
  host.appendChild(root);
  const ui = layout(root, model.get("width"), model.get("height"));
  let scene = null; // until the renderer has loaded, only the charts and readouts draw

  const f32 = (dv) => new Float32Array(dv.buffer, dv.byteOffset, dv.byteLength / 4);
  const u8 = (dv) => new Uint8Array(dv.buffer, dv.byteOffset, dv.byteLength);
  const n = model.get("radii").byteLength / 4;
  const css = (name) => getComputedStyle(root).getPropertyValue(name).trim();
  // Lets the page count drawn updates ("applied") and replayed frames ("played").
  const signal = (kind) => host.dispatchEvent(new CustomEvent("gromos-view-frame", { bubbles: true, detail: kind }));

  let pending = null, pendingColors = null, traj = null, nframes = 0, colorTraj = null;
  let frame = -1, cursor = -1, last = 0, dirty = true, alive = true;
  const showFrame = (k) => {
    frame = cursor = k;
    ui.slider.value = String(k);
    scene?.setPositions(traj.subarray(3 * n * k, 3 * n * (k + 1)));
    if (colorTraj) scene?.setColors(colorTraj.subarray(3 * n * k, 3 * n * (k + 1)));
    dirty = true;
  };
  const setPlaying = (p) => { model.set("playing", p); model.save_changes(); };
  const theme = () => { root.classList.toggle("gv-dark", isDark()); dirty = true; };

  function onTrajectory() {
    const t = model.get("trajectory"), nf = t && t.byteLength ? t.byteLength / 12 / n : 0;
    traj = nf ? f32(t) : null;
    nframes = nf;
    frame = -1;
    ui.slider.max = String(Math.max(nf - 1, 0));
    ui.slider.disabled = !traj;
  }
  function charts() {
    const rows = model.get("energy");
    const k = Math.min(cursor < 0 ? rows.length - 1 : cursor, rows.length - 1);
    chart(ui.energyCanvas, rows, [[1, "--gv-ke"], [2, "--gv-pe"], [3, "--gv-ink"]], k, css);
    chart(ui.tempCanvas, rows, [[4, "--gv-temp"]], k, css);
    if (k < 0) return;
    const [t, kinetic, potential, , T] = rows[k];
    ui.clock.textContent = `t = ${t.toFixed(3)} ps`;
    [`${T.toFixed(1)} K`, fmt(kinetic), fmt(potential)].forEach((s, i) => { ui.readouts[i].textContent = s; });
  }

  const onColorTrajectory = () => {
    const c = model.get("color_trajectory");
    colorTraj = c && c.byteLength ? u8(c) : null;
  };
  model.on("change:positions", () => { pending = f32(model.get("positions")); });
  model.on("change:colors", () => { pendingColors = u8(model.get("colors")); });
  model.on("change:radii", () => scene?.setRadii(f32(model.get("radii"))));
  model.on("change:vmax", () => legend(ui.legend, model.get("ramp"), model.get("vmax")));
  model.on("change:color_trajectory", onColorTrajectory);
  model.on("change:trajectory", onTrajectory);
  model.on("change:energy", () => { dirty = true; });
  model.on("change:playing", () => { ui.button.textContent = model.get("playing") ? "Pause" : "Play"; });
  ui.button.addEventListener("click", () => setPlaying(!model.get("playing")));
  ui.slider.addEventListener("input", () => { setPlaying(false); if (traj) showFrame(+ui.slider.value); });
  const observer = new MutationObserver(theme);
  observer.observe(document.body, { attributes: true, attributeFilter: ["data-jp-theme-light"] });
  const resizer = new ResizeObserver(() => scene?.resize());
  resizer.observe(ui.stage);
  theme();
  legend(ui.legend, model.get("ramp"), model.get("vmax"));
  onColorTrajectory();
  onTrajectory();
  if (model.get("positions")?.byteLength) pending = f32(model.get("positions"));
  ui.button.textContent = model.get("playing") ? "Pause" : "Play";

  function applyPending() { // newest coordinates win: at most one redraw per repaint
    if (pendingColors) { scene.setColors(pendingColors); pendingColors = null; }
    if (!pending) return false;
    scene.setPositions(pending);
    pending = null; cursor = -1; dirty = true;
    signal("applied");
    return true;
  }
  (function loop(now) {
    if (!alive) return;
    requestAnimationFrame(loop);
    if (!(scene && applyPending()) && traj && model.get("playing")) {
      const fps = model.get("fps");
      if (fps <= 0 || now - last >= 1000 / fps) {
        last = now;
        showFrame((frame + 1) % nframes);
        signal("played");
      }
    }
    if (dirty) { dirty = false; charts(); }
  })(0);

  // A failure while loading or creating the scene is shown in the stage and in `status`.
  const setStatus = (text) => { model.set("status", text); model.save_changes(); };
  (async (name, step) => {
    const deps = {
      loadScript: (src) => { step = "loading 3Dmol.js"; return loadScript(src).then(() => { step = "creating the scene"; }); },
    };
    try {
      step = "loading three.js";
      if (name === "three") deps.THREE = await importSource("three", "three", model.get("three_src"));
      step = "loading the renderer module";
      const { create } = await importSource(`renderer:${name}`, "renderer", model.get("renderer_src"));
      step = "creating the scene";
      const created = await create(ui.stage, model, deps);
      if (!alive) return created.dispose();
      step = "drawing the first frame";
      scene = created;
      applyPending();
      setStatus("ready");
    } catch (e) {
      scene = null;
      const text = `${step} failed: ${e?.name ?? "Error"}: ${e?.message ?? e}`;
      console.error(`gromos view: ${text}`, e);
      ui.stage.append(Object.assign(el("pre", "gv-error"), { textContent: text }));
      setStatus(text);
    }
  })(model.get("renderer"), "");

  return () => { alive = false; observer.disconnect(); resizer.disconnect(); scene?.dispose(); };
}

export default { render };
