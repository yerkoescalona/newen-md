// three.js renderer for the gromos.education viewer: instanced, lit spheres with soft shadows,
// a table under a flat (2D) system, the box by its edges for a 3D one. Lengths are Angstrom.
// GROMOS (x, y, z) is drawn at scene (x, z, -y): the scene's y is up, and the map is a rotation.

function orbit(THREE, dom, camera, target, request) {
  // drag to rotate, wheel to zoom
  const sph = new THREE.Spherical().setFromVector3(camera.position.clone().sub(target));
  let drag = null;
  const apply = () => { camera.position.setFromSpherical(sph).add(target); camera.lookAt(target); request(); };
  dom.addEventListener("pointerdown", (e) => { drag = [e.clientX, e.clientY]; dom.setPointerCapture(e.pointerId); });
  dom.addEventListener("pointermove", (e) => {
    if (!drag) return;
    sph.theta -= (e.clientX - drag[0]) * 0.008;
    sph.phi = Math.min(Math.max(sph.phi - (e.clientY - drag[1]) * 0.008, 0.05), Math.PI - 0.05);
    drag = [e.clientX, e.clientY];
    apply();
  });
  dom.addEventListener("pointerup", () => { drag = null; });
  dom.addEventListener("wheel", (e) => { e.preventDefault(); sph.radius *= Math.exp(e.deltaY * 0.001); apply(); }, { passive: false });
  apply();
}

// A small lit room baked into an environment map, for reflections on the spheres.
function environment(THREE, renderer) {
  const room = new THREE.Scene();
  room.add(new THREE.Mesh(new THREE.BoxGeometry(10, 8, 10), new THREE.MeshStandardMaterial({ side: THREE.BackSide, color: 0x7f7f7f })));
  for (const [x, y, z, s] of [[-3, 3.9, 0, 2.5], [3, 3.9, 2, 1.5], [0, 2, -4.9, 3]]) {
    const lamp = new THREE.Mesh(new THREE.BoxGeometry(s, 0.1, s), new THREE.MeshBasicMaterial({ color: new THREE.Color(6, 6, 6) }));
    lamp.position.set(x, y, z);
    room.add(lamp);
  }
  const pmrem = new THREE.PMREMGenerator(renderer);
  const texture = pmrem.fromScene(room, 0.04).texture;
  pmrem.dispose();
  room.traverse((o) => { o.geometry?.dispose(); o.material?.dispose(); });
  return texture;
}

// Region the scene frames, in GROMOS coordinates: the box, else the atoms of the first frame.
function bounds(box, p, radii) {
  if (box.length === 3) return { lo: [0, 0, 0], hi: box };
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  const pad = Math.max(...radii, 1);
  for (let i = 0; i < radii.length; i++) {
    for (let d = 0; d < 3; d++) { lo[d] = Math.min(lo[d], p[3 * i + d] - pad); hi[d] = Math.max(hi[d], p[3 * i + d] + pad); }
  }
  return { lo, hi };
}

export function create(stage, model, { THREE }) {
  const f32 = (dv) => new Float32Array(dv.buffer, dv.byteOffset, dv.byteLength / 4);
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(stage.clientWidth, stage.clientHeight);
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.05;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.domElement.style.display = "block";
  stage.prepend(renderer.domElement);

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(model.get("background"));
  scene.environment = environment(THREE, renderer);
  scene.environmentIntensity = 0.9;
  const key = new THREE.DirectionalLight(0xffffff, 2.4);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  key.shadow.radius = 4;
  key.shadow.bias = -0.0004;
  const rim = new THREE.DirectionalLight(0x8fb4ff, 0.8);
  scene.add(key, key.target, rim, new THREE.HemisphereLight(0xbcd2ff, 0x1a1410, 0.35));

  let radii = f32(model.get("radii")), last = f32(model.get("positions"));
  const n = radii.length, flat = model.get("flat"), box = model.get("box_angstrom");
  const sphere = new THREE.SphereGeometry(1, 40, 24);
  const material = new THREE.MeshPhysicalMaterial({ roughness: 0.34, clearcoat: 0.3, clearcoatRoughness: 0.3 });
  const mesh = new THREE.InstancedMesh(sphere, material, n);
  mesh.castShadow = mesh.receiveShadow = true;
  mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  scene.add(mesh);

  const { lo, hi } = bounds(box, last, radii);
  const L = [0, 1, 2].map((d) => hi[d] - lo[d]);
  const span = flat ? Math.max(L[0], L[1]) : Math.max(...L);
  const c = new THREE.Vector3((lo[0] + hi[0]) / 2, flat ? 0 : (lo[2] + hi[2]) / 2, -(lo[1] + hi[1]) / 2);
  const owned = [sphere, material];
  const edgeMaterial = new THREE.LineBasicMaterial({ color: 0x5c7cfa, transparent: true, opacity: 0.7 });
  owned.push(edgeMaterial);
  if (flat) {
    const floorGeometry = new THREE.PlaneGeometry(L[0] * 1.6, L[1] * 1.6);
    const floorMaterial = new THREE.MeshStandardMaterial({ color: 0x1b2433, roughness: 0.9 });
    const floor = new THREE.Mesh(floorGeometry, floorMaterial);
    floor.rotation.x = -Math.PI / 2;
    floor.position.copy(c);
    floor.receiveShadow = true;
    scene.add(floor);
    owned.push(floorGeometry, floorMaterial);
  }
  if (box.length === 3) {
    // a 3D box by its edges; for a flat system, rails at the periodic edges of the table
    const edges = new THREE.EdgesGeometry(new THREE.BoxGeometry(L[0], flat ? 1e-4 : L[2], L[1]));
    const lines = new THREE.LineSegments(edges, edgeMaterial);
    lines.position.copy(c);
    if (flat) lines.position.y = 0.002 * span;
    scene.add(lines);
    owned.push(edges);
  }
  const sc = key.shadow.camera;
  sc.left = sc.bottom = -span * 1.1;
  sc.right = sc.top = span * 1.1;
  sc.near = span * 0.01;
  sc.far = span * 6;
  sc.updateProjectionMatrix();
  key.position.set(c.x - span * 0.6, c.y + span * 2.5, c.z + span * 0.8);
  key.target.position.copy(c);
  rim.position.set(c.x + span * 1.5, c.y + span * 0.4, c.z - span);
  const camera = new THREE.PerspectiveCamera(35, stage.clientWidth / stage.clientHeight, span * 1e-3, span * 50);
  camera.position.set(c.x + span * (flat ? 0 : 1.55), c.y + span * (flat ? 1.35 : 1.1), c.z + span * (flat ? 1.9 : 1.95));

  // Every change asks for one render, after the current task: several changes, one render.
  let queued = false;
  const request = () => {
    if (queued) return;
    queued = true;
    queueMicrotask(() => { queued = false; renderer.render(scene, camera); });
  };
  orbit(THREE, renderer.domElement, camera, c, request);

  const M = new THREE.Matrix4(), C = new THREE.Color();
  function setPositions(p) {
    last = p;
    for (let i = 0; i < n; i++) {
      const r = radii[i];
      mesh.setMatrixAt(i, M.makeScale(r, r, r).setPosition(p[3 * i], flat ? r : p[3 * i + 2], -p[3 * i + 1]));
    }
    mesh.instanceMatrix.needsUpdate = true;
    request();
  }
  function setColors(rgb) {
    for (let i = 0; i < n; i++) mesh.setColorAt(i, C.setRGB(rgb[3 * i] / 255, rgb[3 * i + 1] / 255, rgb[3 * i + 2] / 255, THREE.SRGBColorSpace));
    mesh.instanceColor.needsUpdate = true;
    request();
  }
  const cb = model.get("colors");
  setColors(new Uint8Array(cb.buffer, cb.byteOffset, cb.byteLength));
  setPositions(last);

  return {
    setPositions,
    setColors,
    setRadii(r) { radii = r; setPositions(last); },
    resize() {
      const w = stage.clientWidth, h = stage.clientHeight;
      if (!w || !h) return;
      renderer.setSize(w, h);
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      request();
    },
    dispose() {
      for (const o of owned) o.dispose();
      mesh.dispose();
      scene.environment.dispose();
      renderer.dispose();
      renderer.domElement.remove();
    },
  };
}
