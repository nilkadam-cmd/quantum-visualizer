// --- Custom mass toggle ---
function toggleCustomMass() {
    const massSelect = document.getElementById("mass");
    const customMassContainer = document.getElementById("custom-mass-container");
    customMassContainer.style.display = (massSelect.value === "custom") ? "block" : "none";
}

// --- Physics ---
function getInputsNormalized() {
    const nx = parseInt(document.getElementById("nx").value);
    const ny = parseInt(document.getElementById("ny").value);
    const nz = parseInt(document.getElementById("nz").value);

    const Lx = parseFloat(document.getElementById("Lx").value) * 1e-9;
    const Ly = parseFloat(document.getElementById("Ly").value) * 1e-9;
    const Lz = parseFloat(document.getElementById("Lz").value) * 1e-9;

    let m;
    const massOpt = document.getElementById("mass").value;
    if (massOpt === "electron") m = 9.109e-31;
    else if (massOpt === "proton") m = 1.673e-27;
    else m = parseFloat(document.getElementById("customMass").value);

    return { nx, ny, nz, Lx, Ly, Lz, m };
}

function energyLevel(nx, ny, nz, Lx, Ly, Lz, m) {
    const hbar = 1.054e-34;
    const pi = Math.PI;
    return (pi * pi * hbar * hbar / (2 * m)) * ((nx * nx) / (Lx * Lx) + (ny * ny) / (Ly * Ly) + (nz * nz) / (Lz * Lz));
}

function psi_at(nx, ny, nz, Lx, Ly, Lz, x, y, z) {
    const pi = Math.PI;
    const pref = Math.sqrt(8 / (Lx * Ly * Lz));
    return pref * Math.sin(nx * pi * x / Lx) * Math.sin(ny * pi * y / Ly) * Math.sin(nz * pi * z / Lz);
}

// --- Calculate ---
function calculate() {
    const { nx, ny, nz, Lx, Ly, Lz, m } = getInputsNormalized();
    const E = energyLevel(nx, ny, nz, Lx, Ly, Lz, m);

    const psi = psi_at(nx, ny, nz, Lx, Ly, Lz, Lx / 2, Ly / 2, Lz / 2);
    const density = psi * psi;

    document.getElementById("wavefunction").innerHTML = `ψ(Lx/2, Ly/2, Lz/2) = ${psi.toExponential(3)}`;
    document.getElementById("probability").innerHTML = `|ψ|² = ${density.toExponential(3)}`;
    document.getElementById("energy").innerHTML = `E = ${E.toExponential(3)} J (${(E / 1.602e-19).toFixed(3)} eV)`;

    MathJax.typeset();
}

// --- Grid & Visualisation ---
function computeGridAndDensity(N = 30) {
    const { nx, ny, nz, Lx, Ly, Lz } = getInputsNormalized();
    const xs = new Float64Array(N), ys = new Float64Array(N), zs = new Float64Array(N);
    for (let i = 0; i < N; i++) { xs[i] = (i + 0.5) / N * Lx; ys[i] = (i + 0.5) / N * Ly; zs[i] = (i + 0.5) / N * Lz; }
    const density = new Float64Array(N * N * N);
    let idx = 0;
    for (let k = 0; k < N; k++) for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
        const p = psi_at(nx, ny, nz, Lx, Ly, Lz, xs[i], ys[j], zs[k]);
        density[idx++] = p * p;
    }
    return { N, xs, ys, zs, density };
}

function plotSlice(grid, axis = 'z', sliceIndex = null) {
    const { N, xs, ys, zs, density } = grid;
    if (sliceIndex === null) sliceIndex = Math.floor(N / 2);
    let xvals, yvals, zslice;
    if (axis === 'z') {
        xvals = xs; yvals = ys;
        zslice = Array.from({ length: N }, (_, j) => Array.from({ length: N }, (_, i) => density[i + j * N + sliceIndex * N * N]));
    }
    else if (axis === 'y') {
        xvals = xs; yvals = zs;
        zslice = Array.from({ length: N }, (_, k) => Array.from({ length: N }, (_, i) => density[i + sliceIndex * N + k * N * N]));
    }
    else {
        xvals = ys; yvals = zs;
        zslice = Array.from({ length: N }, (_, k) => Array.from({ length: N }, (_, j) => density[sliceIndex + j * N + k * N * N]));
    }
    Plotly.newPlot("slicePlot", [{
        z: zslice, x: xvals, y: yvals, type: "heatmap", colorscale: "Viridis"
    }],
        { title: `|ψ|² slice (${axis})` }, { responsive: true });
}

function expandCoords(xs, ys, zs, mode) { const N = xs.length; const arr = []; for (let k = 0; k < N; k++)for (let j = 0; j < N; j++)for (let i = 0; i < N; i++) { if (mode === "x") arr.push(xs[i]); if (mode === "y") arr.push(ys[j]); if (mode === "z") arr.push(zs[k]); } return arr; }

function plot3D(grid, mode = 'isosurface') {
    const { N, xs, ys, zs, density } = grid; const values = Array.from(density); let trace;
    if (mode === 'isosurface') { trace = { type: "isosurface", x: expandCoords(xs, ys, zs, "x"), y: expandCoords(xs, ys, zs, "y"), z: expandCoords(xs, ys, zs, "z"), value: values, isomin: Math.max(...values) * 0.15, isomax: Math.max(...values) * 0.9, caps: { x: { show: false }, y: { show: false }, z: { show: false } } }; }
    else { trace = { type: "volume", x: expandCoords(xs, ys, zs, "x"), y: expandCoords(xs, ys, zs, "y"), z: expandCoords(xs, ys, zs, "z"), value: values, opacity: 0.1, isomin: 0, isomax: Math.max(...values) }; }
    Plotly.newPlot("plot3d", [trace], { margin: { l: 0, r: 0, b: 0, t: 30 } }, { responsive: true });
}

// --- UI Hooks ---
document.getElementById("plotBtn").addEventListener("click", () => {
    const N = parseInt(document.getElementById("gridN").value) || 30;
    const mode = document.getElementById("plotMode").value;
    const grid = computeGridAndDensity(N);
    const axis = document.getElementById("sliceAxis").value;
    const idx = Math.floor(N * parseInt(document.getElementById("sliceSlider").value) / 100);
    plotSlice(grid, axis, idx);
    if (mode === "slices") document.getElementById("plot3d").style.display = "none";
    else { document.getElementById("plot3d").style.display = "block"; plot3D(grid, mode); }
    document.getElementById("sliceSlider").oninput = function () { const idx = Math.floor(N * parseInt(this.value) / 100); plotSlice(grid, document.getElementById("sliceAxis").value, idx); };
    document.getElementById("sliceAxis").onchange = function () { const idx = Math.floor(N * parseInt(document.getElementById("sliceSlider").value) / 100); plotSlice(grid, this.value, idx); };
});

// --- Particle Background ---
const canvas = document.getElementById("particles-bg"); const ctx = canvas.getContext("2d"); let particles = [];
function resizeCanvas() { canvas.width = window.innerWidth; canvas.height = window.innerHeight; }
window.addEventListener("resize", resizeCanvas); resizeCanvas();
for (let i = 0; i < 120; i++) { particles.push({ x: Math.random() * canvas.width, y: Math.random() * canvas.height, r: Math.random() * 2, dx: (Math.random() - 0.5) * 0.5, dy: (Math.random() - 0.5) * 0.5 }); }
function animateParticles() { ctx.clearRect(0, 0, canvas.width, canvas.height); particles.forEach(p => { ctx.beginPath(); ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2); ctx.fillStyle = "cyan"; ctx.shadowBlur = 20; ctx.shadowColor = "cyan"; ctx.fill(); p.x += p.dx; p.y += p.dy; if (p.x < 0 || p.x > canvas.width) p.dx *= -1; if (p.y < 0 || p.y > canvas.height) p.dy *= -1; }); requestAnimationFrame(animateParticles); }
animateParticles();
