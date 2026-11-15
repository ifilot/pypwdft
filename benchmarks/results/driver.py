"""Budgeted, escalating benchmark sweep. Each case runs in its own process
with a hard timeout; a series stops as soon as one case gets expensive."""
import json, subprocess, sys, time, os
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]
BUDGET = float(sys.argv[2]) * 60
START = time.time()
TIMEOUT = {"cpu": 420, "cuda": 240}

def log(msg):
    print(f"[{(time.time()-START)/60:6.1f} min] {msg}", flush=True)

def run(stage, **cfg):
    cfg.setdefault("molecule", "CO"); cfg.setdefault("cell", 10.0)
    cfg.setdefault("xc", "pbe"); cfg.setdefault("eigensolver", "davidson")
    if cfg["device"] == "cpu":
        cfg.setdefault("threads", 8)
    limit = TIMEOUT[cfg["device"]]
    left = BUDGET - (time.time() - START)
    if left < limit:
        log(f"{stage}: budget exhausted, skipping {cfg}")
        return None
    rec = dict(cfg, stage=stage)
    t = time.time()
    try:
        p = subprocess.run([sys.executable, os.path.join(HERE, "worker.py"),
                            json.dumps(cfg)], capture_output=True, text=True,
                           timeout=limit)
        line = [l for l in p.stdout.splitlines() if l.startswith("RESULT ")]
        if line:
            rec.update(json.loads(line[-1][7:])); rec["status"] = "ok"
        else:
            rec["status"] = "error"; rec["error"] = p.stderr.strip().splitlines()[-1:]
    except subprocess.TimeoutExpired:
        rec["status"] = "timeout"; rec["timeout"] = limit
    rec["wall"] = time.time() - t
    with open(OUT, "a") as f:
        f.write(json.dumps(rec) + "\n")
    tag = {k: cfg[k] for k in ("molecule", "cell", "cutoff", "device", "threads",
                               "eigensolver", "fft_backend", "bands", "xc") if cfg.get(k) is not None}
    if rec["status"] == "ok":
        log(f"{stage}: {tag} -> best {rec['best']:.2f} s "
            f"(grid {rec['grid']}^3, {rec['plane_waves']} PW, {rec['iterations']} it, "
            f"E={rec['energy']:.8f})")
    else:
        log(f"{stage}: {tag} -> {rec['status']} {rec.get('error','')}")
    return rec

def series(stage, cases, stop_after, devices=("cuda", "cpu"), **common):
    """Run cases in ascending cost; per device, stop once a case exceeds
    stop_after seconds (or fails), so the next, larger one is never started."""
    for dev in devices:
        for case in cases:
            rec = run(stage, device=dev, **common, **case)
            if rec is None:
                return
            if rec["status"] != "ok" or rec["best"] > stop_after[dev]:
                log(f"{stage}: stopping {dev} series here")
                break

STOP = {"cuda": 30, "cpu": 60}
# 1. cutoff scaling
series("cutoff", [dict(cutoff=c) for c in (10, 15, 20, 30, 40, 50, 60, 80, 100)], STOP)
# 2. CPU thread scaling
for th in (1, 2, 4, 8):
    run("threads", device="cpu", cutoff=40, threads=th)
# 3. cell-size scaling
series("cell", [dict(cell=c) for c in (8, 10, 12, 14, 16, 20)], STOP, cutoff=30)
# 4. molecules (ascending electron count)
mols = ["h2", "he", "lih", "bh3", "ch4", "nh3", "h2o", "co", "bf3", "co2", "ethylene"]
series("molecule", [dict(molecule=m) for m in mols], STOP, cutoff=30, cell=12)
# 5. eigensolver comparison (ARPACK is the expensive one; small cases only)
for cutoff in (10, 20):
    for dev in ("cuda", "cpu"):
        for es in ("davidson", "arpack"):
            rec = run("eigensolver", device=dev, cutoff=cutoff, eigensolver=es)
# 6. CPU FFT backends
for be, th in (("pyfftw", 8), ("scipy", 8), ("numpy", 1)):
    run("backend", device="cpu", cutoff=30, fft_backend=be, threads=th)
# 7. functional
for xc in ("lda", "pbe"):
    for dev in ("cuda", "cpu"):
        run("xc", device=dev, cutoff=40, xc=xc)
# 8. extra (unoccupied) bands
series("bands", [dict(bands=b) for b in (5, 10, 20, 40)], STOP, cutoff=30)
# 9. benzene, the largest bundled system
series("benzene", [dict(cutoff=c) for c in (10, 20, 30, 40)],
       {"cuda": 60, "cpu": 100}, molecule="benzene", cell=15)
log("sweep finished")
