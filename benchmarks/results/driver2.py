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
            if (rec["status"] != "ok" or rec["best"] > stop_after[dev]
                    or rec.get("gpu_pool_mb", 0) > 14000):
                log(f"{stage}: stopping {dev} series here")
                break


STOP = {"cuda": 45, "cpu": 100}
# A. timing noise: six repeats of one case
for dev in ("cuda", "cpu"):
    run("variance", device=dev, cutoff=40, repeats=6)
# B. thread scaling on larger problems (single FFT thread leaves BLAS threaded)
for th in (1, 2, 4, 8):
    run("threads-large", device="cpu", cutoff=30, cell=16, threads=th)
for th in (1, 8):
    run("threads-benzene", device="cpu", molecule="benzene", cell=15, cutoff=20, threads=th)
# C. convergence threshold
for conv in (1e-4, 1e-6, 1e-8):
    for dev in ("cuda", "cpu"):
        run("convergence", device=dev, cutoff=30, convergence=conv)
# D. ARPACK against Davidson at larger sizes
for es in ("davidson", "arpack"):
    series("eigensolver-large", [dict(cutoff=c) for c in (30, 40, 60)], STOP, eigensolver=es)
# E. large GPU-only cases
G = ("cuda",)
series("gpu-cutoff", [dict(cutoff=c) for c in (150, 200, 300, 400)], STOP, devices=G)
series("gpu-cell", [dict(cell=c) for c in (24, 28, 32, 40)], STOP, devices=G, cutoff=30)
series("gpu-benzene", [dict(cutoff=c) for c in (60, 80, 100, 150)], STOP, devices=G, molecule="benzene", cell=15)
series("gpu-benzene-cell", [dict(cell=c) for c in (20, 25, 30)], STOP, devices=G, molecule="benzene", cutoff=40)
series("gpu-benzene-bands", [dict(bands=b) for b in (30, 60, 120)], STOP, devices=G, molecule="benzene", cell=15, cutoff=30)
# F. one step further on the CPU for the largest system
series("cpu-benzene", [dict(cutoff=c) for c in (60, 80)], STOP, devices=("cpu",), molecule="benzene", cell=15)
log("sweep finished")
