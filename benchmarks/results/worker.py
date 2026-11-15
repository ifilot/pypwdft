"""Run one benchmark case and print a single JSON line with the result."""
import json, sys, time
cfg = json.loads(sys.argv[1])
from pypwdft import PWDFT, Structure

t0 = time.perf_counter()
calc = PWDFT(Structure.from_name(cfg["molecule"], cell=cfg["cell"]),
             cutoff=cfg["cutoff"], xc=cfg.get("xc", "pbe"),
             device=cfg["device"], fft_backend=cfg.get("fft_backend"),
             threads=cfg.get("threads"))
setup = time.perf_counter() - t0
gpu = cfg["device"] == "cuda"
times = []
# GPU: first run is cold (start-up + kernel compilation) and reported separately.
nruns = cfg.get("repeats", 2) + (1 if gpu else 0)
for i in range(nruns):
    t = time.perf_counter()
    r = calc.run(convergence=cfg.get("convergence", 1e-6), bands=cfg.get("bands"),
                 eigensolver=cfg.get("eigensolver", "davidson"))
    times.append(time.perf_counter() - t)
    if times[-1] > cfg.get("no_repeat_above", 20):
        break
out = dict(cfg, setup=setup, times=times, iterations=r.iterations,
           h_applications=r.scf.hamiltonian_applications,
           energy=r.energy.total, grid=r.basis.density_grid,
           plane_waves=r.basis.plane_waves, orbitals=len(r.orbitals.energies))
if gpu:
    import cupy
    out["gpu_pool_mb"] = cupy.get_default_memory_pool().total_bytes() / 2**20
    out["cold"] = times[0]
    out["warm"] = min(times[1:]) if len(times) > 1 else None
    out["best"] = out["warm"] if out["warm"] is not None else times[0]
else:
    out["best"] = min(times)
print("RESULT " + json.dumps(out), flush=True)
