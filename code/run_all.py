"""Run the whole PPL-115 numerics campaign in parallel, memory-aware and resumable.

    python3 run_all.py            # everything
    python3 run_all.py --dry      # show the plan and resource detection, run nothing
Results land in ./results ; a summary is written to results/report.txt and packed into results.tar.gz
"""
import os, sys, time, subprocess, json, glob, tarfile

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SVN_OUT", os.path.join(ROOT, "..", "results", "rerun"))
LOG = os.path.join(OUT, "run.log")


def cpus():
    n = len(os.sched_getaffinity(0))
    try:
        q, p = open("/sys/fs/cgroup/cpu.max").read().split()
        if q != "max":
            n = min(n, max(1, int(int(q) / int(p))))
    except Exception:
        try:
            q = int(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
            p = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
            if q > 0:
                n = min(n, max(1, q // p))
        except Exception:
            pass
    return n


def mem_avail_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 2**20
    return 0.0


def external_jobs():
    """job.py processes already running (e.g. launched by a previous scheduler instance)"""
    out = []
    for pid in os.listdir("/proc"):
        if not pid.isdigit() or int(pid) == os.getpid():
            continue
        try:
            argv = open(f"/proc/{pid}/cmdline", "rb").read().split(b"\0")
        except Exception:
            continue
        argv = [a.decode(errors="ignore") for a in argv if a]
        if len(argv) >= 3 and argv[1].endswith("job.py"):
            out.append((int(pid), argv[2:]))
    return out


def mem_gb():
    lim = None
    for fn in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            v = open(fn).read().strip()
            if v != "max" and int(v) < 1 << 60:
                lim = int(v) / 2**30
                break
        except Exception:
            pass
    avail = None
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            avail = int(line.split()[1]) / 2**20
    return min(x for x in (lim, avail) if x is not None)


def est_study(N):            # GB, peak of one SuperLU saddle solve (calibrated at N=64,96)
    return max(1.0, CAL_A * (N / 96.0) ** CAL_P)


def est_thresh(N, m):
    ntri = 2 * N * m
    return max(1.0, 0.6 * CAL_A * (ntri / 4608.0) ** (CAL_P / 2))


CAL_A, CAL_P = 6.0, 2.6      # overwritten by calibration.json if present
try:
    c = json.load(open(os.path.join(ROOT, "calibration.json")))
    CAL_A, CAL_P = c["A"], c["P"]
except Exception:
    pass


def plan():
    jobs = []
    for N in (16, 24, 32, 48, 64, 96, 128, 192):
        for kind in ("A", "B1", "B100"):
            for mu in (100.0, 1000.0, 10000.0):
                jobs.append(("study", N, kind, mu))
    for kind in ("A", "B1"):                       # finest level: key cases only
        jobs.append(("study", 256, kind, 100.0))
    for N in (8, 12, 16, 20, 24, 32, 48, 64):    # threshold vs refinement (fixed rho)
        jobs.append(("thresh", N, 2.5, N // 4))
    for N in (16, 32):                            # threshold vs rho (enlarged outer box)
        for L in (2.5, 5.0, 10.0, 20.0, 40.0):
            m = max(N // 4, int(round(N / 4 * (L - 1) / 1.5)))
            jobs.append(("thresh", N, L, m))
    out = []
    for j in jobs:
        if j[0] == "study":
            _, N, kind, mu = j
            fn = os.path.join(OUT, "study", f"N{N}_{kind}_mu{mu:g}.jsonl")
            out.append(dict(key=f"study N={N} {kind} mu={mu:g}", cmd=["study", str(N), kind, str(mu)],
                            fn=fn, mem=est_study(N)))
        else:
            _, N, L, m = j
            fn = os.path.join(OUT, "thresh", f"N{N}_L{L:g}_m{m}.json")
            out.append(dict(key=f"thresh N={N} L={L:g} m={m}", cmd=["thresh", str(N), str(L), str(m)],
                            fn=fn, mem=est_thresh(N, m)))
    return out


def log(msg):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    ncpu = cpus()
    ext = external_jobs()
    ext_keys = {tuple(a) for _, a in ext}
    allp = [j for j in plan() if not os.path.exists(j["fn"])]
    ext_jobs = [j for j in allp if tuple(j["cmd"]) in ext_keys]
    mtot = mem_gb() + sum(j["mem"] for j in ext_jobs) if not os.path.exists("/sys/fs/cgroup/memory.max") else mem_gb()
    try:
        v = open("/sys/fs/cgroup/memory.max").read().strip()
        if v != "max":
            mtot = int(v) / 2**30
    except Exception:
        pass
    budget = 0.85 * mtot
    jobs = [j for j in allp if tuple(j["cmd"]) not in ext_keys]
    log(f"detected {ncpu} usable CPUs, {mtot:.1f} GB available; memory budget {budget:.1f} GB; "
        f"{len(jobs)} jobs to run (A={CAL_A}, P={CAL_P})")
    too_big = [j for j in jobs if j["mem"] > budget]
    for j in too_big:
        log(f"SKIP (est {j['mem']:.0f} GB > budget): {j['key']}")
    jobs = sorted([j for j in jobs if j["mem"] <= budget], key=lambda j: -j["mem"])
    if "--dry" in sys.argv:
        for j in jobs:
            print(f"  {j['mem']:7.1f} GB  {j['key']}")
        return
    running, retried, failed = [], set(), []
    class _Ext:
        def __init__(self, pid): self.pid = pid
        def poll(self): return None if os.path.exists(f"/proc/{self.pid}") else 0
    for pid, a in ext:
        for j in ext_jobs:
            if tuple(j["cmd"]) == tuple(a):
                running.append(dict(p=_Ext(pid), job=j, t=time.time()))
                log(f"adopt running {j['key']} (pid {pid}, reserve {j['mem']:.1f} GB)")
    t0 = time.time()
    while jobs or running:
        # reap
        for r in running[:]:
            rc = r["p"].poll()
            if rc is None:
                continue
            running.remove(r)
            j = r["job"]
            if rc == 0 and os.path.exists(j["fn"]):
                log(f"done  {j['key']}  ({time.time() - r['t']:.0f}s)")
            elif j["key"] not in retried:
                retried.add(j["key"])
                j["mem"] = min(budget, 2 * j["mem"])
                jobs.append(j); jobs.sort(key=lambda x: -x["mem"])
                log(f"FAIL rc={rc} {j['key']} -> retry with {j['mem']:.0f} GB reserved (see results/err)")
            else:
                failed.append(j["key"])
                log(f"FAIL rc={rc} {j['key']} (gave up)")
        # launch
        used = sum(r["job"]["mem"] for r in running)
        for j in jobs[:]:
            if len(running) >= ncpu:
                break
            if (used + j["mem"] <= budget and mem_avail_gb() > j["mem"] + 4.0) or not running:
                os.makedirs(os.path.join(OUT, "err"), exist_ok=True)
                errf = open(os.path.join(OUT, "err", j["key"].replace(" ", "_").replace("=", "") + ".txt"), "w")
                p = subprocess.Popen([sys.executable, os.path.join(ROOT, "job.py")] + j["cmd"],
                                     stdout=errf, stderr=subprocess.STDOUT, cwd=ROOT)
                running.append(dict(p=p, job=j, t=time.time()))
                jobs.remove(j); used += j["mem"]
                log(f"start {j['key']}  (reserve {j['mem']:.1f} GB; running {len(running)}, queued {len(jobs)})")
        time.sleep(3)
    log(f"all jobs finished in {(time.time() - t0) / 60:.1f} min; failures: {failed or 'none'}")
    rep = subprocess.run([sys.executable, os.path.join(ROOT, "report.py"), OUT], capture_output=True, text=True, cwd=ROOT)
    open(os.path.join(OUT, "report.txt"), "w").write(rep.stdout + rep.stderr)
    with tarfile.open(os.path.join(OUT, "..", "results.tar.gz"), "w:gz") as tf:
        tf.add(OUT, arcname="results")
    log("wrote results/report.txt and results.tar.gz")


if __name__ == "__main__":
    main()
