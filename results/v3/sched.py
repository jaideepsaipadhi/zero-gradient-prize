"""Memory-aware job scheduler for job_k.py:  python3 sched.py JOBFILE LOG
At most 2 jobs at once; a job starts only when the memory-cgroup headroom exceeds its estimated peak RSS (+0.7 GB margin),
estimated from the measured N = 64 peak scaled by (N/64)^2.  Each job runs under ulimit -v 5.5 GB.
Skips jobs whose output file already exists."""
import os, sys, time, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "..", "..", "code")
RSS64 = {("std", 5): 0.92, ("std", 6): 1.6, ("std", 4): 0.8, ("ct", 2): 0.19, ("ct", 3): 0.55, ("ct", 4): 1.0}  # GB, fitted to measured peaks


CG = "/sys/fs/cgroup/memory/process_api/01a0e678-92d4-730e-a52f-a4606e8329e7/claude-code-bash/"


def memavail():
    """headroom in the shared bash memory cgroup (limit ~5.8 GiB, shared with other agents' jobs)"""
    try:
        lim = int(open(CG + "memory.limit_in_bytes").read()); use = int(open(CG + "memory.usage_in_bytes").read())
        return (lim - use) / 1024**3
    except Exception:
        for l in open("/proc/meminfo"):
            if l.startswith("MemAvailable"):
                return int(l.split()[1]) / 1024**2


def est(job):
    a = job.split()
    k, mesh, N = int(a[1]), a[2], int(a[3])
    return max(0.3, RSS64.get((mesh, k), 1.5) * (N / 64) ** 2 * 1.1) * (1.0 if a[0] == "study" else 0.5)


def outpath(job):
    a = job.split()
    if a[0] == "study":
        return os.path.join(HERE, "study", f"{a[2]}_k{a[1]}_N{a[3]}_{a[4]}_mu{float(a[5]):g}.jsonl")
    return os.path.join(HERE, "thresh", f"{a[2]}_k{a[1]}_N{a[3]}.json")


jobs = [l.strip() for l in open(sys.argv[1]) if l.strip() and not os.path.exists(outpath(l.strip()))]
log = open(sys.argv[2], "a")
env = dict(os.environ, PYPARDISO_MKL_RT="/usr/local/lib/libmkl_rt.so.3")
running = []
while jobs or running:
    for (p, j, t0) in running[:]:
        if p.poll() is not None:
            running.remove((p, j, t0))
            log.write(f"{time.strftime('%H:%M:%S')} rc={p.returncode} {time.time() - t0:.0f}s {j}\n"); log.flush()
    if jobs and len(running) < 2:
        j = jobs[0]
        if memavail() > est(j) + 0.7 or (not running and memavail() > est(j) + 0.2):
            jobs.pop(0)
            p = subprocess.Popen(["bash", "-c", f"ulimit -v 5500000; exec python3 job_k.py {j}"], cwd=CODE, env=env,
                                 stdout=subprocess.DEVNULL, stderr=open(os.path.join(HERE, "err.log"), "a"))
            running.append((p, j, time.time()))
            time.sleep(20)            # let it reach its footprint before judging memory again
            continue
    time.sleep(5)
