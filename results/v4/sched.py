"""Memory-aware job scheduler for code/job_v4.py:   python3 sched.py JOBFILE LOG
At most 2 jobs at once.  A job starts only when the headroom of the shared bash memory cgroup exceeds
its estimated peak RSS (+0.6 GB margin); a job estimated above 2.4 GB never runs next to another one.
Estimates: 0.26 GB per 1000 triangles (+0.15 GB) for Stokes (measured: ellipse N=192, 13402 triangles,
3.40 GB), x1.25 for Navier-Stokes.  Each job runs under ulimit -v 6 GB.  Skips jobs whose output exists."""
import os, sys, time, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "..", "..", "code")
# triangle counts of the mesh families (mesh_general, fixed parameters in job_v4.py)
NTRI = {"twodisk": {16: 266, 24: 450, 32: 734, 48: 1532, 64: 2552, 96: 5436, 128: 9482},
        "ellipse": {16: 98, 24: 180, 32: 288, 48: 544, 64: 930, 96: 1774, 128: 2988, 192: 6400, 256: 11244},
        "star": {16: 110, 24: 160, 32: 208, 48: 426, 64: 674, 96: 1270, 128: 2100, 192: 4420, 256: 7678},
        "st": {16: 240, 24: 481, 32: 726, 48: 1538, 64: 2557, 80: 4012, 96: 5726, 128: 10036}}


def _cg():
    for l in open("/proc/self/cgroup"):
        a = l.strip().split(":")
        if a[1] == "memory":
            return "/sys/fs/cgroup/memory" + a[2] + "/"


def memavail():
    try:
        cg = _cg()
        lim = int(open(cg + "memory.limit_in_bytes").read()); use = int(open(cg + "memory.usage_in_bytes").read())
        return (lim - use) / 1024**3
    except Exception:
        for l in open("/proc/meminfo"):
            if l.startswith("MemAvailable"):
                return int(l.split()[1]) / 1024**2


def est(job):
    a = job.split()
    if a[0] == "gd":
        nt = NTRI[a[1]][int(a[2])]; f = 1.0
    elif a[0] == "ns":
        nt = int(a[1]) ** 2 // 2; f = 1.25
    elif a[0] == "st":
        nt = NTRI["st"][int(a[1])]; f = 1.25
    else:
        nt = NTRI.get(a[1], {}).get(int(a[2]), 4000); f = 0.6
    return f * (0.15 + 0.26 * nt / 1000)


def outpath(job):
    a = job.split()
    if a[0] == "gd":
        return os.path.join(HERE, "gd", f"{a[1]}_N{a[2]}_{a[3]}_mu{float(a[4]):g}.jsonl")
    if a[0] == "ns":
        return os.path.join(HERE, "ns", f"circle_N{a[1]}_{a[2]}_mu{float(a[3]):g}_nu{float(a[4]):g}_{a[5]}.jsonl")
    if a[0] == "st":
        return os.path.join(HERE, "st", f"st_N{a[1]}_mu{float(a[2]):g}.jsonl")
    return os.path.join(HERE, "thresh", f"{a[1]}_N{a[2]}.json")


if __name__ == "__main__":
    jobs = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")
            and not os.path.exists(outpath(l.strip()))]
    log = open(sys.argv[2], "a")
    env = dict(os.environ, PYPARDISO_MKL_RT="/usr/local/lib/libmkl_rt.so.3")
    running = []
    log.write(f"{time.strftime('%H:%M:%S')} start {len(jobs)} jobs\n"); log.flush()
    while jobs or running:
        for (p, j, t0) in running[:]:
            if p.poll() is not None:
                running.remove((p, j, t0))
                log.write(f"{time.strftime('%H:%M:%S')} rc={p.returncode} {time.time() - t0:.0f}s {j}\n"); log.flush()
        if jobs and len(running) < 2:
            j = jobs[0]; e = est(j)
            big_running = any(est(r[1]) > 2.4 for r in running)
            ok = (not running and memavail() > e + 0.2) or \
                 (running and not big_running and e <= 2.4 and memavail() > e + 0.6)
            if ok:
                jobs.pop(0)
                p = subprocess.Popen(["bash", "-c", f"ulimit -v 6000000; exec python3 job_v4.py {j}"], cwd=CODE,
                                     env=env, stdout=subprocess.DEVNULL, stderr=open(os.path.join(HERE, "err.log"), "a"))
                running.append((p, j, time.time()))
                log.write(f"{time.strftime('%H:%M:%S')} start est={e:.2f}GB avail={memavail():.2f}GB {j}\n"); log.flush()
                time.sleep(15)
                continue
        time.sleep(5)
    log.write(f"{time.strftime('%H:%M:%S')} done\n"); log.flush()
