# Top-level targets.  All commands run in the foreground.
.PHONY: pdf report figures smoke certificates lemma-tests campaign local clean

pdf:
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

# Regenerate the tables from the shipped raw results (seconds)
report:
	cd code && python3 report.py ../results/pod > ../results/pod/report.regen.txt && diff -q ../results/pod/report.regen.txt ../results/pod/report.txt && echo "report.txt reproduced exactly"
	cd code && python3 report_local.py > ../results/local/report_local.regen.txt && diff -q ../results/local/report_local.regen.txt ../results/local/report_local.txt && echo "report_local.txt reproduced exactly"
	cp results/v3/report_v3.txt results/v3/report_v3.shipped.txt && cd code && python3 report_v3.py ../results/v3 > /dev/null && diff -q ../results/v3/report_v3.txt ../results/v3/report_v3.shipped.txt && echo "report_v3.txt reproduced exactly"

figures:
	cd code && python3 plots.py

# < 1 min: one GS/CNS* solve pair at N=16 for tests B1 and P, printed to stdout
smoke:
	cd code && SVN_OUT=/tmp/svn_smoke python3 job.py study 16 B1 100 && SVN_OUT=/tmp/svn_smoke python3 job.py study 16 P1 100 && cat /tmp/svn_smoke/study/*.jsonl

# All logs: exact certificates, float lambda*, lemma dual-norm tests, and the v3 logs (about an hour)
certificates lemma-tests:
	cd code && ./regen_logs.sh

# Full A/B campaign + thresholds (74 + 16 jobs, ~36 min on 32 cores / 125 GB; writes results/rerun)
campaign:
	cd code && python3 run_all.py

# Tests P and C up to N=96 (2 cores, 7 GB is enough; writes results/rerun)
local:
	cd code && for N in 16 24 32 48 64 96; do for MU in 100 1000 10000; do python3 job.py study $$N P1 $$MU; done; python3 job.py study $$N C 100; done

clean:
	cd paper && latexmk -c
