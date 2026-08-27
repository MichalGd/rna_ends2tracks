# Recovery and troubleshooting

## See what the workflow is doing

The authoritative human-readable progress log is `OUTPUT_DIR/rna_ends2tracks.log`:

```bash
tail -F /path/to/OUTPUT_DIR/rna_ends2tracks.log
rna-ends2tracks status /path/to/project/config/config.conf
```

The status command can also receive `OUTPUT_DIR` directly. It reports the workflow PID/process state, ordered stage states, free disk space, and counts of expected principal outputs. Add `--json` for the complete machine-readable snapshot. The master log reports stages and native-command outcomes in chronological order and points to the detailed log for each native tool. `logs/events.jsonl` remains available for programmatic parsing.

If the status is `failed`, the last message identifies the failed stage. Read the referenced detailed log before using `--force-step`; a normal rerun is sufficient when valid receipts already exist.

## Safe restart

Re-run the same command. Matching module, sample and contrast receipts skip complete outputs. Use `--from-step STEP` to avoid scanning earlier stages or `--force-step STEP` after deliberately changing/rebuilding that stage. A changed sample set must use a new output directory because it defines a new condition-blind PAS universe.

## Run lock

`.checkpoints/workflow.lock` records PID, host and start time. If a run is active, do not remove it. If the recorded process has ended and no workflow process targets that output directory, remove only that exact lock file and resume.

## Resource failure

Preflight fails if any pool’s parallel jobs multiplied by per-job CPU/RAM exceeds `MAX_TOTAL_THREADS` or `MAX_TOTAL_MEMORY_GB`. Reduce parallel jobs first; reduce per-tool threads only when appropriate. Check `00_metadata/resource_plan.tsv` and `.checkpoints/timings/`.

## Orientation failure

A low reverse-compatible fraction usually indicates the wrong library protocol, primer/orientation, sample identity or reference. Do not suppress it blindly. Review STAR counts and Lexogen library details; only lower `ORIENTATION_MIN_FRACTION` for a documented low-information exception.

## STAR overhang warning

STAR indices are not intrinsically tied to exact read length. A 150-overhang index can usually align 101-nt reads; the warning requests review rather than rebuilding one index per read length. FASTA/GTF/contig compatibility matters more.

## Cleanup

Cleanup requires successful receipts for all enabled final modules and the report. If interrupted before cleanup, resume through `cleanup`. `provenance/cleanup/cleanup_manifest.tsv` is cumulative. Final outputs are not cleanup targets. Retained intermediates can be controlled in `config.conf` for diagnostic runs.

## Failed environment installation

The stable launcher is changed only after tests. A failed new environment therefore leaves the prior stable release intact. Inspect the deployment log, correct the dependency/tag issue, and install a new versioned target; do not mutate a frozen production environment.
