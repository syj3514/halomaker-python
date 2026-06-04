# Patch Branches

Patch branches are optional developer/testing overrides. They let you build or
run an experimental implementation without editing the canonical release
sources in the repository root.

The default release path is:

```bash
bash build.sh
bash run.sh
```

To test a local patch branch, create a directory under `patch/` and place only
the files you want to override:

```text
patch/my_experiment/
  compute_adaptahop.f90
  compute_adaptahop_zoomin.f90
```

Then build or run with:

```bash
BRANCH=my_experiment bash build.sh
BRANCH=my_experiment bash run.sh
```

Files that are absent from the patch branch fall back to the repository root.
If a public Fortran signature changes, put the matching `.pyf` file in the
same patch branch.

Local patch directories are ignored by git by default. Commit a patch branch
only when it is intentionally part of the public release.
