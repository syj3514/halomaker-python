# HaloMaker Python

RAMSES snapshot을 위한 HaloMaker / AdaptaHOP의 Python + Fortran 구현입니다.
Fortran extension은 메모리 사용이 큰 neighbor 및 structure-tree 루틴을
담당합니다. full-box와 zoom-in workflow를 모두 포함합니다.

이 release는 dev3 memory-optimized 구현을 기준으로 준비되었습니다. 공개
파일명은 기존 script가 development branch 이름에 의존하지 않도록 canonical
이름(`compute_adaptahop.f90`, `compute_adaptahop_zoomin.f90`)을 유지합니다.

개발 중 사용한 local reference 구현과 비교했을 때, dev3는 39990 full-box
test의 peak RSS를 약 104 GiB에서 약 24 GiB로 줄였고 catalog-level 결과를
보존했습니다. 테스트한 workload에서 runtime은 reference path와 비슷한
수준을 유지했습니다. experimental dev4 zoom-in read-compact path는 이
release staging copy에 포함되어 있지 않습니다.

## Requirements

- Linux
- Python 3.10 이상
- `gfortran`
- OpenMP runtime (`libgomp`)
- RAMSES snapshot files

## Python Environment

`uv`, conda, 또는 다른 Python environment manager 중 하나를 선택합니다.
Python 3.10 이상을 지원합니다. 아래 예시는 Python 3.12를 사용하지만, 지원되는
다른 version으로 바꿔도 됩니다.

### uv

원하는 Python version으로 project environment를 만듭니다.

```bash
uv venv --python 3.12
source .venv/bin/activate
uv sync
```

`uv sync`는 local `uv.lock` 파일을 생성합니다. 이 lock file은 배포에
포함하지 않습니다. 사용자가 선택한 Python version과 platform에 맞춰
dependency를 직접 resolve할 수 있도록 하기 위함입니다.

### conda

원하는 Python version으로 environment를 만들고 activate합니다.

```bash
conda create -n halomaker-python python=3.12
conda activate halomaker-python
python -m pip install -e .
```

또는 포함된 environment file을 사용해 conda가 지원되는 Python version을
resolve하게 할 수 있습니다.

```bash
conda env create -f environment.yml
conda activate halomaker-python
```

## Build

두 f2py extension을 compile합니다.

```bash
bash build.sh
```

이 명령은 다음 파일을 생성합니다.

- `compute_adaptahop*.so`
- `compute_adaptahop_zoomin*.so`

Build와 run은 같은 activated Python environment에서 수행해야 합니다.
Environment를 activate하지 않으려면 interpreter를 명시적으로 전달합니다.

```bash
PYTHON=.venv/bin/python bash build.sh
PYTHON=.venv/bin/python bash run.sh
```

## Configure

Example input을 repository root로 복사합니다.

```bash
cp examples/input_HaloMaker.dat.example input_HaloMaker.dat
cp examples/inputfiles_HaloMaker.dat.example inputfiles_HaloMaker.dat
```

`inputfiles_HaloMaker.dat`를 편집하여 각 active line이 실제 RAMSES snapshot을
가리키도록 설정합니다. Zoom-in processing에는 `input_HaloMaker.dat`에서
`zoomin = .true.`를 사용하고, periodic full-box processing에는
`zoomin = .false.`를 사용합니다.

RAMSES snapshot의 경우 `lbox`, `omega_f`, `lambda_f`는 optional입니다.
Code는 `read_data()` 중 RAMSES AMR header에서 authoritative box size와
snapshot cosmology를 읽습니다. 이 값들이 생략된 경우 box-size-dependent 및
cosmology-dependent quantity는 snapshot header를 읽은 뒤 확정됩니다.

권장 print level은 다음과 같습니다.

- `verbose = .false.`, `megaverbose = .false.`: compact production log.
- `verbose = .true.`, `megaverbose = .false.`: 주요 count와 timing을 포함한
  normal diagnostic log.
- `megaverbose = .true.`: 자세한 Fortran 및 memory tracking output을 포함한
  development log.

## Run

```bash
bash run.sh
```

주요 HDF5 catalog output은 `tree_bricks*.h5`로 기록됩니다.

Run이 Ctrl-C로 interrupt되었거나 scheduler에 의해 kill되었거나 Python
`forkserver` / `resource_tracker` process가 남은 경우, 다음 명령으로 runtime
leftover를 확인하고 정리할 수 있습니다.

```bash
bash clean_runtime.sh
bash clean_runtime.sh --force
```

첫 번째 명령은 dry run입니다. `--force` 명령은 이 repository의 HaloMaker
runtime process 중 matching되는 항목을 종료하고, 현재 user가 소유한 matching
shared memory file을 제거합니다.

## Files

- `HaloMaker.py`: command-line entry point
- `compute_halo_props.py`: HaloMaker workflow 및 halo properties
- `input_output.py`: RAMSES reader 및 HDF5 catalog writer
- `halo_defs.py`: shared runtime state 및 utility
- `num_rec.py`: numerical helper
- `compute_neiKDtree_mod.py`: Python-to-Fortran bridge
- `compute_adaptahop.f90`: optimized full-box AdaptaHOP extension
- `compute_adaptahop_zoomin.f90`: optimized zoom-in AdaptaHOP extension
- `compute_adaptahop*.pyf`: portable build를 위한 explicit f2py interface
- `hdf_output_example.py`: 간단한 HDF5 catalog reader example
- `clean_runtime.sh`: interrupted run을 위한 dry-run / cleanup helper

## Release Checklist

Public release 전에 license를 선택하고 `LICENSE` file을 추가해야 합니다. 또한
redistribution rights가 허용한다면 작은 redistributable RAMSES fixture와
automated smoke test를 추가하는 것이 좋습니다.
