# HaloMaker Python

RAMSES snapshot을 위한 HaloMaker / AdaptaHOP의 Python + Fortran 구현입니다.
Fortran extension은 메모리 사용이 큰 neighbor 및 structure-tree 루틴을
담당합니다. 파이프라인은 **full-box 전용**이며, 레거시 zoom-in 경로는 제거되었습니다.

이 release는 dev3 memory-optimized 구현을 기준으로 준비되었습니다. 공개
파일명은 기존 script가 development branch 이름에 의존하지 않도록 canonical
이름(`compute_adaptahop.f90`)을 유지합니다. 파이프라인은 이제 **full-box
전용**이며, 레거시 zoom-in Fortran 모듈은 제거되었습니다.

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
가리키도록 설정합니다. 각 line은 `'<snapshot dir>' <format> <worker count>
<output number> [prefix]` 형식이며 5번째 `prefix`는 선택입니다. line은 순서대로
처리되므로 여러 snapshot(다른 시뮬레이션이라도)을 한 run에서 체이닝할 수 있습니다.
파이프라인은 periodic full-box 모드로 동작하며, 레거시 zoom-in 모드는 제거되었으므로
`input_HaloMaker.dat`에 `zoomin` 키가 남아 있어도 무시됩니다.

RAMSES snapshot의 경우 `lbox`, `omega_f`, `lambda_f`, `H_f`는 `read_data()` 중
RAMSES header/`info_*.txt`에서 **snapshot별로** 읽으므로 `input_HaloMaker.dat`에서
optional입니다. `H_f`는 config 값이 snapshot의 `H0`와 `rtol=1e-6` 내로 일치하면
config 값을, 아니면 snapshot의 authoritative `H0`를 사용합니다(출력 header에
`H_f_source`/`info_H0`로 기록) — 덕분에 한 run에서 `H0`가 다른 시뮬레이션도
체이닝할 수 있습니다. config 값이 생략되면 box-size·cosmology 의존 quantity는
snapshot header를 읽은 뒤 확정됩니다.

체이닝한 snapshot들이 output number를 공유하면 출력 파일명이 자동으로
구분됩니다(`tree_bricks{n}_{tag}.h5`) — 덮어쓰기 없음. 단일 snapshot은 기존
`tree_bricks{n:05d}.h5` 이름을 유지합니다.

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

### 병렬성 (`nbPes`)

`inputfiles_HaloMaker.dat` 각 라인의 **3번째 필드가 `nbPes`** 로, Python 단계(입자 읽기·density·
halo 속성)와 Fortran OpenMP 단계(이웃 탐색·mean density·saddle 연결·node tree) 양쪽의 worker
프로세스/스레드 수를 정합니다. `nbPes`를 키우면 대부분 단계가 빨라지지만 전체 speedup은 **포화**합니다
(1900만 입자 박스에서 32코어 기준 ≈12×): 구조 트리의 `create nodes` 단계가 현재 **직렬**이라
`nbPes`에 비례하지 않으며, 대형 박스에서는 이 단계가 wall time을 지배합니다.

### 출력 단위 (breaking change: `halomaker_units_v2`)

이제 catalog와 GasMaker 출력은 단일 단위계를 씁니다: **질량 `Msun`**, **위치/반경/
형상축 RAMSES code unit `[0,1)`**, 각운동량 `Msun Mpc km/s`, 에너지 `Msun (km/s)^2`,
`rho_0` `Msun/kpc^3`, 속도 `km/s`, 나이 `Gyr`. 기존 `10^11 Msun` 질량 스케일과
physical Mpc 위치/반경은 HDF5 출력에서 **더 이상 쓰지 않습니다.** 파일에는
`/header.units_version = "halomaker_units_v2"`, `box_comoving_mpc`/`box_physical_mpc`
attr, 필드별 `field_units` JSON attr이 붙습니다. physical Mpc 복원은
`x_phys = x_code * box_physical_mpc`. 옛 단위를 가정한 분석 스크립트는 `units_version`
으로 분기해야 합니다. 전체 필드 단위표는 **`CATALOG_FORMAT.md`** 참고.

per-element 항성 화학조성을 저장하는 Ra4 스냅샷의 경우, catalog에 질량가중 항성 원소비
`H_star, O_star, Fe_star, Mg_star, C_star, N_star, Si_star, S_star, D_star`(mass fraction)도
포함됩니다(항성 chem 없는 스냅샷은 `NaN`). 전 멤버 pos/vel/mass를 함께 내보내려면
`dump_members`를 켭니다(전 멤버 dump 의미; 옛 키 `dump_DMs`/`dump_stars`도 alias로 동작).

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

## GasMaker (가스 후처리기)

GasMaker는 기존 HaloMaker/GalaxyMaker catalog와 RAMSES AMR/hydro 데이터를 읽어
halo·galaxy별 **가스 물리량**을 추가하는 별도 도구입니다: total / cold(T<10⁴ K)
/ dense 가스 질량, 가스 금속·원소별 화학, 가스 운동학·각운동량(r\*, r50, r90,
r_vir 내부), 그리고 구형 overdensity(r200/m200, r500/m500 및 내부 DM/star/gas
질량). 재시작(restart)을 지원합니다.

```bash
# 특정 root만:
python GasMaker.py <catalog.h5> <ramses_repo> <iout> --root-ids 3,11,15
# 또는 전체 top-level halo:
python GasMaker.py <catalog.h5> <ramses_repo> <iout> --roots all
```

출력 기본 이름은 `gas_bricks{iout:05d}.h5`입니다(`--output`으로 변경 가능). 카탈로그와
row 정렬되어 `id`로 join되며, 두 출력의 전체 필드 목록은 **`CATALOG_FORMAT.md`** 참고.

스냅샷 reader는 **교체 가능(pluggable)** 합니다. 기본 reader
(`gasmaker/readers/rur.py`)는 `rur` 패키지를 쓰지만 **lazy import** 되므로
GasMaker core는 `rur`에 의존하지 않습니다(없어도 설치·import 가능). 기본적으로 설치된
`rur`(또는 `$RUR_PATH`)를 쓰며, checkout은 `--rur-path`로 지정하거나, 다른 시뮬/포맷은
`gasmaker.readers.base.CellReader` 인터페이스(`read_boxes`·`hydro_fields` 포함)를
구현하면 됩니다.

> 상태: stratified NH2 표본에서 RUR reference와 머신정밀로 검증됨(가스/입자 질량,
> 금속, 화학 — `WHATS_NEW.md` 참고). `r200/r500`은 threshold-crossing 보간을
> 쓰며 RUR의 nearest-shell 선택과는 의도적으로 다릅니다.

## Files

- `HaloMaker.py`: command-line entry point
- `compute_halo_props.py`: HaloMaker workflow 및 halo properties
- `input_output.py`: RAMSES reader 및 HDF5 catalog writer
- `halo_defs.py`: shared runtime state 및 utility
- `num_rec.py`: numerical helper
- `compute_neiKDtree_mod.py`: Python-to-Fortran bridge
- `compute_adaptahop.f90`: optimized full-box AdaptaHOP extension
- `compute_adaptahop.pyf`: portable build를 위한 explicit f2py interface
- `hdf_output_example.py`: 간단한 HDF5 catalog reader example
- `ssp_photometry.py`: compact SSP table interpolation
- `halomaker_data/ssp_tables`: build 시 생성되는 SSP runtime table; Git에서 제외
- `clean_runtime.sh`: interrupted run을 위한 dry-run / cleanup helper
- `GasMaker.py`: 가스 후처리기 command-line entry point
- `gasmaker/`: GasMaker 패키지 (`pipeline`, `catalog`, `geometry`, `overlap`)
- `gasmaker/readers/`: 교체 가능한 스냅샷 reader (`base` 인터페이스 + `rur` 어댑터)

## SSP table 준비

BC03, CB07, FSPS compact table은 HaloMaker와 함께 재배포하지 않습니다.
첫 build 전에 원본 model data 경로를 지정해야 합니다. Local development
copy가 있는 경우 표준 위치는 Git에서 제외되는 `assets/ssp_originals/bc03` 및
`assets/ssp_originals/cb07`입니다. 명시적 경로가 이 기본값보다 우선합니다.

- `BC03_PATH`: BC03 Chabrier/Padova 1994 source tarball 또는 extracted directory
- `CB07_PATH`: RUR에서 사용한 CB07 source-table directory
- `FSPS_PATH`: FSPS source/data installation (`SPS_HOME`도 허용)

BC03 원본은 Bruzual & Charlot 2003 original release page에서 받을 수 있습니다.
`https://www.bruzual.org/bc03/Original_version_2003/`

FSPS를 생성하려면 optional generator dependency가 필요합니다.

```bash
uv sync --extra ssp-generation
BC03_PATH=/path/to/bc03 \
CB07_PATH=/path/to/cb07 \
FSPS_PATH=/path/to/fsps \
PYTHON=.venv/bin/python bash build.sh
```

`build.sh`는 없는 table만 `halomaker_data/ssp_tables/` 아래에 생성하고,
이후 build에서는 기존 파일을 재사용한 뒤 Fortran extension을 compile합니다.
SSP table 준비 로직은 `tools/prepare_ssp_tables.sh`에 분리되어 있습니다.
생성된 npz 파일은 Git에서 제외됩니다. SSP table만 준비하거나 다시 만들려면
다음처럼 실행합니다.

```bash
PYTHON=.venv/bin/python bash tools/prepare_ssp_tables.sh
HALOMAKER_TABLES_ONLY=1 PYTHON=.venv/bin/python bash build.sh
```

개별 table을 명시적으로 다시 만들려면 다음처럼 실행합니다.

```bash
uv run python tools/generate_bc03_table.py \
    --bc03-path /path/to/bc03 --force
uv run python tools/generate_cb07_table.py \
    --cb07-path /path/to/cb07 --force
uv run python tools/generate_fsps_table.py \
    --fsps-path /path/to/fsps --force
```

## Release Checklist

Public release 전에 license를 선택하고 `LICENSE` file을 추가해야 합니다. 또한
redistribution rights가 허용한다면 작은 redistributable RAMSES fixture와
automated smoke test를 추가하는 것이 좋습니다.
