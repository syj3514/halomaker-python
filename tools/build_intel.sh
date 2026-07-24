#!/usr/bin/env bash
set -euo pipefail

# Build the f2py extension without Meson, which can misidentify Intel Fortran.
PYTHON="${PYTHON:-python}"
FC="${FC:-${F90:-ifx}}"
CC="${CC:-icx}"
BRANCH="${BRANCH:-main}"
SCRIPT="$(readlink -f "${BASH_SOURCE[0]}")"
ROOT="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
MODULE="compute_adaptahop"

if [[ "$BRANCH" == "main" ]]; then
    SOURCE_DIR="$ROOT/src"
else
    SOURCE_DIR="$ROOT/patch/$BRANCH"
    if [[ ! -d "$SOURCE_DIR" ]]; then
        echo "Unknown patch branch: $BRANCH" >&2
        exit 1
    fi
fi
OUTPUT_DIR="${OUTPUT_DIR:-$SOURCE_DIR}"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

SOURCE="$SOURCE_DIR/$MODULE.f90"
INTERFACE="$SOURCE_DIR/$MODULE.pyf"
[[ -f "$SOURCE" ]] || SOURCE="$ROOT/fortran/$MODULE.f90"
[[ -f "$INTERFACE" ]] || INTERFACE="$ROOT/fortran/$MODULE.pyf"

for command in "$PYTHON" "$FC" "$CC" nm ldd; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Required command not found: $command" >&2
        exit 1
    fi
done

resolve_command() {
    local resolved
    resolved="$(command -v "$1")"
    if [[ "$resolved" == /* ]]; then
        printf '%s\n' "$resolved"
    else
        printf '%s/%s\n' "$PWD" "$resolved"
    fi
}
PYTHON="$(resolve_command "$PYTHON")"
FC="$(resolve_command "$FC")"
CC="$(resolve_command "$CC")"

[[ -f "$SOURCE" ]] || { echo "Missing source: $SOURCE" >&2; exit 1; }
[[ -f "$INTERFACE" ]] || { echo "Missing interface: $INTERFACE" >&2; exit 1; }

case "$(basename "$FC")" in
    ifort)
        visibility_flags=()
        ;;
    ifx)
        # ifx may otherwise hide the f2py module initializer at link time.
        visibility_flags=(-fvisibility=default)
        ;;
    *)
        echo "tools/build_intel.sh requires FC=ifort or FC=ifx (got $FC)." >&2
        exit 1
        ;;
esac

if [[ -n "${BUILD_DIR:-}" ]]; then
    mkdir -p "$BUILD_DIR"
    BUILD_DIR="$(cd "$BUILD_DIR" && pwd)"
    cleanup=false
else
    BUILD_DIR="$(mktemp -d "${TMPDIR:-/tmp}/halomaker-intel-build.XXXXXX")"
    cleanup=true
fi
if $cleanup; then
    trap 'rm -rf "$BUILD_DIR"' EXIT
fi

PYINC="$("$PYTHON" -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
NPINC="$("$PYTHON" -c 'import numpy; print(numpy.get_include())')"
F2PYINC="$("$PYTHON" -c \
    'import pathlib, numpy.f2py; print(pathlib.Path(numpy.f2py.__file__).parent / "src")')"
EXT_SUFFIX="$("$PYTHON" -c \
    'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
ARTIFACT="$BUILD_DIR/$MODULE$EXT_SUFFIX"
INSTALLED_ARTIFACT="$OUTPUT_DIR/$MODULE$EXT_SUFFIX"

cd "$BUILD_DIR"

# The checked-in .pyf defines the public API. Generate wrappers only, then use
# Intel compilers directly so Meson cannot drop or reinterpret OpenMP flags.
"$PYTHON" -m numpy.f2py "$INTERFACE"
[[ -f "$MODULE"module.c ]] || {
    echo "f2py did not generate ${MODULE}module.c" >&2
    exit 1
}

"$FC" -c -fPIC -O3 -qopenmp -fpp -free -heap-arrays -m64 \
    "${visibility_flags[@]}" \
    "$SOURCE" -o "$MODULE.o"

wrapper_objects=()
shopt -s nullglob
for wrapper in "$MODULE"-f2pywrappers.f "$MODULE"-f2pywrappers2.f90; do
    [[ -f "$wrapper" ]] || continue
    object="${wrapper%.*}.o"
    wrapper_flags=(-c -fPIC -O3 -qopenmp "${visibility_flags[@]}" -I"$BUILD_DIR")
    if [[ "$wrapper" == *.f90 ]]; then
        wrapper_flags+=(-free)
    fi
    "$FC" "${wrapper_flags[@]}" "$wrapper" -o "$object"
    wrapper_objects+=("$object")
done
if [[ ${#wrapper_objects[@]} -eq 0 ]]; then
    echo "f2py did not generate any Fortran wrapper source." >&2
    exit 1
fi

"$CC" -c -fPIC -O2 -std=c99 \
    -I"$PYINC" -I"$NPINC" -I"$F2PYINC" \
    "$MODULE"module.c -o "$MODULE"module.o
"$CC" -c -fPIC -O2 -std=c99 \
    -I"$PYINC" -I"$NPINC" -I"$F2PYINC" \
    "$F2PYINC/fortranobject.c" -o fortranobject.o

intel_library_flags=()
intel_libdir="${INTEL_LIBDIR:-${CMPLR_ROOT:+$CMPLR_ROOT/lib}}"
if [[ -n "$intel_libdir" ]]; then
    intel_library_flags=(-L"$intel_libdir" -Wl,-rpath,"$intel_libdir")
fi

"$FC" -shared -qopenmp "${visibility_flags[@]}" \
    "$MODULE.o" "${wrapper_objects[@]}" "$MODULE"module.o fortranobject.o \
    "${intel_library_flags[@]}" -liomp5 -o "$ARTIFACT"

"$ROOT/tools/check_openmp_extension.sh" "$ARTIFACT"
PYTHONPATH="$BUILD_DIR" "$PYTHON" -c \
    "import $MODULE as module; print('Import check passed:', module.__file__)"

mkdir -p "$OUTPUT_DIR"
install -m 755 "$ARTIFACT" "$INSTALLED_ARTIFACT"
"$ROOT/tools/check_openmp_extension.sh" "$INSTALLED_ARTIFACT"
echo "Intel build completed: $INSTALLED_ARTIFACT"
