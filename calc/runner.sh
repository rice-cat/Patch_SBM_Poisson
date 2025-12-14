#!/usr/bin/env bash
# runner.sh - simplified runner adapted for Patch_SBM_Poisson
# Usage: ./runner.sh [SETTINGS_FILE]
set -euo pipefail

# Default settings file
SETTINGS_FILE="${1:-settings}"

if [ ! -f "$SETTINGS_FILE" ]; then
    cat >&2 <<EOF
Settings file not found: $SETTINGS_FILE

Create a settings file with the following (example `calc/settings`):
RESULTS=results_2d
MODE=2d
PREFIX=/path/to/executables
TEMPLATE=mg_param.prm
SHYNESS=3
FEDEGREE=1,2,3
REFINEMENTS=4,5

Usage: $0 [SETTINGS_FILE]
EOF
    exit 1
fi

# Read settings file, ignoring comments and empty lines
read_settings() {
    while IFS= read -r line || [ -n "$line" ]; do
        line=${line%$'\r'}
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        [[ -z "$line" || "${line:0:1}" = "#" ]] && continue
        key=${line%%=*}
        value=${line#*=}
        key="${key#"${key%%[![:space:]]*}"}"
        key="${key%"${key##*[![:space:]]}"}"
        value="${value#"${value%%[![:space:]]*}"}"
        value="${value%"${value##*[![:space:]]}"}"
        case "$key" in
            MODE|PREFIX|TEMPLATE|RESULTS|SHYNESS|FEDEGREE|REFINEMENTS|N_JOBS)
                declare -g "$key"="$value"
                ;;
        esac
    done < "$SETTINGS_FILE"
}

read_settings

# Defaults
RESULTS="${RESULTS:-results}"

# Validate required parameters
if [[ -z "${MODE:-}" || -z "${PREFIX:-}" || -z "${TEMPLATE:-}" || -z "${FEDEGREE:-}" || -z "${REFINEMENTS:-}" ]]; then
    echo "Missing required parameters in $SETTINGS_FILE" >&2
    echo "Required: MODE, PREFIX, TEMPLATE, FEDEGREE, REFINEMENTS" >&2
    exit 2
fi

case "$MODE" in
    2d) EXE="sbm_mg_solver_2d" ;;
    3d) EXE="sbm_mg_solver_3d" ;;
    *) echo "MODE must be '2d' or '3d'"; exit 2 ;;
esac

EXE_PATH="${PREFIX%/}/$EXE"
if [ ! -x "$EXE_PATH" ]; then
    echo "Executable not found or not executable: $EXE_PATH" >&2
    exit 3
fi
if [ ! -f "$TEMPLATE" ]; then
    echo "Template parameter file not found: $TEMPLATE" >&2
    exit 4
fi

echo "Running with settings from: $SETTINGS_FILE"
echo "MODE=$MODE, PREFIX=$PREFIX, TEMPLATE=$TEMPLATE"
echo "FEDEGREE=$FEDEGREE, REFINEMENTS=$REFINEMENTS, SHYNESS=${SHYNESS:-(unset)}"

IFS=',' read -r -a fedegs <<< "$FEDEGREE"
IFS=',' read -r -a refinements <<< "$REFINEMENTS"
IFS=',' read -r -a shy_list <<< "${SHYNESS:-}"

mkdir -p "${RESULTS%/}"

for fedeg in "${fedegs[@]}"; do
    fedeg=$(echo "$fedeg" | xargs)
    if ! [[ "$fedeg" =~ ^[0-9]+$ ]]; then
        echo "Invalid FEDEGREE value: $fedeg" >&2
        exit 5
    fi
    for ref in "${refinements[@]}"; do
        ref=$(echo "$ref" | xargs)
        if ! [[ "$ref" =~ ^[0-9]+$ ]]; then
            echo "Invalid REFINEMENTS value: $ref" >&2
            exit 5
        fi
        # support empty SHYNESS (use template value) or comma-separated list
        if [ -z "${SHYNESS:-}" ]; then
            shy_vals=("")
        else
            shy_vals=("${shy_list[@]}")
        fi
        for shy in "${shy_vals[@]}"; do
            shy=$(echo "$shy" | xargs)
            MODE_UPPER=$(echo "$MODE" | tr '[:lower:]' '[:upper:]')
            folder="${RESULTS%/}/${MODE_UPPER}_p${fedeg}_ref${ref}"
            [ -n "$shy" ] && folder+="_shy${shy}"
            mkdir -p "$folder"

            cp "$TEMPLATE" "$folder/$(basename "$TEMPLATE")"
            PARAM_FILE="$folder/$(basename "$TEMPLATE")"

            # Replace relevant parameters in the parameter file
            sed -E \
                -e "s/^[[:space:]]*set[[:space:]]+fe_degree[[:space:]]*=.*/  set fe_degree = ${fedeg}/" \
                -e "s/^[[:space:]]*set[[:space:]]+n_refinements[[:space:]]*=.*/  set n_refinements = ${ref}/" \
                -e "s/^[[:space:]]*set[[:space:]]+shyness[[:space:]]*=.*/  set shyness = ${shy}/" \
                "$PARAM_FILE" > "$PARAM_FILE.tmp" && mv "$PARAM_FILE.tmp" "$PARAM_FILE" || {
                echo "Failed to update parameter file: $PARAM_FILE" >&2
                exit 6
            }

            echo "Running $EXE_PATH with params $PARAM_FILE -> $folder/output.txt"
            (cd "$folder" && "$EXE_PATH" "$(basename "$PARAM_FILE")") >"$folder/output.txt" 2>&1 || {
                echo "Execution failed for $folder (see $folder/output.txt)" >&2
                exit 7
            }
        done
    done
done

echo "All runs finished."