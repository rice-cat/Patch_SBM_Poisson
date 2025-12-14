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
            MODE|PREFIX|TEMPLATE|RESULTS|SHYNESS|FEDEGREE|REFINEMENTS|SMOOTHING|N_JOBS)
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
echo "FEDEGREE=$FEDEGREE, REFINEMENTS=$REFINEMENTS, SHYNESS=${SHYNESS:-(unset)}, SMOOTHING=${SMOOTHING:-(unset)}"

IFS=',' read -r -a fedegs <<< "$FEDEGREE"
IFS=',' read -r -a refinements <<< "$REFINEMENTS"
IFS=',' read -r -a shy_list <<< "${SHYNESS:-}"
IFS=',' read -r -a smooth_list <<< "${SMOOTHING:-}"

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
            if [ -z "${SMOOTHING:-}" ]; then
                smooth_vals=("")
            else
                smooth_vals=("${smooth_list[@]}")
            fi
            for smooth in "${smooth_vals[@]}"; do
                smooth=$(echo "$smooth" | xargs)
                if [ -n "$smooth" ] && ! [[ "$smooth" =~ ^[0-9]+$ ]]; then
                    echo "Invalid SMOOTHING value: $smooth" >&2
                    exit 5
                fi
            MODE_UPPER=$(echo "$MODE" | tr '[:lower:]' '[:upper:]')
            folder="${RESULTS%/}/${MODE_UPPER}_p${fedeg}_ref${ref}"
            [ -n "$shy" ] && folder+="_shy${shy}"
                [ -n "$smooth" ] && folder+="_smo${smooth}"
            mkdir -p "$folder"

            cp "$TEMPLATE" "$folder/$(basename "$TEMPLATE")"
            PARAM_FILE="$folder/$(basename "$TEMPLATE")"

            # Replace relevant parameters in the parameter file using awk (handles optional smoothing)
            awk -v fe="${fedeg}" -v refv="${ref}" -v shyv="${shy}" -v smo="${smooth}" '
            BEGIN{
                re_fe="^[[:space:]]*set[[:space:]]+fe_degree[[:space:]]*=.*$";
                re_ref="^[[:space:]]*set[[:space:]]+n_refinements[[:space:]]*=.*$";
                re_shy="^[[:space:]]*set[[:space:]]+shyness[[:space:]]*=.*$";
                re_smo="^[[:space:]]*set[[:space:]]+n_smoothing_steps[[:space:]]*=.*$";
                done_fe=0; done_ref=0; done_shy=0; done_smo=0;
            }
            {
                if(!done_fe && $0 ~ re_fe){ sub("=.*$","= " fe); done_fe=1 }
                if(!done_ref && $0 ~ re_ref){ sub("=.*$","= " refv); done_ref=1 }
                if(!done_shy && $0 ~ re_shy && shyv!=""){ sub("=.*$","= " shyv); done_shy=1 }
                if(!done_smo && $0 ~ re_smo && smo!=""){ sub("=.*$","= " smo); done_smo=1 }
                print
            }
            ' "$PARAM_FILE" > "$PARAM_FILE.tmp" && mv "$PARAM_FILE.tmp" "$PARAM_FILE" || {
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
done

echo "All runs finished."