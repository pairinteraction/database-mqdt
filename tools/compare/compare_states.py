from pathlib import Path
from typing import Literal

import pandas as pd

from utils import add_global_uid_column

COLUMNS_COMPARE_EXACT = ["id", "n", "f", "is_j_total_momentum", "is_calculated_with_mqdt", "parity"]
COLUMNS_COMPARE_NUMERICAL = [
    "energy",
    "nu",
    "exp_nui",
    "std_nui",
    "exp_l",
    "std_l",
    "exp_j",
    "std_j",
    "exp_s",
    "std_s",
    "exp_l_ryd",
    "std_l_ryd",
    "exp_j_ryd",
    "std_j_ryd",
    "underspecified_channel_contribution",
]
VERBOSE_COLUMNS: list[str] = []


def main() -> None:
    # CHANGE THESE PATHS, TO THE FOLDERS YOU WANT TO COMPARE
    db_dir = Path(__file__).parent.parent.parent / "database"
    species = "Yb174_mqdt"
    old_path = db_dir / f"{species}_v1.2"
    new_path = db_dir / f"{species}_v1.3"

    compare_states_table(new_path, old_path, compare_id=False, verbosity="none")


def compare_states_table(  # noqa: C901, PLR0912, PLR0915
    new_path: Path,
    old_path: Path,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    *,
    min_nu: float = 0,
    max_nu: float = float("inf"),
    compare_id: bool = False,
    verbosity: Literal["none", "all", "diff", "states"],
) -> None:
    """Compare the states table of two versions of the database.

    Given two parquet files containing the states table, with the following columns:
    - n: principal quantum number
    - exp_l: experimental orbital angular momentum
    - exp_j: experimental total angular momentum
    - energy: energy of the state
    - etc.
    This function will compare the tables by
    i) checking if both tables include the same states (a state is defined by the combination of n, exp_l, and exp_j)
    ii) checking if the energy values of the states are equal within a given tolerance (rtol and atol).
    iii) checking if all other columns are exactly equal.

    """
    print(f"Comparing states tables:\n  New: {new_path}\n  Old: {old_path}\n")

    species = new_path.name.split("_v")[0]
    if species != old_path.name.split("_v")[0]:
        raise ValueError(f"Cannot compare different species: {species} vs {old_path.name.split('_v')[0]}")

    paths = {"new": new_path, "old": old_path}
    states_dict = {key: pd.read_parquet(path / "states.parquet") for key, path in paths.items()}

    for key, states in states_dict.items():
        print(f"{key.capitalize()} states table:")
        shape_pre = states.shape
        states_dict[key] = states[(states["nu"] >= min_nu) & (states["nu"] <= max_nu)]
        print(f"  Table shape: Pre-filtering: {shape_pre}; After-filtering: {states.shape}")
        print()

    all_columns = COLUMNS_COMPARE_EXACT + COLUMNS_COMPARE_NUMERICAL
    for key, states in states_dict.items():
        missing_cols = [col for col in all_columns if col not in states.columns]
        if len(missing_cols) > 0:
            print(f"WARNING: {key.capitalize()} states table is missing columns: {missing_cols}")
        extra_cols = [col for col in states.columns if col not in all_columns and col != "id"]
        if len(extra_cols) > 0:
            print(f"WARNING: {key.capitalize()} states table has extra columns: {extra_cols}")

    for states in states_dict.values():
        add_global_uid_column(species, states)

    # Check if both tables have the same states
    only_in_new = set(states_dict["new"].index) - set(states_dict["old"].index)
    only_in_old = set(states_dict["old"].index) - set(states_dict["new"].index)

    if only_in_new or only_in_old:
        print("States don't match between tables:")
        if only_in_new:
            print(f"  {len(only_in_new)} states only in new table:")
            if verbosity in ["all", "states"]:
                for state in sorted(only_in_new):
                    print(f"  State({state})")
        if only_in_old:
            print(f"  {len(only_in_old)} states only in old table:")
            if verbosity in ["all", "states"]:
                for state in sorted(only_in_old):
                    print(f"  State({state})")

        # Remove non-matching states from both tables
        states_dict["new"] = states_dict["new"].drop(index=list(only_in_new), errors="ignore")
        states_dict["old"] = states_dict["old"].drop(index=list(only_in_old), errors="ignore")

    if len(states_dict["new"]) == 0:
        raise ValueError("No matching states to compare.")
    if len(states_dict["new"]) != len(states_dict["old"]):
        raise ValueError("This should not happen: After filtering, both tables should have the same number of states.")
    print(f"Continuing comparison with {len(states_dict['new'])} matching states ...\n")

    columns_without_diff: list[str] = []
    new, old = states_dict["new"], states_dict["old"]

    # Compare all columns that should be exactly equal
    for col in COLUMNS_COMPARE_EXACT:
        if not compare_id and col == "id":
            continue

        differences = new[col].ne(old[col])
        if not differences.any():
            columns_without_diff.append(col)
            continue

        print(f"Found {differences.sum()} differences in column '{col}':")
        if verbosity in ["all", "diff"] or col in VERBOSE_COLUMNS:
            diff_states = differences.loc[differences].index
            for state in diff_states:
                print(f"  State({state})")
                print(f"    New value: {new.loc[state, col]}")
                print(f"    Old value: {old.loc[state, col]}")
    print()

    # Compare numeric values within tolerance
    for col in COLUMNS_COMPARE_NUMERICAL:
        differences = (new[col] - old[col]).abs()
        tolerance = atol + rtol * old[col].abs()
        mask = differences.gt(tolerance)

        if mask.sum() == 0 and differences.max() < 1e-15:  # noqa: PLR2004
            columns_without_diff.append(col)
            continue

        print(f"Found {mask.sum()} {col} differences outside tolerance:")
        if (verbosity in ["all", "diff"] and mask.any()) or col in VERBOSE_COLUMNS:
            diff_states = mask.loc[mask].index
            for state in diff_states:
                new_val = new.loc[state, col]
                old_val = old.loc[state, col]
                diff = differences.loc[state]
                print(f"  State({state}):")
                print(f"    New {col}: {new_val:.12f}")
                print(f"    Old {col}: {old_val:.12f}")
                print(f"    Absolute difference: {diff:.2e}")
                print(f"    Relative difference: {diff / abs(old_val):.2e}")

        max_rdiff = (differences / old[col].abs()).max()
        print(f"  Maximum absolute {col} difference: {differences.max():.2e}")
        print(f"  Maximum relative {col} difference: {max_rdiff:.2e}")
    print()

    print(f"Columns without differences: {', '.join(columns_without_diff)}")


if __name__ == "__main__":
    main()
