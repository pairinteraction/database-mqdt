from pathlib import Path

import pandas as pd

from utils import add_global_uid_column

TABLE_NAMES: list[str] = [
    "matrix_elements_d",
    "matrix_elements_q",
    "matrix_elements_o",
    "matrix_elements_q0",
    "matrix_elements_mu",
]


def main() -> None:
    # CHANGE THESE PATHS, TO THE FOLDERS YOU WANT TO COMPARE
    db_dir = Path(__file__).parent.parent.parent / "database"
    species = "Yb174_mqdt"
    old_path = db_dir / f"{species}_v1.2"
    new_path = db_dir / f"{species}_v1.3"

    print(f"Comparing matrix elements tables:\n  New: {new_path}\n  Old: {old_path}")
    for table_name in TABLE_NAMES:
        if not (new_path / f"{table_name}.parquet").exists() or not (old_path / f"{table_name}.parquet").exists():
            print(f"\nSkipping {table_name} as it does not exist in either the new or the old path.\n")
            continue
        compare_matrix_elements_table(table_name, new_path, old_path)


def compare_matrix_elements_table(  # noqa: C901, PLR0915
    table_name: str,
    new_path: Path,
    old_path: Path,
    rtol: float = 1e-2,
    atol: float = 1e-5,
    *,
    min_nu: float = 0,
    max_nu: float = float("inf"),
    max_delta_nu: int = 3,
    only_compare_absolute_values: bool = False,
    verbose: bool = False,
) -> None:
    """Compare the matrix elements table of two versions of the database.

    Given the path to states and matrix elements parquet tables,
    this function will
    1) create a new index mapping from a state (defined by n, exp_l, and exp_j) to a unique identifier
    2) replace the old id_initial and id_final column in the matrix elements table with this new index
    3) compare the column "val" of the two matrix elements tables.

    """
    print(f"Comparing matrix elements table:     {table_name}")

    species = new_path.name.split("_v")[0]
    if species != old_path.name.split("_v")[0]:
        raise ValueError(f"Cannot compare different species: {species} vs {old_path.name.split('_v')[0]}")

    paths = {"new": new_path, "old": old_path}
    states_dict = {key: pd.read_parquet(path / "states.parquet") for key, path in paths.items()}
    table_dict = {key: pd.read_parquet(path / f"{table_name}.parquet") for key, path in paths.items()}
    print(f"  Table shape pre-filtering:    New: {table_dict['new'].shape}; Old: {table_dict['old'].shape}")

    all_columns = ["id_initial", "id_final", "val"]
    for key, table in table_dict.items():
        missing_cols = [col for col in all_columns if col not in table.columns]
        if len(missing_cols) > 0:
            print(f"WARNING: {key.capitalize()} table is missing columns: {missing_cols}")
        extra_cols = [col for col in table.columns if col not in all_columns and col != "id"]
        if len(extra_cols) > 0:
            print(f"WARNING: {key.capitalize()} table has extra columns: {extra_cols}")

    # Filter states and matrix elements tables by nu and delta_nu
    for key, states in states_dict.items():
        add_global_uid_column(species, states)
        id_to_global_uid = dict(zip(states["id"], states["global_uid"], strict=True))

        # Add columns global_uid_initial and global_uid_final to the matrix elements table
        table = table_dict[key]
        for which in ["initial", "final"]:
            table[f"global_uid_{which}"] = table[f"id_{which}"].map(id_to_global_uid)
            if table[f"global_uid_{which}"].isna().any():  # this should not happen
                raise ValueError(f"Some {which} ids in {key} table dont have a entry in the states table.")

        # Set new column nu_... and filter by it
        id_to_nu = dict(zip(states["id"], states["nu"], strict=True))
        for which in ["initial", "final"]:
            table[f"nu_{which}"] = table[f"id_{which}"].map(id_to_nu)
            table_dict[key] = table = table[table[f"nu_{which}"] >= min_nu].copy()
            table_dict[key] = table = table[table[f"nu_{which}"] <= max_nu].copy()

        # Set new column with delta_nu = abs(nu_final - nu_initial) and filter by it
        table["delta_nu"] = (table["nu_final"] - table["nu_initial"]).abs()
        table_dict[key] = table = table[table["delta_nu"] <= max_delta_nu].copy()

        # Index the matrix elements by the global state uid pair
        table = table.set_index(["global_uid_initial", "global_uid_final"], drop=False)
        table_dict[key] = table.sort_index()
    print(f"  Table shape after-filtering:  New: {table_dict['new'].shape}; Old: {table_dict['old'].shape}")

    # Only keep entries that are present in both tables
    new, old = table_dict["new"], table_dict["old"]
    common_index = new.index.intersection(old.index)
    for table in table_dict.values():
        table.drop(index=common_index.symmetric_difference(table.index), inplace=True, errors="ignore")  # noqa: PD002
        table.sort_index(inplace=True)  # noqa: PD002
    print(f"  Common entries: {new.shape[0]}")

    # Compare val values within tolerance
    val_str = "val"
    if only_compare_absolute_values:
        new["val"] = new["val"].abs()
        old["val"] = old["val"].abs()
        val_str = "|val|"

    val_diff = (new["val"] - old["val"]).abs()
    tolerance = atol + rtol * old["val"].abs()
    val_mask = val_diff.gt(tolerance)

    print(f"  Found {val_mask.sum()}/{val_mask.shape[0]} differences outside tolerance:")
    if verbose and val_mask.any():
        diff_uids = val_mask.loc[val_mask].index
        for uid in diff_uids:
            new_val = new.loc[uid, "val"]
            old_val = old.loc[uid, "val"]
            diff_val = val_diff.loc[uid]
            print(f"    State initial: {uid[0]}, State final: {uid[1]}")
            print(f"      New {val_str}: {new_val:.5f}, Old {val_str}: {old_val:.5f}")
            print(f"      Absolute difference: {diff_val:.2e}, Relative difference: {diff_val / abs(old_val):.2e}")

    rdiff = val_diff / old["val"].abs()
    print(f"  Maximum absolute {val_str} difference: {val_diff.max():.2e}")
    print(f"  Maximum relative {val_str} difference: {rdiff.max():.2e}")
    print()


if __name__ == "__main__":
    main()
