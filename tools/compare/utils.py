import pandas as pd

SQDT_INDEX_COLUMNS: dict[str, int | None] = {
    "n": None,
    "exp_l": None,
    "exp_j": None,
    "exp_s": None,
}
MQDT_INDEX_COLUMNS: dict[str, int | None] = {
    "nu": 2,
    "f": None,
    "exp_l": 1,
    "exp_j": 1,
    "exp_j_ryd": 1,
}


def add_global_uid_column(species: str, states: pd.DataFrame) -> None:
    """Return a copy of the states table indexed by a database-independent state key."""
    index_columns = get_index_columns(species)
    missing_columns = [col for col in index_columns if col not in states.columns]
    if missing_columns:
        raise ValueError(f"States table is missing index columns: {missing_columns}")

    states_copy = states.copy()
    for col, decimals in index_columns.items():
        if decimals is not None:
            states_copy[col] = states_copy[col].round(decimals)
    if "mqdt" in species:
        states_copy["exp_j_ryd"] = (2 * states_copy["exp_j_ryd"]).round(0) / 2

    global_uid = [
        ", ".join(f"{col}={state[i]}" for i, col in enumerate(index_columns))
        for state in states_copy[list(index_columns)].itertuples(index=False, name=None)
    ]
    global_index = pd.Index(global_uid, name="global_uid")
    if not global_index.is_unique:
        duplicate_descriptions = list(global_index[global_index.duplicated(keep=False)].drop_duplicates())[:5]
        details = "; ".join(duplicate_descriptions)
        if global_index.duplicated(keep=False).sum() > len(duplicate_descriptions):
            details += "; ..."
        raise ValueError(f"Duplicate global state index: {details}")

    states["global_uid"] = global_uid
    states.index = global_index


def get_index_columns(species: str) -> dict[str, int | None]:
    """Return the index columns for a given species."""
    if "mqdt" in species:
        return MQDT_INDEX_COLUMNS
    return SQDT_INDEX_COLUMNS
