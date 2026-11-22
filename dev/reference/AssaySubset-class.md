# An S4 class to manage subset representation.

Used inside the `Subset*` classes; not for direct use.

## Slots

- `subsetName`:

  `character(1)` Name of the subset.

- `rowIndices`:

  `vector("numeric")` Indices of the rows to include in the subset.

- `colIndices`:

  `vector("numeric")` Indices of the columns to include in the subset.

- `parentAssay`:

  `character(1)` Name of the parent of this subset.

- `internalAssay`:

  `SummarizedExperiment` An internal experiment object to store
  additional subset data.

- `altExpName`:

  `character(1)` or `NULL`. For a subset made with
  `createSubset(altExp =)`, the name of the altExp whose features its
  rows index; `NULL` for other subsets (ADR 0007).
