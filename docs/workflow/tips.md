# tips/tips

This module implements the `tipsFilter` process, which uses the [TIPS
CLI](../cli/filter.md) to perform different filtering of the dataset. The `inp`
key from the inputs is parsed to the CLI to specify the algorithm for filtering.

## tipsFilter

| Inputs | Default | Description      |
|--------|---------|------------------|
| inp    | `''`    | TIPS CLI options |
| ds     | `null`  | dataset          |

Returns the filtered dataset.
