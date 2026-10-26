# popmix populations

List every population of a dataset with its number of samples, optionally
filtered by a wildcard pattern. Use it to find the exact population names
other commands need.

## Synopsis

```text
popmix populations PREFIX [--filter PATTERN] [--format text|json|tsv]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `PREFIX` | text (path) | required | Genotype prefix, or the full path of a VCF. |
| `--filter`, `-f` | text | none (all) | Shell-style wildcard pattern (`fnmatch`: `*`, `?`, `[...]`), case-sensitive, matched against the whole name: `"French*"`, `"*_N"`, `"*.DG"`. Quote it so the shell does not expand it. |
| `--format` | `text`, `json`, `tsv` | `text` | `text`: an aligned table and tips; `json`: a list of `{"population", "samples"}`; `tsv`: `population` and `samples` columns. |

## Inputs and outputs

Reads the `.ind` / `.fam` / VCF header. Writes to stdout. The populations
are sorted by sample count, largest first.

## Examples

```bash
# runs
popmix populations family
popmix populations family --format tsv
popmix populations family --filter "Sim*" --format json
```

```bash
# not run: needs the AADR
popmix populations /references/aadr/v62.0_1240k_public --filter "*Yamnaya*"
```

## Reading the result

The text format ends with tips that show how to pass names to other
commands. For scripts use `--format tsv` or `json`.

## Caveats

- A misspelled population in any other command already prints the closest
  names; `populations --filter` is for browsing.
- The pattern is case-sensitive (`french*` does not match `French`).

## Reference

None; the equivalent of counting the third column of a `.ind` file.
