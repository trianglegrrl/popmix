# popmix install-skill

Install this agent skill (SKILL.md and one file per command, as shipped
inside the installed popmix package) where Claude Code finds personal
skills, `~/.claude/skills/popmix`, or into any directory.

## Synopsis

```text
popmix install-skill [--dest DIR] [--force]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--dest` | text (path) | `~/.claude/skills/popmix` | Directory to install into (`~` is expanded). For a project skill, use `.claude/skills/popmix` in the project. |
| `--force` | flag | off | Replace an earlier copy of the popmix skill at `--dest` (for example after upgrading popmix). Without it an existing copy is left alone and the command exits 1. |

## Inputs and outputs

Copies the skill files from the package (`popmix/skill/`) into `--dest`
and prints `Installed the popmix <version> skill (<n> files) to DEST`. An empty
directory at `--dest` is filled; a directory holding something other than a
popmix skill is never replaced, even with `--force`.

## Examples

```bash
# runs
popmix install-skill --dest skills/popmix
popmix install-skill --dest skills/popmix --force
ls skills/popmix/commands
```

```bash
# not run: writes to your home directory
popmix install-skill
popmix install-skill --force
popmix install-skill --dest .claude/skills/popmix
```

## Reading the result

Claude Code loads the skill in new sessions; it states the popmix version
it documents (`This skill documents popmix X.` in SKILL.md), which should
match `popmix --version`.

## Caveats

- The skill is versioned with popmix: after `pip install -U popmix`, run
  `popmix install-skill --force` so the documentation matches the CLI.
- `--force` replaces the whole directory, including files added to it.

## Reference

None; Claude Code's skills directory layout (`~/.claude/skills/<name>/SKILL.md`).
