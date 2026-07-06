# Building an MDTF POD with Claude on NCAR HPC — a workflow guide

A companion to [`natl_ocean_notebook_user_guide.md`](./natl_ocean_notebook_user_guide.md). That
guide tells you how to *run* the `natl_ocean` POD. This one documents how the POD was *developed*
with **Claude Code** on NCAR Casper — the setup, the guardrails that keep the assistant honest, and
a library of reusable prompts you can adapt for your own POD.

> This reflects one researcher's setup as of mid-2026. Claude Code and its install steps change over
> time — verify commands against Anthropic's current Claude Code documentation rather than assuming
> the exact incantations below are still current.

---

## 1. What this is / getting started with Claude

**Claude Code** is Anthropic's coding agent. It reads and edits files in your repo, runs shell
commands, and works from a conversation. It runs in several places:

- a **terminal CLI**,
- an **IDE extension** (VS Code, JetBrains),
- the **desktop app** and **web** (claude.ai/code).

For POD development the two that matter are the **VS Code extension** and the **terminal CLI**, both
running against your checkout of `MDTF-diagnostics` on Casper.

Getting started, at a high level (check the official docs for exact steps):

1. Install Claude Code (CLI and/or the VS Code extension) and authenticate with your Anthropic
   account.
2. Launch it **from inside the project directory** you want to work in (this matters — see the
   per-project setup in §3). For this repo that is your `MDTF-diagnostics` checkout.
3. Describe the task. Claude explores the code, proposes a plan, and — once you approve — makes the
   edits and runs the checks.

The rest of this guide is about making that loop *reliable* on HPC.

---

## 2. Working on NCAR HPC (Casper)

Casper has a few properties that trip up a naive assistant. These are encoded once in a global
`CLAUDE.md` (§3) so Claude respects them every session:

- **Filesystem layout.** Home `/glade/u/home/<user>` (small, backed up), work
  `/glade/work/<user>` (code, envs), scratch `/glade/derecho/scratch/<user>` (large/temporary data).
  Keep large or throwaway outputs on scratch, not home.
- **Conda envs.** Prefer conda over global pip. Load with `module load conda`.
- **Casper conda gotcha (important).** `conda run -n <env> python|pip` is unreliable here — `PATH`
  can resolve `python`/`pip` from a different env that was prepended earlier, so a `-n` command
  silently runs in the wrong env. **Use absolute interpreter paths instead**, e.g.
  `/glade/u/home/<user>/miniconda3/envs/<env>/bin/python`. The `natl_ocean` work uses
  `/glade/u/home/taydral/miniconda3/envs/mdtf_nb/bin/python`.
- **Login node vs. compute node.** Login-node work is limited to editing, static checks
  (`nbformat` parse, `grep`, syntax), and doc generation. Anything that touches real data — running
  the POD, executing the notebooks, the CMORize step — must go through a **compute node**
  (PBS / `qsub`, or a `qvscode` interactive session). A good habit is to have Claude *flag* when a
  command needs a compute node instead of running it on the login node.

**Worked example — CMORize on a compute node.** The `cmorize_timeslice.py` step (native POP →
CMORized NetCDF) is real-data work, so it does not run on a login node. The pattern used here was:
Claude wrote/edited the script and validated it with `--dry-run` on the login node, then the actual
1-year build ran as a **PBS job** (`qsub`) on a compute node. Claude inspected the resulting files
and catalog afterward. See Step 3 (Option B) of the user guide for the science details.

### Using VS Code (and Claude inside it)

**VS Code** is Microsoft's free code editor. On HPC you typically run it on your laptop and connect
to Casper with the **Remote-SSH** extension (or an NCAR `qvscode` interactive session), so the editor
UI is local but the files, terminal, and Python all live on Casper.

Two ways to run Claude inside VS Code:

- the **Claude Code VS Code extension** — a side panel where you converse with Claude; it can open a
  **plan-review pane** (you approve/comment on a plan before any edits), and file references in its
  replies are **clickable** links straight to the line in your editor.
- the **integrated terminal** — run the Claude Code CLI in a VS Code terminal tab; same agent,
  terminal UI.

Why the IDE integration helps for this POD specifically: notebooks (`.ipynb`) are easier to review
cell-by-cell in the editor, the plan-review pane is useful for multi-file changes (e.g. stripping
both notebooks to a single data source), and clickable `file:line` references make it fast to jump to
exactly what Claude changed.

---

## 3. How Claude is configured: `CLAUDE.md`, notes, and memory

Claude Code reads a `CLAUDE.md` file at startup and treats it as standing instructions. This setup
uses **two tiers**:

- **Global** `~/CLAUDE.md` — user-wide conventions: the HPC layout and conda gotchas above, a
  response-style preference, a **projects index**, and the truth rules in §4. Applies everywhere.
- **Per-project** `CLAUDE.md` (e.g. in the repo root / project root) — project-specific layout,
  commands, and working rules that override or extend the global file where they conflict.

Two more mechanisms keep sessions oriented:

- **Per-project notes + a SessionStart hook.** Each project has a `notes/` directory of dated
  session summaries. A small **SessionStart hook** (in that project's
  `.claude/settings.local.json`) surfaces the newest notes automatically when a session starts in
  that project — so you resume with the prior decisions in view. Because it keys off the launch
  directory, launch Claude from the project you intend to work in.
- **Project-switch slash commands.** For moving between projects in one workspace, small global
  slash commands (`/mdtf`, `/osnap`, `/wmt`) load the right project's newest notes on demand — type
  the command and Claude reads that project's `notes/` before continuing.
- **Persistent memory.** Claude keeps a small file-based memory per project (under
  `~/.claude/projects/.../memory/`) for durable facts that aren't in the code or git history. It is
  updated deliberately, and recalled facts are treated as *possibly stale* — Claude re-verifies a
  named file/flag still exists before relying on it.

---

## 4. Keeping Claude honest (reduce hallucination)

The single most important part of this setup is a block of **truth rules** in the global
`~/CLAUDE.md`. They instruct Claude to prioritize accuracy over helpfulness and to flag uncertainty
instead of guessing. Reproduced here (adapt to your own file):

1. **Uncertainty** — if not fully certain, say so ("I am not certain, but…", "verify this…"). Never
   state guesses as facts.
2. **Sources** — do not invent paper titles, authors, URLs, or references. If no real, verifiable
   source can be named, say so.
3. **Statistics** — flag any number you are not fully confident in; say "approximately" and
   recommend verifying against a primary source.
4. **Recent events** — note when a topic may have changed since the knowledge cutoff; don't present
   outdated info as current.
5. **People & quotes** — never attribute a quote unless certain; otherwise say the quote can't be
   confirmed.
6. **Code & technical** — never invent function names, library methods, or API syntax; if unsure a
   function exists, say to verify it in current docs.
7. **Logic gaps** — don't fill missing context with assumptions; ask a clarifying question first.

In practice these show up as: Claude **verifying before asserting** (reading the file/running the
grep rather than recalling), **asking a clarifying question** when a request is ambiguous, and
**flagging** any command that should run on a compute node or any claim it can't confirm. Pairing the
rules with a habit of *running the check and showing the output* (rather than describing what the
code "should" do) is what keeps the POD work trustworthy.

A lightweight `/session-summary` command writes a detailed summary to `notes/` at the end of a
working session, which feeds the SessionStart hook next time — closing the loop.

---

## 5. Prompt library

Reusable prompts in a fixed shape:

```
Role:      who Claude should act as
Resources: the concrete files/docs to ground the work (this is what reduces hallucination)
Task:      the specific thing to do
Context:   constraints, environment, gotchas
Output:    what a good result looks like
```

Fill **Resources** with real paths — the more specific, the less Claude has to guess.

### 5.1 — CMORize the native timeslice on a compute node
```
Role:      An HPC-aware Python developer working on the MDTF natl_ocean POD.
Resources: diagnostics/natl_ocean/cmorize_timeslice.py; diagnostics/natl_ocean/POD_utils.py
           (load_native_timeslice); the native catalog CESM_timeslice_ocean_monthly_001.json;
           ~/CLAUDE.md (Casper conda + login/compute-node rules).
Task:      Build a 1-year CMORized timeslice bundle from the native POP output.
Context:   Casper. Real-data work → must run on a compute node via PBS/qsub, NOT the login node.
           Use the absolute interpreter /glade/u/home/<user>/miniconda3/envs/mdtf_nb/bin/python.
           Validate with --dry-run on the login node first.
Output:    A --dry-run that passes on the login node, then a qsub job script; after it runs,
           confirm the file count/size and that the intake-esm catalog opens the bundle by CMIP name.
```

### 5.2 — Review the POD for MDTF publish-readiness
```
Role:      An MDTF POD reviewer.
Resources: doc/sphinx/pod_requirements.rst, pod_settings.rst, dev_guidelines.rst,
           ref_catalogs.rst; the example_multicase and example_notebook reference PODs;
           the whole diagnostics/natl_ocean/ directory.
Task:      List what still blocks natl_ocean from being submittable, ranked by importance.
Context:   Target the NOAA-GFDL main branch conventions. Naming: driver/html/settings named after
           the POD short name. Flag unused/commented-out code and missing/incorrect docstrings.
Output:    A checklist of concrete fixes with file:line references; no edits yet — propose first.
```

### 5.3 — Add a data source to the notebook without hardcoding paths
```
Role:      A notebook developer for the natl_ocean POD.
Resources: diagnostics/natl_ocean/natl_ocean_esnb.ipynb (Section 1 case_info/DATA_CATALOG pattern);
           the intake-esm catalog JSON/CSV pair; ref_catalogs.rst.
Task:      Add a new model as a data_source option, resolved via the catalog (case_info['DATA_CATALOG']).
Context:   Do NOT hardcode file paths in cells — the catalog carries them. Keep the mode
           (interactive/prod) scaffolding intact. Match the existing per-variable dispatch.
Output:    A minimal cell diff that adds the branch; a note on which catalog/vars it needs.
```

### 5.4 — Set up a brand-new POD from scratch (no existing POD code)
```
Role:      An MDTF POD author starting from nothing.
Resources: doc/sphinx/pod_requirements.rst and pod_settings.rst; the example_multicase and
           example_notebook PODs as templates; ref_catalogs.rst.
Task:      Scaffold a new POD named <short_name> under diagnostics/<short_name>/.
Context:   Follow MDTF conventions: the driver, HTML template, and settings.jsonc are named after
           the POD short name; a doc .rst goes under the POD's doc/ dir. Nothing science-specific
           yet — just a correct, runnable skeleton.
Output:    diagnostics/<short_name>/ with <short_name>.py (or .ipynb), settings.jsonc, <short_name>.html,
           doc/<short_name>.rst — each matching the example POD's structure, with TODO markers.
```

### 5.5 — Convert an existing `.py` POD into notebook format
```
Role:      A notebook developer converting a driver-script POD to a demonstrable notebook.
Resources: the example_notebook POD (the canonical notebook template); this POD's
           natl_ocean_esnb.ipynb as a worked, section-structured example; the POD's existing
           driver .py; pod_requirements.rst (note: the framework runs notebooks via nbconvert,
           so notebook and .py implementations are near-identical).
Task:      Produce an .ipynb that mirrors the driver's logic, organized into the standard sections.
Context:   Keep the interactive/prod mode split so it runs both standalone and under the framework.
           Load data through a catalog, not hardcoded paths. Match the example notebook's section
           layout (Section 1 settings → Section 2 load → Section 3 compute → Section 4 plots).
Output:    A section-structured notebook that reproduces the driver's outputs; a short note on any
           logic that didn't translate cleanly.
```

### 5.6 — Finalize / clean before a PR
```
Role:      A careful reviewer preparing a branch to push.
Resources: git diff of the feature branch; ~/CLAUDE.md (git rules); the notebooks + guides in
           diagnostics/natl_ocean/.
Task:      Make the deliverables self-consistent (comments match code, guides match notebooks),
           then stage focused commits.
Context:   Don't push or force-push without explicit approval. Snapshot risky changes on an archive
           branch first. Verify notebooks parse (nbformat) and carry no dead references.
Output:    Clean, single-purpose commits with clear messages; a summary of what changed and any
           open items; nothing pushed.
```

---

## 6. See also

- [`natl_ocean_notebook_user_guide.md`](./natl_ocean_notebook_user_guide.md) — how to install,
  get the data, configure, and run the POD notebook.
- MDTF POD-development docs: `doc/sphinx/pod_requirements.rst`, `pod_settings.rst`,
  `dev_guidelines.rst`, `ref_catalogs.rst`.
- Reference PODs: `diagnostics/example_multicase/`, `diagnostics/example_notebook/`.
