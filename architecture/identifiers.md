# Identifiers in gDRutils

Identifiers map gDR's logical fields (`drug_name`, `duration`, …) to the column names a dataset
actually uses (`DrugName`, `Duration`, …). Always read them through these functions — never
hardcode a column name.

Source of truth: `R/identifiers_list.R`, `R/global_cache.R`, `R/identifiers.R`, `R/se_metadata.R`,
`R/merge_SE.R`, `R/standardize_MAE.R` and `R/validate_identifiers.R` here, plus the callers named
below in gDRimport and gDRcore. Checked against `main` of gDRutils `d48580a0`, gDRimport `52e7a93`
and gDRcore `bce0cc58`.

## Where the values live

```mermaid
flowchart LR
  subgraph pkg["gDRutils (package constants)"]
    IL["IDENTIFIERS_LIST<br/>defaults, e.g. drug_name = DrugName"]
    SYN["SYNONYMS_LIST<br/>alternative column names"]
    EO["EXPECT_ONE_IDENTIFIERS"]
  end

  subgraph session["R session"]
    GC[("global_cache$identifiers_list<br/>empty until first read")]
  end

  subgraph obj["Data objects"]
    SE[("SE metadata<br/>'identifiers'")]
    MAE[("MAE = one SE per experiment")]
  end

  IL -- "first get_env_identifiers()" --> GC
  IL -- "reset_env_identifiers()" --> GC
  U(["caller"]) -- "set_env_identifier(k, v)" --> GC
  GC -- "get_env_identifiers(k)" --> U
  GC -- "set_SE_identifiers() at create_SE" --> SE
  SE --> MAE
  MAE -- "update_env_idfs_from_mae()" --> GC
  SE -- "merge_SE(): identifiers of the first SE" --> GC
  SE -- "get_SE_identifiers(se, k)" --> U
  GC -. "fallback with a warning when<br/>k is missing on the SE" .-> SE
```

## How the suite uses them

```mermaid
flowchart LR
  subgraph utils["gDRutils"]
    SYN["SYNONYMS_LIST"]
    GC[("session identifiers")]
    STD["standardize_se()"]
  end

  subgraph imp["gDRimport"]
    LF["load_files<br/>update_idfs_synonyms(colnames)<br/>validate_identifiers(manifest)"]
    PR["PRISM converters<br/>convert_LEVEL5/6_prism_to_gDR_input()"]
  end

  subgraph core["gDRcore"]
    MD["merge_data<br/>validate_identifiers()"]
    CS["create_SE<br/>set_SE_identifiers()"]
    NS["normalize_SE, average_SE, ...<br/>mostly get_env_identifiers(k)"]
    PL["runDrugResponseProcessingPipeline()"]
  end

  SE[("SE metadata")]

  SYN --> LF
  GC --> LF
  GC --> MD --> CS --> SE
  GC --> NS
  SE -. "a few get_SE_identifiers(se, k)" .-> NS
  PR == "reset_env_identifiers()<br/>on every call" ==> GC
  PL == "reset at the end of every run" ==> GC
  PL -- "MAE input" --> STD
  STD == "reset_env_identifiers()" ==> GC
```

## Things the diagrams make visible

- **Processing reads the session, not the SE.** gDRcore calls `get_env_identifiers()` about 130
  times and `get_SE_identifiers()` 8 times. `average_SE` reads no SE identifiers, and
  `normalize_SE` takes the duration from the SE at `gDRcore/R/normalize_SE.R:105` and from the
  session at `:355`. An SE with custom identifiers keeps them in its metadata, but most of the
  processing still uses the session values.
- **Four places reset the session and do not restore it:** the PRISM converters
  (`gDRimport/R/prism_to_gdrDF.R:37`, `:275`), the end of every pipeline run
  (`gDRcore/R/runDrugResponseProcessingPipeline.R:443`) and `standardize_se()`
  (`R/standardize_MAE.R:44`), which the pipeline calls for MAE input
  (`runDrugResponseProcessingPipeline.R:291`). Identifiers a caller set with
  `set_env_identifier()` are gone after any of them. (`setEnvForPSet()` resets them too, but
  setting them is its documented purpose.) GDR-3862 replaces these resets with
  `local_identifiers()`, one package at a time; the gDRcore and gDRimport ones are next.
- **`merge_SE()` writes to the session:** it copies the identifiers of the first SE into it
  (`R/merge_SE.R:173-176`). GDR-3864 replaces this with a scoped, harmonized merge.
- **Fallback, not failure.** `get_SE_identifiers()` falls back to the session value with a warning
  when an identifier is missing on the SE (`R/se_metadata.R:200-203`), so an old SE can be read
  with today's defaults without anyone noticing.

## Scoped identifiers (gDRutils 1.11.32)

`local_identifiers(ids)` makes `ids` the active identifiers until the calling function returns,
on success or error, then restores the previous ones. `with_identifiers()` does the same around
one expression, and `with_SE_identifiers(se, code)` uses the identifiers stored on `se`.

```mermaid
flowchart LR
  SE[("SE metadata<br/>'identifiers'")]
  GC[("session identifiers")]
  EP["SE entry point<br/>convert_se_assay_to_dt(), validate_SE(), apply_fit(), ..."]
  R["readers beneath it<br/>get_env_identifiers(k)"]
  SE -- ".local_SE_identifiers(se)<br/>on entry" --> GC
  GC --> R
  EP --> R
  EP -. "on exit: previous identifiers restored" .-> GC
```

- **The SE's identifiers take precedence** in the entry points listed below. A key the SE does not
  store keeps the session's value, as `get_SE_identifiers()` always did, and `same_identifiers()`
  compares a missing key as the session's value too. `set_SE_identifiers()` corrects an object;
  `set_env_identifier()` no longer overrides what an SE stores in these functions.
- **Scoped entry points** (gDRutils): `convert_se_assay_to_dt()`, `convert_se_assay_to_custom_dt()`,
  `convert_mae_assay_to_dt()` (through the former), `validate_SE()`, `apply_fit()`, `apply_fits()`.
  Everything else still reads the session.
- **Mixed reads remain** where a function calls a scoped entry point and then reads the session
  itself, e.g. `merge_assay()` (gDRutils, GDR-3864) and gDRcore's `.compute_growth_rates()`
  (`fit_SE.timecourse.R:578-584`). They agree for objects with default identifiers; the gDRcore
  follow-up of GDR-3862 scopes them.
- **`validate_SE()` checks** that the stored `cellline`, `cellline_name`, `drug`, `drug_name` and
  `concentration` name columns the SE has, so a stale identifier fails loudly instead of being
  read silently with the wrong column.

## Keeping it current

A change to how identifiers are stored, set or read updates these diagrams in the same pull
request, together with the commits above. Diagrams are [Mermaid](https://mermaid.js.org/); GitHub
renders them, and `npx @mermaid-js/mermaid-cli -i identifiers.md -o out.md` checks the syntax
locally.
