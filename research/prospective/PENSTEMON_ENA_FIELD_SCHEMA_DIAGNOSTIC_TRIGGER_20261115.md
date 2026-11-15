# ENA field schema only, no biological data

The public ENA 'returnFields?result=sample' endpoint was requested but
could not be parsed as JSON or ordinary TSV. The next action is restricted
to a one-shot diagnosis of the schema response envelope (content-type,
first two field-description lines) and does not open sample or phenotype
records. No candidate model is modified.
