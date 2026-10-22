# Story: specvec

## Claim
Learns abundance-aware vectors for species and for vegetation plots from a long plot-species-cover table, placing both in one space, with a harness that scores embedding methods on neutral ecological tasks.

## Audience
Vegetation ecologists with plot data (EVA-style long tables, percent or Braun-Blanquet cover) who know ordination; R users.

## Relation
- to: none (correspondence analysis and GloVe are compared as methods, not as packages)

## Arcs

### A1 Species neighbours from co-occurrence
- question: which species does my data treat as ecological partners?
- answer: `species_embedding()` with the cover-weighted PMI default, then `nearest_species()`
- evidence: vignettes/specvec-quickstart.Rmd; dev_notes/demo_neophytes.R (EVA neighbour lists)
- home: vignettes/specvec-quickstart.Rmd (Getting started with specvec)
- also: README Usage, README neighbours section

### A2 Communities in the same space
- question: how similar or unusual are my plots?
- answer: `community_embedding()`, `community_similarity()`, `community_novelty()`
- home: vignettes/specvec-communities.Rmd (Community embeddings and novelty)

### A3 Change through time
- question: how do species and communities move, and is a neophyte integrating?
- answer: `species_trajectory()`, `community_trajectory()`, `integration_trajectory()`
- home: vignettes/specvec-temporal.Rmd, vignettes/specvec-integration.Rmd

### A4 Which method to trust
- question: is the default method the right one for my data?
- answer: `compare_embeddings()` under a fixed protocol and win criterion
- evidence: bakeoff/ results at 10k and 200k EVA plots; IMPLEMENTATION_PLAN.md section 13
- limits: GloVe needs text2vec; results shown in the README come from access-controlled EVA data
- home: vignettes/specvec-benchmark.Rmd (Compare embedding methods on your own data)

### A5 New methods by registration
- home: vignettes/specvec-extending.Rmd

## Documents
- README / landing: claim, toy example, one section per arc, routing
- articles: quickstart (tutorial), data, methods, communities, benchmark, temporal, integration, extending
