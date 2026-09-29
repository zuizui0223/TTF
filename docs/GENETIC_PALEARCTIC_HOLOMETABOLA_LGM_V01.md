# Palearctic terrestrial insects × LGM refugia — response-blind program v0.1

## Question

> Do terrestrial Palearctic insect species that shared Last Glacial Maximum refugial geography tend to show genetic differentiation in the same places today?

This is a new ecological question derived from the failure mode exposed by the broad 211-species transfer analysis. It is **not** a rescue or subset retest of the original confirmatory claim.

## Why this is biologically narrower

The parent 211-species panel mixes many continents and life histories. The new program restricts the design before constructing any new pairwise genetic response:

- class: Insecta;
- orders: Coleoptera, Hymenoptera, Lepidoptera;
- at least 80% of the frozen Phase-1 COI localities inside the Ecoregions2017 Palearctic realm;
- no backfilling after the candidate count is known.

A response-blind conservative pre-screen found 48 plausible species (29 Lepidoptera, 11 Coleoptera, 8 Hymenoptera). This is feasibility evidence only; the official realm polygon decides the final panel.

## Historical predictor

For every eligible species, an external GBIF occurrence set is thinned at 10 km and used to fit a current-climate SDM from four CHELSA climate variables. The fitted current niche is projected to CHELSA-TraCE21k at 21 ka BP.

The primary historical relation between two species is continuous Schoener overlap of their normalized LGM suitability surfaces. A thresholded refugial Jaccard index is sensitivity analysis only.

## Why TTF-Q comes first

Before constructing the pairwise genetic response, TTF-Q asks whether LGM refugial similarity actually contains usable pairwise information:

1. how much survives exact species endpoint identity;
2. how much survives current niche similarity, current geography, taxonomy and sampling imbalance;
3. whether residual information is concentrated in only a few species;
4. whether the intended two-endpoint inference is null-qualified on the actual subset geometry;
5. which effect sizes are evaluable conditional on null qualification.

Only a TTF-Q GO result licenses construction of the new genetic-concordance response.

## Response if qualified

For an unordered pair of species, genetic-geography concordance is the average of two directional single-source TTF scores:

- train the frozen post-IBD spatial kernel on species i and score species j;
- train on j and score i;
- average the two directional scores.

This asks directly whether one species' post-IBD differentiation geography predicts another species' geography.

The hypothesis test is then:

> LGM refugial similarity positively predicts symmetric pairwise genetic-geography concordance, conditional on current niche similarity, present geographic opportunity, order and sampling geometry.

## Inference boundary

The aggregate Phase-4 result for the 211-species parent dataset is already known. Therefore this cannot be described as a never-opened confirmatory dataset.

What remains shielded is the **new species-pair response matrix** for this historical hypothesis. Panel selection, SDM construction and TTF-Q qualification are frozen before that matrix is constructed or inspected.

A failed qualification closes the program as NOT_EVALUABLE. It is not a biological null.
