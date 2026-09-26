args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("usage: build_wcvp_hosts_contemporary_sidecar.R <wcvp_repo> <hosts_repo> <output_dir>")
}
wcvp_repo <- args[[1]]
hosts_repo <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_names.rda"), envir=e)
load(file.path(wcvp_repo, "data", "wcvp_distributions.rda"), envir=e)
stopifnot(exists("wcvp_names", envir=e), exists("wcvp_distributions", envir=e))
names_df <- get("wcvp_names", envir=e)
dist_df <- get("wcvp_distributions", envir=e)

hosts <- read.csv(
  file.path(hosts_repo, "resource.csv"),
  stringsAsFactors=FALSE,
  check.names=FALSE,
  na.strings=c("")
)
required_hosts <- c("Hostplant Genus","Hostplant Species")
stopifnot(all(required_hosts %in% names(hosts)))

genus <- trimws(ifelse(is.na(hosts[["Hostplant Genus"]]), "", hosts[["Hostplant Genus"]]))
species <- trimws(ifelse(is.na(hosts[["Hostplant Species"]]), "", hosts[["Hostplant Species"]]))
keep <- nzchar(genus) & nzchar(species)
host_names <- sort(unique(paste(genus[keep], species[keep])))

required_names <- c(
  "plant_name_id","taxon_rank","taxon_status","family","taxon_name",
  "accepted_plant_name_id"
)
stopifnot(all(required_names %in% names(names_df)))

rank_species <- tolower(trimws(as.character(names_df$taxon_rank))) == "species"
exact <- names_df[
  rank_species & as.character(names_df$taxon_name) %in% host_names,
  required_names,
  drop=FALSE
]
exact$plant_name_id <- as.character(exact$plant_name_id)
exact$accepted_plant_name_id <- as.character(exact$accepted_plant_name_id)

status_lower <- tolower(trimws(as.character(exact$taxon_status)))
missing_acc <- is.na(exact$accepted_plant_name_id) | !nzchar(exact$accepted_plant_name_id)
exact$accepted_plant_name_id[missing_acc & status_lower=="accepted"] <-
  exact$plant_name_id[missing_acc & status_lower=="accepted"]

id_to_name <- setNames(
  as.character(names_df$taxon_name),
  as.character(names_df$plant_name_id)
)
id_to_family <- setNames(
  as.character(names_df$family),
  as.character(names_df$plant_name_id)
)

split_ids <- split(exact$accepted_plant_name_id, as.character(exact$taxon_name))
resolved <- lapply(host_names, function(nm) {
  ids <- unique(split_ids[[nm]])
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids)==0) {
    return(data.frame(
      input_host_name=nm, match_status="unmatched",
      accepted_plant_name_id="", accepted_name="", family="",
      stringsAsFactors=FALSE
    ))
  }
  if (length(ids)>1) {
    return(data.frame(
      input_host_name=nm, match_status="ambiguous",
      accepted_plant_name_id=paste(sort(ids), collapse=";"),
      accepted_name="", family="", stringsAsFactors=FALSE
    ))
  }
  id <- ids[[1]]
  data.frame(
    input_host_name=nm, match_status="matched",
    accepted_plant_name_id=id,
    accepted_name=ifelse(is.na(id_to_name[[id]]),"",id_to_name[[id]]),
    family=ifelse(is.na(id_to_family[[id]]),"",id_to_family[[id]]),
    stringsAsFactors=FALSE
  )
})
resolved <- do.call(rbind, resolved)

insect_genus <- trimws(ifelse(is.na(hosts[["Insect Genus"]]), "", hosts[["Insect Genus"]]))
insect_species_epithet <- trimws(ifelse(is.na(hosts[["Insect Species"]]), "", hosts[["Insect Species"]]))
insect_keep <- keep & nzchar(insect_genus) & nzchar(insect_species_epithet)
interaction <- data.frame(
  insect_species=paste(insect_genus[insect_keep], insect_species_epithet[insect_keep]),
  input_host_name=paste(genus[insect_keep], species[insect_keep]),
  stringsAsFactors=FALSE
)
interaction <- unique(interaction)
interaction <- merge(
  interaction,
  resolved[resolved$match_status=="matched",
           c("input_host_name","accepted_plant_name_id","accepted_name","family"),
           drop=FALSE],
  by="input_host_name",
  all=FALSE,
  sort=FALSE
)
interaction <- unique(interaction)
interaction <- interaction[
  order(interaction$insect_species, interaction$accepted_plant_name_id),
  ,
  drop=FALSE
]
write.csv(
  interaction,
  file.path(output_dir, "insect_host_accepted.csv"),
  row.names=FALSE, na=""
)

matched_ids <- unique(resolved$accepted_plant_name_id[resolved$match_status=="matched"])
required_dist <- c(
  "plant_name_id","area_code_l3","introduced","extinct","location_doubtful"
)
stopifnot(all(required_dist %in% names(dist_df)))

# Contemporary resource opportunity deliberately differs from the primary S2
# native-only layer in exactly one dimension: introduced ranges are retained.
# Extinct and location-doubtful records remain excluded.
d <- dist_df[
  as.character(dist_df$plant_name_id) %in% matched_ids &
  as.integer(dist_df$extinct)==0L &
  as.integer(dist_df$location_doubtful)==0L,
  required_dist,
  drop=FALSE
]
d$plant_name_id <- as.character(d$plant_name_id)
d$area_code_l3 <- as.character(d$area_code_l3)
d <- unique(d[c("plant_name_id","area_code_l3")])
d <- d[nzchar(d$area_code_l3),,drop=FALSE]
names(d)[1] <- "accepted_plant_name_id"
d$accepted_name <- unname(id_to_name[d$accepted_plant_name_id])
d <- d[order(d$accepted_plant_name_id,d$area_code_l3),]
write.csv(
  d,
  file.path(output_dir, "extant_nondoubtful_wgsrpd3.csv"),
  row.names=FALSE, na=""
)

native <- dist_df[
  as.character(dist_df$plant_name_id) %in% matched_ids &
  as.integer(dist_df$introduced)==0L &
  as.integer(dist_df$extinct)==0L &
  as.integer(dist_df$location_doubtful)==0L,
  required_dist,
  drop=FALSE
]
native <- unique(native[c("plant_name_id","area_code_l3")])

counts <- data.frame(
  key=c(
    "hosts_resource_rows",
    "distinct_species_level_host_binomials",
    "resolved_insect_host_pairs",
    "contemporary_host_x_wgsrpd3_rows",
    "native_host_x_wgsrpd3_rows",
    "contemporary_minus_native_rows"
  ),
  value=c(
    nrow(hosts),
    length(host_names),
    nrow(interaction),
    nrow(d),
    nrow(native),
    nrow(d)-nrow(native)
  )
)
write.table(
  counts,
  file.path(output_dir, "counts.tsv"),
  row.names=FALSE, quote=FALSE, sep="\t"
)
