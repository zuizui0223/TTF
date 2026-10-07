args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("usage: export_historical_host_wcvp_names.R <panel_csv> <wcvp_repo> <output_csv>")
}
panel_path <- args[[1]]
wcvp_repo <- args[[2]]
out_path <- args[[3]]

panel <- read.csv(panel_path, stringsAsFactors=FALSE, check.names=FALSE)
ids <- unique(unlist(strsplit(as.character(panel$accepted_host_ids), ";", fixed=TRUE)))
ids <- sort(ids[nzchar(ids)])

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_names.rda"), envir=e)
stopifnot(exists("wcvp_names", envir=e))
n <- get("wcvp_names", envir=e)

n$plant_name_id <- as.character(n$plant_name_id)
n$accepted_plant_name_id <- as.character(n$accepted_plant_name_id)
status <- tolower(trimws(as.character(n$taxon_status)))
missing_acc <- is.na(n$accepted_plant_name_id) | !nzchar(n$accepted_plant_name_id)
n$accepted_plant_name_id[missing_acc & status=="accepted"] <- n$plant_name_id[missing_acc & status=="accepted"]

accepted <- n[
  n$plant_name_id %in% ids &
  tolower(trimws(as.character(n$taxon_rank)))=="species",
  c("plant_name_id","taxon_name","taxon_status","accepted_plant_name_id"),
  drop=FALSE
]
if (nrow(accepted) != length(ids)) {
  stop(sprintf("WCVP host ID coverage mismatch: expected %d got %d", length(ids), nrow(accepted)))
}
if (any(duplicated(accepted$plant_name_id))) stop("duplicate WCVP plant_name_id")
if (any(accepted$accepted_plant_name_id != accepted$plant_name_id)) {
  stop("frozen accepted host ID is no longer self-accepted in exact WCVP snapshot")
}
out <- data.frame(
  accepted_host_id=as.character(accepted$plant_name_id),
  accepted_host_name=as.character(accepted$taxon_name),
  stringsAsFactors=FALSE
)
out <- out[order(out$accepted_host_id),]
dir.create(dirname(out_path), recursive=TRUE, showWarnings=FALSE)
write.csv(out, out_path, row.names=FALSE, na="")
cat(sprintf("hosts=%d\n", nrow(out)))
