args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("usage: export_historical_host_native_ranges.R <panel_csv> <wcvp_repo> <output_geojson>")
}
panel_path <- args[[1]]
wcvp_repo <- args[[2]]
out_path <- args[[3]]

suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

panel <- read.csv(panel_path, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("accepted_host_ids","accepted_host_names") %in% names(panel)))

ids <- unique(unlist(strsplit(as.character(panel$accepted_host_ids), ";", fixed=TRUE)))
ids <- ids[nzchar(ids)]

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wcvp_distributions.rda"), envir=e)
load(file.path(wcvp_repo, "data", "wcvp_names.rda"), envir=e)
load(file.path(wcvp_repo, "data", "wgsrpd3.rda"), envir=e)
stopifnot(exists("wcvp_distributions", envir=e), exists("wcvp_names", envir=e), exists("wgsrpd3", envir=e))
dist <- get("wcvp_distributions", envir=e)
names_df <- get("wcvp_names", envir=e)
areas <- get("wgsrpd3", envir=e)
names_df$plant_name_id <- as.character(names_df$plant_name_id)
id_to_name <- setNames(as.character(names_df$taxon_name), names_df$plant_name_id)

required_dist <- c("plant_name_id","area_code_l3","introduced","extinct","location_doubtful")
stopifnot(all(required_dist %in% names(dist)))
stopifnot(all(c("LEVEL3_COD","geometry") %in% names(areas)))

d <- dist[
  as.character(dist$plant_name_id) %in% ids &
  as.integer(dist$introduced)==0L &
  as.integer(dist$extinct)==0L &
  as.integer(dist$location_doubtful)==0L,
  required_dist,
  drop=FALSE
]
d$plant_name_id <- as.character(d$plant_name_id)
d$area_code_l3 <- as.character(d$area_code_l3)
d <- unique(d[c("plant_name_id","area_code_l3")])

rows <- list()
for (id in sort(ids)) {
  codes <- sort(unique(d$area_code_l3[d$plant_name_id == id]))
  codes <- codes[nzchar(codes)]
  if (length(codes)==0) stop(paste("host lacks frozen primary-native range:", id))
  p <- areas[as.character(areas$LEVEL3_COD) %in% codes, , drop=FALSE]
  if (nrow(p)==0) stop(paste("WGSRPD3 polygons missing for host:", id))
  geom <- st_union(st_make_valid(st_geometry(p)))
  rows[[length(rows)+1]] <- st_sf(
    accepted_host_id=id,
    native_area_count=length(codes),
    geometry=st_sfc(geom, crs=st_crs(areas))
  )
}
out <- do.call(rbind, rows)
dir.create(dirname(out_path), recursive=TRUE, showWarnings=FALSE)
st_write(out, out_path, driver="GeoJSON", delete_dsn=TRUE, quiet=TRUE)
cat(sprintf("hosts=%d\n", nrow(out)))
