args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("usage: export_historical_host_native_ranges_from_sidecar.R <panel_csv> <frozen_native_csv> <wcvp_repo> <output_geojson>")
}
panel_path <- args[[1]]
native_path <- args[[2]]
wcvp_repo <- args[[3]]
out_path <- args[[4]]

suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

panel <- read.csv(panel_path, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot("accepted_host_ids" %in% names(panel))
ids <- unique(unlist(strsplit(as.character(panel$accepted_host_ids), ";", fixed=TRUE)))
ids <- sort(trimws(ids[nzchar(trimws(ids))]))

native <- read.csv(native_path, stringsAsFactors=FALSE, check.names=FALSE)
required_native <- c("accepted_plant_name_id","area_code_l3","accepted_name")
stopifnot(all(required_native %in% names(native)))
native$accepted_plant_name_id <- trimws(as.character(native$accepted_plant_name_id))
native$area_code_l3 <- trimws(as.character(native$area_code_l3))
native$accepted_name <- trimws(as.character(native$accepted_name))
native <- native[
  native$accepted_plant_name_id %in% ids &
  nzchar(native$area_code_l3),
  required_native,
  drop=FALSE
]
seen <- sort(unique(native$accepted_plant_name_id))
if (!identical(seen, ids)) {
  missing <- setdiff(ids, seen)
  extra <- setdiff(seen, ids)
  stop(sprintf(
    "exact sidecar host-ID coverage mismatch: missing=%s extra=%s",
    paste(missing,collapse=";"), paste(extra,collapse=";")
  ))
}

e <- new.env(parent=emptyenv())
load(file.path(wcvp_repo, "data", "wgsrpd3.rda"), envir=e)
stopifnot(exists("wgsrpd3", envir=e))
areas <- get("wgsrpd3", envir=e)
stopifnot(all(c("LEVEL3_COD","geometry") %in% names(areas)))
areas$LEVEL3_COD <- trimws(as.character(areas$LEVEL3_COD))

rows <- list()
for (id in ids) {
  d <- native[native$accepted_plant_name_id == id,,drop=FALSE]
  codes <- sort(unique(d$area_code_l3[nzchar(d$area_code_l3)]))
  names_i <- sort(unique(d$accepted_name[nzchar(d$accepted_name)]))
  if (length(codes)==0) stop(paste("sidecar host lacks primary-native range:",id))
  if (length(names_i)!=1) stop(paste("sidecar host canonical-name ambiguity:",id))
  p <- areas[areas$LEVEL3_COD %in% codes,,drop=FALSE]
  if (nrow(p)==0) stop(paste("WGSRPD3 polygons missing for sidecar host:",id))
  geom <- st_union(st_make_valid(st_geometry(p)))
  rows[[length(rows)+1]] <- st_sf(
    accepted_host_id=id,
    accepted_host_name=names_i[[1]],
    native_area_count=length(codes),
    geometry=st_sfc(geom,crs=st_crs(areas))
  )
}
out <- do.call(rbind,rows)
dir.create(dirname(out_path),recursive=TRUE,showWarnings=FALSE)
st_write(out,out_path,driver="GeoJSON",delete_dsn=TRUE,quiet=TRUE)
cat(sprintf("hosts=%d\n",nrow(out)))
