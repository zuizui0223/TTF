#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sf)
  library(digest)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: export_wgsrpd3_support.R <wgsrpd3.rda> <output_dir>")
rda <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

env <- new.env(parent=emptyenv())
load(rda, envir=env)
if (!exists("wgsrpd3", envir=env, inherits=FALSE)) stop("wgsrpd3 object missing")
x <- get("wgsrpd3", envir=env)
if (!inherits(x, "sf")) stop("wgsrpd3 is not sf")
if (is.na(st_crs(x))) st_crs(x) <- 4326
x <- st_transform(x, 4326)

candidates <- c("LEVEL3_COD","LEVEL3_CODE","level3_cod","level3_code","area_code_l3")
code_col <- candidates[candidates %in% names(x)][1]
if (is.na(code_col)) stop(paste("cannot find WGSRPD3 code column; fields:", paste(names(x), collapse=",")))
codes <- trimws(as.character(x[[code_col]]))
if (any(!nzchar(codes)) || anyDuplicated(codes)) stop("WGSRPD3 codes missing or duplicated")

lon <- seq(-179.75, 179.75, by=0.5)
lat <- seq(-89.75, 89.75, by=0.5)
grid_df <- expand.grid(longitude=lon, latitude=lat, KEEP.OUT.ATTRS=FALSE)
grid <- st_as_sf(grid_df, coords=c("longitude","latitude"), crs=4326, remove=FALSE)
inside <- st_intersects(x, grid)

fmt <- function(v) sprintf("%.2f", as.numeric(v))
rows <- vector("list", nrow(x))
for (i in seq_len(nrow(x))) {
  code <- codes[[i]]
  idx <- inside[[i]]
  if (length(idx) == 0L) {
    p <- suppressWarnings(st_point_on_surface(st_geometry(x[i,])))
    xy <- st_coordinates(p)[1,c("X","Y")]
    rows[[i]] <- data.frame(
      area_code_l3=code,
      point_index=1L,
      latitude=as.numeric(xy[["Y"]]),
      longitude=as.numeric(xy[["X"]]),
      selection_mode="point_on_surface_fallback",
      selection_sha256="",
      stringsAsFactors=FALSE
    )
    next
  }
  pts <- grid_df[idx,,drop=FALSE]
  keys <- paste(
    "hhm-wgsrpd-sample-v0.1", code, fmt(pts$latitude), fmt(pts$longitude),
    sep="|"
  )
  hashes <- vapply(keys, digest, "", algo="sha256", serialize=FALSE)
  ord <- order(hashes, pts$latitude, pts$longitude)
  if (length(ord) > 25L) ord <- ord[seq_len(25L)]
  pts <- pts[ord,,drop=FALSE]
  hashes <- hashes[ord]
  rows[[i]] <- data.frame(
    area_code_l3=code,
    point_index=seq_len(nrow(pts)),
    latitude=pts$latitude,
    longitude=pts$longitude,
    selection_mode="fixed_0p5deg_hash_cap",
    selection_sha256=hashes,
    stringsAsFactors=FALSE
  )
}
out <- do.call(rbind, rows)
out <- out[order(out$area_code_l3, out$point_index),]
write.csv(out, file.path(outdir,"wgsrpd3_support_points_v0.1.csv"), row.names=FALSE, na="")

summary <- data.frame(
  key=c("wgsrpd3_units","support_points","fallback_units","min_points_per_unit","median_points_per_unit","max_points_per_unit"),
  value=c(
    nrow(x),
    nrow(out),
    sum(out$selection_mode=="point_on_surface_fallback"),
    min(table(out$area_code_l3)),
    median(as.numeric(table(out$area_code_l3))),
    max(table(out$area_code_l3))
  )
)
write.table(summary, file.path(outdir,"counts.tsv"), row.names=FALSE, quote=FALSE, sep="\t")
cat(paste(code_col, nrow(x), nrow(out), sep="\t"),"\n")
