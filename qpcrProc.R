library(data.table)

platePos <- fread("../randomised_plate_positions.csv")

qpcrQuants <- list.files(pattern = "*quant.csv")

geneList <- c("AOA_amoA", "AOB_amoA", "Bac_16S", "Anammox_hzo", "Anammox_hzsA", "nirK_C1", "nirK_C2", "nirK_C3", "nirS_C1", "nirS_C2", "nosZ_C1", "nosZ_C2", "Nitrospira_nxrB", "Thaum_ureC")

plate1 <- grep(pattern = "plate1", qpcrQuants, value = T)
plate2 <- grep(pattern = "plate2", qpcrQuants, value = T)

plate1Dat <- lapply(plate1, fread)
plate2Dat <- lapply(plate2, fread)

plate1Dat <- lapply(plate1Dat, function(x) {colnames(x)[13] <- "copies"; x})
plate2Dat <- lapply(plate2Dat, function(x) {colnames(x)[13] <- "copies"; x})

plate1Wells <- lapply(plate1Dat, function(x)
  x[grepl("Unkn", Content),
  .(well = Well[1], copies = copies[1]),
  by = Content])

plate2Wells <- lapply(plate2Dat, function(x)
  x[grepl("Unkn", Content),
  .(well = Well[1], copies = copies[1]),
  by = Content])

platePos[, ":="(row = toupper(sapply(strsplit(well, "-"), "[[", 1)),
  col = sapply(strsplit(well, "-"), "[[", 2))]

platePos[as.integer(col) < 10, col := paste0("0", col)]

platePos[, newWell := paste0(row, col)]

plate1Wells <- lapply(plate1Wells, function(x) x[, plate := 1])
plate2Wells <- lapply(plate2Wells, function(x) x[, plate := 2])

plate1Wells <- lapply(plate1Wells, function(x) x[, ":="(
  col = as.integer(substr(well, start = 2, stop = 3)),
  row = substr(well, start = 1, stop = 1))])

plate2Wells <- lapply(plate2Wells, function(x) x[, ":="(
  col = as.integer(substr(well, start = 2, stop = 3)),
  row = substr(well, start = 1, stop = 1))])

# get original well position
plate1Wells <- lapply(plate1Wells, function(x)
  x[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)])

plate1Wells <- lapply(plate1Wells, function(x)
  x[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)])

plate1Wells <- lapply(plate1Wells, function(x)
  x[origCol < 10, origWell := paste0(origRow, "0", origCol)])

plate1Wells <- lapply(plate1Wells, function(x)
  x[origCol >= 10, origWell := paste0(origRow, origCol)])

plate2Wells <- lapply(plate2Wells, function(x)
  x[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)])

plate2Wells <- lapply(plate2Wells, function(x)
  x[origCol < 10, origWell := paste0(origRow, "0", origCol)])

plate2Wells <- lapply(plate2Wells, function(x)
  x[origCol >= 10, origWell := paste0(origRow, origCol)])

allqPCR <- lapply(1:length(plate1Wells), function(x)
  rbindlist(list(plate1Wells[[x]], plate2Wells[[x]])))

qPCRDat <- lapply(1:length(geneList), function(x)
  platePos[allqPCR[[x]],
    geneList[x] := i.copies,
    on = c(newWell = "origWell", plate = "plate")])

qPCRDat[[1]][, (geneList) := lapply(.SD, function(x) round(x * 100, 0)), .SD = geneList]

fwrite(qPCRDat[[1]], "mergedqPCRData.csv")
