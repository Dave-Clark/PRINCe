library(data.table)

### PLATE 2

plate <- fread("PRINCe_Bac_16S_plate2.csv")

pools <- plate[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

fwrite(pools, "PRINCe_Bac_16S_pool_plate2.csv")
# rows
# B = A
# E = B
# H = C
# K = D
# N = E

#1 1
#2 3
#3 5
#4 7

# PLATE 1
plate <- fread("PRINCe_Bac_16S_plate1.csv")

pools <- plate[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

fwrite(pools, "PRINCe_Bac_16S_pool_plate1.csv")

### hzo plate1 pooling
plate1  <- fread("PRINCe_hzo_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_hzo_pool_plate1.csv")

### hzo plate2 pooling

plate2 <- fread("PRINCe_hzo_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_hzo_pool_plate2.csv")

### nirS C1 plate1 pooling
plate1  <- fread("PRINCe_nirS_C1_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_nirS_C1_pool_plate1.csv")

### hzo plate2 pooling
plate2 <- fread("PRINCe_nirS_C1_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_nirS_C1_pool_plate2.csv")


### Archaea 16S plate1 pooling
plate1  <- fread("PRINCe_Arch_16S_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_arch_16S_pool_plate1.csv")

### hzo plate2 pooling
plate2 <- fread("PRINCe_Arch_16S_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_Arch_16S_pool_plate2.csv")

### COMAMMOX plate 1
plate1  <- fread("PRINCe_comammox_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_comammox_pool_plate1.csv")

### hzo plate2 pooling
plate2 <- fread("PRINCe_comammox_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]


# row A from plate 1
pools[col %% 2 == 0, ":="(origRow = "A",
  origCol = col / 2,
  plate = 1)]


pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(plate, origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_comammox_pool_plate2.csv")

## AOB pools
plate1  <- fread("PRINCe_AOB_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_AOB_pool_plate1.csv")

### hzo plate2 pooling
plate2 <- fread("PRINCe_AOB_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_AOB_pool_plate2.csv")

## AOA pools
plate1  <- fread("PRINCe_AOA_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_AOA_pool_plate1.csv")

### hzo plate2 pooling
plate2 <- fread("PRINCe_AOA_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_AOA_pool_plate2.csv")

## nxrB pools
plate1  <- fread("PRINCe_nxrB_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_nxrB_pool_plate1.csv")

### nxrB plate2 pooling
plate2 <- fread("PRINCe_nxrB_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_nxrB_pool_plate2.csv")

# nosZ Clade 1.

plate1  <- fread("PRINCe_nosZ_C1_plate1.csv")

pools <- plate1[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[col %% 2 != 0, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C", "D"))),
  origCol = (col + 1)/2)]

pools[col %% 2 == 0, ":="(origRow = ifelse(row == "B", "E",
  ifelse(row == "E", "F",
    ifelse(row == "H", "G", "H"))),
  origCol = col / 2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_nosZ_C1_pool_plate1.csv")

### nxrB plate2 pooling
plate2 <- fread("PRINCe_nosZ_C1_plate2.csv")

pools <- plate2[, .(row = row[1], col = col[1],
  quant = mean(ng_ul)), by = sample]

pools[quant <= 0, quant := 0.0001]

pools[, pool := round(max(quant)/quant, 2)]

# add original well position
pools[, ":="(origRow = ifelse(row == "B", "A",
  ifelse(row == "E", "B",
    ifelse(row == "H", "C",
      ifelse(row == "K", "D", "E")))),
  origCol = (col + 1)/2)]

pools[, c("sample", "row", "col") := NULL]

pools <- pools[order(origRow, origCol)]

setcolorder(pools, c("origRow", "origCol", "quant", "pool"))

fwrite(pools, "PRINCe_nosZ_C1_pool_plate2.csv")
