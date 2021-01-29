REWISAV2 <- function(FPKM.IP, FPKM.INPUT, Gene.index, Condition.cell, Optimize,
                      thr.IoU, thr.row.optimize, thr.col.optimize, max.epoch,
                      thr.row.start, thr.row.end, thr.row.step,
                      thr.col.start, thr.col.end, thr.col.step) {
  # Define subfunctions.
  # Global normalize.
  global_normalize <- function(data.process) {
    data.result <- (data.process - min(data.process)) / (max(data.process) - min(data.process))
    return(data.result)
  }

  # Row min-max normalize.
  row_normalize <- function(data.process) {
    gene.max <- apply(data.process, 1, max)
    gene.min <- apply(data.process, 1, min)
    gene.diff <- gene.max - gene.min
    data.row <- sweep(data.process, 1, gene.min, FUN = "-")
    data.row <- sweep(data.row, 1, gene.diff, FUN = "/")
    return(data.row)
  }

  # Column min-max normalize.
  col_normalize <- function(data.process) {
    condition.max <- apply(data.process, 2, max)
    condition.min <- apply(data.process, 2, min)
    condition.diff <- condition.max - condition.min
    data.col <- sweep(data.process, 2, condition.min, FUN = "-")
    data.col <- sweep(data.col, 2, condition.diff, FUN = "/")
    return(data.col)
  }

  # Eliminate global effects.
  eliminate_global_effect <- function(data.process) {
    data.process <- data.process - mean(data.process)
    gene.mean <- apply(data.process, 1, mean)
    condition.mean <- apply(data.process, 2, mean)
    data.process <- sweep(data.process, 1, gene.mean, FUN = "-")
    data.result <- sweep(data.process, 2, condition.mean, FUN = "-")
    return(data.result)
  }

  # Calculate the evaluation index WE_score.
  cal_we_score <- function(gene.vector) {
    gene.vector <- as.vector(unique(na.omit(gene.vector)))
    go <- enrichGO(gene.vector, 'org.Hs.eg.db', ont = "ALL", pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID', readable = TRUE)
    if (is.null(go)) {
      we.score <- 0
    } else {
      if (nrow(go@result) == 0) {
        we.score <- 0
      } else {
        result <- matrix(nrow = nrow(go@result), ncol = 2)
        result[, 1] <- go@result[, 10]
        result[, 2] <- as.numeric(go@result[, 6])

        gene.get <- go@result[, 9]
        gene.get <- paste(gene.get[1:length(gene.get)], collapse = "/")
        gene.use.num <- length(unique(unlist(strsplit(gene.get, "[/]"))))

        result[, 1] <- result[, 1] / length(gene.vector)
        result[, 2] <- -log(result[, 2], 10)
        denominator <- (apply(result, 2, sum))[1]
        denominator <- denominator + ((length(gene.vector) - gene.use.num) / length(gene.vector))
        molecule <- sum(result[, 1] * result[, 2])
        we.score <- molecule / denominator
      }
    }
    return(we.score)
  }

  # Data preprocessing.
  FPKM.IP <- apply(FPKM.IP, 2, as.numeric)
  FPKM.INPUT <- apply(FPKM.INPUT, 2, as.numeric)
  data.sum <- FPKM.IP + FPKM.INPUT
  Gene.index.unique <- unique(na.omit(Gene.index))

  # Generate methylation level matrix.
  FPKM.IP <- FPKM.IP + 0.001
  FPKM.INPUT <- FPKM.INPUT + 0.001
  data.sum <- data.sum + 0.002
  Methylation.level <- FPKM.IP / data.sum
  Methylation.rec <- Methylation.level
  len.row <- nrow(Methylation.level)
  len.col <- ncol(Methylation.level)

  # Generate expression level matrix
  Expression.level <- log2(data.sum + 1 - 0.002)
  Expression.level <- global_normalize(Expression.level)
  # Eliminate global effects
  Methylation.level <- eliminate_global_effect(Methylation.level)
  # Normalization of methylation level matrix and expression level matrix
  Methylation.level <- global_normalize(Methylation.level)
  # Expression.level <- global_normalize(Expression.level)
  # Row min-max normalize
  Methylation.row <- row_normalize(Methylation.level)
  # Column min-max normalize
  Methylation.col <- col_normalize(Methylation.level)
  # Record data for similarity calculation
  Expression.rec <- Expression.level

  if (Optimize == FALSE) {
    thr.row.start <- thr.row.optimize
    thr.col.start <- thr.col.optimize
    thr.row.end <- thr.row.optimize
    thr.col.end <- thr.col.optimize
    thr.row.step <- 2 * thr.row.optimize
    thr.col.step <- 2 * thr.col.optimize
  }

  # Preset hyperparameters
  seed.num <- round(nrow(Methylation.rec))
  thr.IoU <- thr.IoU
  thr.row.start <- thr.row.start
  thr.row.end <- thr.row.end
  thr.row.step <- thr.row.step
  thr.col.start <- thr.col.start
  thr.col.end <- thr.col.end
  thr.col.step <- thr.col.step

  thr.row.loc <- 0
  bicluster.rec <- data.frame()
  we.final <- data.frame()
  we.score.original <- data.frame()
  sites.num.rec <- data.frame()
  genes.num.rec <- data.frame()
  # random.score.rec <- data.frame()
  random.genes.result <- data.frame()
  long.index.row <- 1

  # Main function
  for (thr.row in seq(thr.row.start, thr.row.end, thr.row.step)) {
    thr.col.loc <- 0
    thr.row.loc <- thr.row.loc + 1
    for (thr.col in seq(thr.col.start, thr.col.end, thr.col.step)) {
      cat("TR:", thr.row, "\t")
      cat("TC:", thr.col, "\n")
      thr.col.loc <- thr.col.loc + 1
      find.bicluster <- TRUE
      bicluster.num <- 0
      bicluster.mean <- mean(Methylation.rec * Expression.rec) + 1
      RowxNumber <- data.frame()
      NumberxCol <- data.frame()

      while (bicluster.mean >= mean(Methylation.rec * Expression.rec) && find.bicluster) {
        cat("Current bicluster ID:", bicluster.num + 1, "\n")
        if (bicluster.num != 0) {
          if (bicluster.num == 1) {
            data.temp <- Methylation.level
            # weight.temp <- Expression.level
          }
          data.temp[row.indicator, col.indicator] <- data.temp[row.indicator, col.indicator] - mean(data.temp[row.indicator, col.indicator])
          data.temp[data.temp < 0] <- 0
          data.temp.row <- row_normalize(data.temp)
          data.temp.col <- col_normalize(data.temp)
          weight.temp[row.indicator, col.indicator] <- weight.temp[row.indicator, col.indicator] - mean(weight.temp[row.indicator, col.indicator])
          weight.temp[weight.temp < 0] <- 0
        } else {
          data.temp.row <- Methylation.row
          data.temp.col <- Methylation.col
          weight.temp <- Expression.level
        }

        epoch <- 1
        area.IoU <- -1
        find.bicluster <- FALSE
        IoU.rec <- vector()
        repeat.distance <- vector()
        repeat.count <- 0
        IoU.index.diff <- 0
        while (area.IoU <= thr.IoU) {
          # Set the maximum number of search rounds.
          if (epoch >= max.epoch) {
            cat("The number of searches reached the upper limit, so quit the search!", "\n")
            find.bicluster <- FALSE
            break
          }
          cat("Current epoch:", epoch, "\t")
          time.start <- proc.time()[[1]]

          # initial row indicator.
          if (epoch == 1) {
            row.indicator <- sort(sample(nrow(Methylation.level), seed.num, FALSE))
            col.indicator <- seq(1, ncol(Methylation.level), 1)
          }

          # select col based on row.indicator.
          if (length(row.indicator) == 0) {
            cat("The row indicator values are all zero, exit the search!", "\n")
            find.bicluster <- FALSE
            break
          } else {
            if (length(row.indicator) == 1) {
              data.min <- t(as.matrix(data.temp.row[row.indicator, ] * weight.temp[row.indicator, ]))
              data.pearson <- t(as.matrix(Methylation.rec[row.indicator, ] * Expression.rec[row.indicator, ]))
            } else {
              data.min <- as.matrix(data.temp.row[row.indicator, ] * weight.temp[row.indicator, ])
              data.pearson <- as.matrix(Methylation.rec[row.indicator, ] * Expression.rec[row.indicator, ])
            }
            # Calculate the Pearson correlation coefficient of the column.
            Condition.pearson <- vector()
            for (var in 1:ncol(Methylation.level)) {
              col.select <- which(Condition.cell == Condition.cell[var])
              if (length(col.select) != 1) {
                if (length(row.indicator) == 1) {
                  condition.mean <- apply(t(as.matrix(data.pearson[1:nrow(data.pearson), col.select])), 1, mean)
                } else {
                  condition.mean <- apply(as.matrix(data.pearson[1:nrow(data.pearson), col.select]), 1, mean)
                }
                if (length(condition.mean) == 1) {
                  Condition.pearson[var] <- 1 - abs(data.pearson[1:nrow(data.pearson), var] - condition.mean)
                } else {
                  if (sd(data.pearson[1:nrow(data.pearson), var]) != 0 && sd(condition.mean) != 0) {
                    Condition.pearson[var] <- abs(cor(data.pearson[1:nrow(data.pearson), var],
                                                      condition.mean, method = "pearson"))
                  } else {
                    Condition.pearson[var] <- 0
                  }
                }
              }
            }
            Condition.pearson[which(is.na(Condition.pearson))] <- mean(na.omit(Condition.pearson))
            col.score <- Condition.pearson * apply(data.min, 2, mean)
            col.score <- col.score - mean(col.score)
            col.indicator.last <- col.indicator
            col.indicator <- which(col.score >= (thr.col / sqrt(length(row.indicator))))
          }


          # select row based on col.indicator.
          if (length(col.indicator) == 0) {
            cat("The col indicator values are all zero, exit the search!", "\n")
            find.bicluster <- FALSE
            break
          } else {
            data.min <- as.matrix(data.temp.col[, col.indicator] * weight.temp[, col.indicator])
            data.pearson <- as.matrix(Methylation.rec[, col.indicator] * Expression.rec[, col.indicator])
            # Calculate the Pearson correlation coefficient of the row.
            Gene.pearson <- vector()
            for (var in 1:nrow(Methylation.level)) {
              row.select <- which(Gene.index == Gene.index[var])
              if (length(row.select) != 1) {
                gene.mean <- apply(as.matrix(data.pearson[row.select, 1:ncol(data.pearson)]), 2, mean)
                if (length(gene.mean) == 1) {
                  Gene.pearson[var] <- 1 - abs(data.pearson[var, 1:ncol(data.pearson)] - gene.mean)
                } else {
                  if (sd(data.pearson[var, 1:ncol(data.pearson)]) != 0 && sd(gene.mean) != 0) {
                    Gene.pearson[var] <- abs(cor(data.pearson[var, 1:ncol(data.pearson)],
                                                 gene.mean, method = "pearson"))
                  } else {
                    Gene.pearson[var] <- 0
                  }
                }
              }
            }
            Gene.pearson[which(is.na(Gene.pearson))] <- mean(na.omit(Gene.pearson))
            row.score <- Gene.pearson * apply(data.min, 1, mean)
            row.score <- row.score - mean(row.score)
            row.indicator.last <- row.indicator
            row.indicator <- which(row.score >= (thr.row / sqrt(length(col.indicator))))
          }

          # Calculate the IoU between the previous round and the current round.
          inter.row <- length(intersect(row.indicator.last, row.indicator))
          inter.col <- length(intersect(col.indicator.last, col.indicator))
          inter.area <- inter.row * inter.col
          area.last <- length(row.indicator.last) * length(col.indicator.last)
          area.current <- length(row.indicator) * length(col.indicator)
          union.area <- area.last + area.current - inter.area
          area.IoU <- inter.area / union.area
          epoch <- epoch + 1
          cat("Area IoU:", sprintf("%0.4f", area.IoU), "\t")
          if (area.IoU >= thr.IoU) {
            find.bicluster <- TRUE
          }

          if (length(IoU.rec) != 0) {
            IoU.index <- which(IoU.rec == area.IoU)
            if (length(IoU.index) != 0 && area.IoU != 0) {
              if (length(IoU.index) == 1) {
                if (IoU.index == length(IoU.rec)) {
                  cat("IoU enters an endless cycle!", "\t")
                  break
                } else {
                  if (repeat.count == 0) {
                    IoU.index.last <- IoU.index
                  } else {
                    IoU.index.diff <- IoU.index - IoU.index.last
                  }

                  IoU.rec[length(IoU.rec) + 1] <- area.IoU
                  repeat.count <- repeat.count + 1
                  repeat.distance[repeat.count] <- length(IoU.rec) - IoU.index
                  if (length(repeat.distance) == 2) {
                    if (repeat.distance[1] == repeat.distance[2] && IoU.index.diff == 1) {
                      cat("IoU falls into interval endless cycle!", "\t")
                      break
                    } else {
                      repeat.count <- 1
                      IoU.index.last <- IoU.index
                      IoU.index.diff <- 0
                      repeat.distance <- vector()
                      repeat.distance[repeat.count] <- length(IoU.rec) - IoU.index
                    }
                  }
                }
              } else {
                iou.rec.add <- length(IoU.rec) + 1
                if ((iou.rec.add - IoU.index[2]) == (IoU.index[2] - IoU.index[1])) {
                  cat("IoU falls into interval endless cycle!", "\t")
                  break
                }
              }
            } else {
              IoU.rec[length(IoU.rec) + 1] <- area.IoU
            }
          } else {
            IoU.rec[length(IoU.rec) + 1] <- area.IoU
          }




          # Calculate the running time of FBCwPlaid to find the bicluster.
          time.end <- proc.time()[[1]]
          time.use <- as.character(sprintf("%0.4f", time.end - time.start))
          time.use <- paste(time.use, " s", sep = "")
          cat("Running time:", time.use, "\n")
        }

        # Store information about the found bicluster.
        if (find.bicluster) {
          bicluster.num <- bicluster.num + 1
          bicluster.mean <- mean(Methylation.rec[row.indicator, col.indicator] * Expression.rec[row.indicator, col.indicator])
          RowxNumber[1:nrow(Methylation.rec), bicluster.num] <- FALSE
          NumberxCol[bicluster.num, 1:ncol(Methylation.rec)] <- FALSE
          RowxNumber[row.indicator, bicluster.num] <- TRUE
          NumberxCol[bicluster.num, col.indicator] <- TRUE
        }
      }

      RowxNumber.sum <- apply(RowxNumber, 2, sum)
      RowxNumber.drop <- which(RowxNumber.sum == 1)
      if (length(RowxNumber.drop) != 0) {
        RowxNumber <- RowxNumber[, -RowxNumber.drop]
        NumberxCol <- NumberxCol[-RowxNumber.drop, ]
        bicluster.num <- bicluster.num - length(RowxNumber.drop)
      }

      # Record the number of bicluster
      bicluster.rec[thr.row.loc, thr.col.loc] <- bicluster.num

      if (Optimize == TRUE) {
        # Calculate the relative growth rate of WE_score.
        if (bicluster.num > 0) {
          we.get <- vector()
          we.random <- vector()
          score.temp <- 0
          post <- 1

          gene_all <- data.frame()
          for (bicluster.var in 1:bicluster.num) {
            gene.temp <- Gene.index[which(RowxNumber[, bicluster.var] == TRUE)]
            sites.num.rec[long.index.row, bicluster.var] <- length(gene.temp)
            gene_all[1:length(gene.temp), bicluster.var] <- gene.temp
          }

          for (var in 1:ncol(gene_all)) {
            gene_data <- gene_all[, var]
            gene_random_len <- length(unique(na.omit(gene_data)))
            we.get[var] <- cal_we_score(gene_data) / gene_random_len
            we.score.original[long.index.row, var] <- we.get[var] * gene_random_len
            genes.num.rec[long.index.row, var] <- gene_random_len

            random.storage <- FALSE
            if (nrow(random.genes.result) == 0) {
              random.genes.result[1, 1] <- gene_random_len
              random.index <- 1
            } else {
              if (length(which(random.genes.result[, 1] == gene_random_len)) == 0) {
                random.genes.result[(nrow(random.genes.result) + 1), 1] <- gene_random_len
                random.index <- nrow(random.genes.result)
              } else {
                random.index <- which(random.genes.result[, 1] == gene_random_len)
                random.storage <- TRUE
              }
            }

            if (!random.storage) {
              we.random.temp <- vector()
              for (random.cycle in 1:8) {
                gene.sample.index <- sort(sample(length(Gene.index.unique), gene_random_len, FALSE))
                gene.sample <- as.vector(unique(na.omit(Gene.index.unique[gene.sample.index])))
                we.random.temp[random.cycle] <- cal_we_score(gene.sample) / gene_random_len
              }
              # we.random.temp <- sort(we.random.temp)[-1]
              # we.random.temp <- we.random.temp[-length(we.random.temp)]
              # Avoid Inf in calculation results.
              while(mean(we.random.temp) == 0) {
                gene.sample.index <- sort(sample(length(Gene.index.unique), gene_random_len, FALSE))
                gene.sample <- as.vector(unique(na.omit(Gene.index.unique[gene.sample.index])))
                we.random.temp[1] <- cal_we_score(gene.sample) / gene_random_len
              }
              we.random[var] <- mean(we.random.temp)
              random.genes.result[random.index, 2] <- mean(we.random.temp)
            } else {
              cat("Successful use of information in the random.genes.result!", "\n")
              we.random[var] <- random.genes.result[random.index, 2]
            }
            # random.score.rec[long.index.row, var] <- we.random[var]
          }
          we.final[thr.row.loc, thr.col.loc] <- mean((we.get - we.random) / we.random)
          cat("The mean of WE_score:", sprintf("%0.4f", we.final[thr.row.loc, thr.col.loc]), "\n")
        } else {
          we.final[thr.row.loc, thr.col.loc] <- 0
          we.score.original[long.index.row, 1] <- 0
          sites.num.rec[long.index.row, 1] <- 0
          genes.num.rec[long.index.row, 1] <- 0
          # random.score.rec[long.index.row, 1] <- 0
          cat("The mean of WE_score:", sprintf("%0.4f", 0), "\n")
        }
      }
      long.index.row <- long.index.row + 1
    }
  }

  if (Optimize == TRUE) {
    row.index <- which(bicluster.rec == max(bicluster.rec), arr.ind = T)[1, 1]
    col.index <- which(bicluster.rec == max(bicluster.rec), arr.ind = T)[1, 2]
    thr.row <- thr.row.start + (row.index - 1) * thr.row.step
    thr.col <- thr.col.start + (col.index - 1) * thr.col.step
    out.result <- list(the.row.optimize = thr.row, thr.col.optimize = thr.col)
    return(out.result)
  } else {
    setClass("Bi_clust", slots = list(RowxNumber = "matrix", NumberxCol = "matrix",
                                      Number = "numeric"))
    RowxNumber <- as.matrix(RowxNumber)
    NumberxCol <- as.matrix(NumberxCol)
    bc <- new("Bi_clust", RowxNumber = RowxNumber, NumberxCol = NumberxCol,
              Number = bicluster.num)
    return(bc)
  }
}

