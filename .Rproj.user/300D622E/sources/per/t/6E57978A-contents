\name{REWISAV2}
\alias{REW-ISA V2}
\title{REW-ISA: RNA Expression Weighted Iterative Signature Algorithm}
\description{
REW-ISA is used to find potential local functional blocks (LFB) based on MeRIP-Seq data, where sites are hyper-methylated simultaneously across the specific conditions.
}
\usage{
### Given the range and step size of row threshold and column threshold, optimize the selection of thresholds in the above range.
result <- REWISAV2(FPKM.IP, FPKM.INPUT, Gene.index, Condition.cell, Optimize,
                      thr.IoU, thr.row.optimize, thr.col.optimize, max.epoch,
                      thr.row.start, thr.row.end, thr.row.step,
                      thr.col.start, thr.col.end, thr.col.step)

### Run REW-ISA under the selected optimal row and column threshold combination
REWISA_bicluster <- REWISAV2(FPKM.IP, FPKM.INPUT, Gene.index, Condition.cell,
                             optimize, thr.IoU, thr.row.optimize, thr.col.optimize,
                             max.epoch)
}
\arguments{
  \item{FPKM.IP}{IP sample data of m6A epi-transcriptome.}
  \item{FPKM.INPUT}{input sample data of m6A epi-transcriptome.}
  \item{Gene.index}{Vector. The genes corresponding to the m6A sites.}
  \item{Condition.cell}{Vector. The cell lines corresponding to the conditions.}
  \item{optimize}{Logical. If it is TRUE, start looking for the best threshold values.}
  \item{max.epoch}{Update the maximum number of iterations for a subset of rows and columns.}

  \item{fixed_side}{Logical. If TRUE, one of the row or column thresholds is fixed. In general, it is recommended to set it to FALSE.}
  \item{thr_col_interval}{Range of col threshold.}
  \item{thr.col.step}{The step size of the col threshold within its range.}
  \item{thr_row_interval}{Range of row threshold.}
  \item{thr.row.step}{The step size of the row threshold within its range.}
  \item{thr.row.optimize}{The optimal row threshold found by grid search.}
  \item{thr.col.optimize}{The optimal col threshold found by grid search.}
}
\value{
  \item{result}{List. Optimized row and col thresholds.}
  \item{REWISA_bicluster}{S4 class. A series of LFBs.}
}
\references{
Lin Zhang, Shutao Chen, Jiani Ma, Zhaoyang Liu, Hui Liu. REW-ISA V2: a biclustering method fusing homologous information for analyzing and mining epi-transcriptome data.
}
\author{
Shutao Chen <shutao.chen@cumt.edu.cn>, Lin Zhang, Jiani Ma, Zhaoyang Liu and Hui Liu.
}
\note{
Give a set of FPKM.IP, FPKM.INPUT to run REW-ISA
}
\examples{
### Load data
data(FPKM.IP, package = "REWISAV2")
data(FPKM.INPUT, package = "REWISAV2")

data <- as.matrix(FPKM.IP)
weight <- as.matrix(FPKM.INPUT)
### Using grid search to find the optimal row threshold, column threshold and the number of LFBs
thr.optimize <- REWISAV2(FPKM.IP = FPKM.IP, FPKM.INPUT = FPKM.INPUT,
                         Gene.index = Gene.index, Condition.cell = Condition.cell,
                         Optimize = TRUE, thr.IoU = 0.98, max.epoch = 30,
                         thr.row.start = 0.1, thr.col.start = 0.1,
                         thr.row.end = 3.0, thr.col.end = 3.0,
                         thr.row.step = 0.1, thr.col.step = 0.1)

### The final LFBs are determined according to the optimal row threshold, column threshold and the number of LFBs.
LFBs <- REWISAV2(FPKM.IP = FPKM.IP, FPKM.INPUT = FPKM.INPUT, Gene.index = Gene.index,
                 Condition.cell = Condition.cell, Optimize = FALSE, thr.IoU = 0.98,
                 max.epoch = 30, thr.row.optimize = 0.1, thr.col.optimize = 0.1)
}
