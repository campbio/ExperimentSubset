library(SingleCellExperiment)
.ExperimentSubset <- setClass("ExampleClass",
                          slots= representation(
                            rowVec="integer",
                            colVec="integer",
                            rowToRowMat="matrix",
                            colToColMat="matrix",
                            rowToColMat="matrix",
                            colToRowMat="matrix"
                          ),
                          contains="SingleCellExperiment"
)

ExperimentSubset <- function(
  rowVec=integer(0),
  colVec=integer(0),
  rowToRowMat=matrix(0,0,0),
  colToColMat=matrix(0,0,0),
  rowToColMat=matrix(0,0,0),
  colToRowMat=matrix(0,0,0),
  ...)
{
  se <- SingleCellExperiment(...)
  .ExperimentSubset(se, rowVec=rowVec, colVec=colVec,
                rowToRowMat=rowToRowMat, colToColMat=colToColMat,
                rowToColMat=rowToColMat, colToRowMat=colToRowMat)
}

