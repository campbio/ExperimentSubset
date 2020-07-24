ExperimentSubset <- setClass(
  Class = "ExperimentSubset",
  slots = c(
    rowIndex = "numeric",
    colIndex = "numeric"
  )
)