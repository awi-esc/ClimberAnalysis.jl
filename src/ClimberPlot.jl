module ClimberPlot

const MODEL_MEMBER_DELIM = "#"

include("plot-utils.jl")
include("plot-data.jl")

export loadData, loadPreprocData, getIndependenceWeights, getPerformanceWeights, getOverallWeights

end