module PlateReaderCore

using CSV, Statistics, SmoothingSplines, DataFrames, DataFramesMeta, Dates, Printf, LsqFit, XLSX, Hyperscript
using Distributions, Random, GZip
import Setfield

include("ReaderCurves.jl")
include("functions.jl")
include("app_fit.jl")

export ReaderCurve, ReaderCurveFit, ReaderPlate, ReaderPlateFit, ReaderFile, ReaderRun, geometry
export linreg_trim, smooth_spline, max_slope
export rc_fit, Q, well_names, well, well96
export xlsx

# Write your package code here.

end
