"""
    ReaderCurve: Datastructure for holding reader curves
    Fields:
    readerplate_well::String = "well"
    kinetic_time::Array
    reader_value::Array{Union{Missing, Real}}
    reader_temperature::Array{Union{Missing, Real}} = [missing]
    time_unit::String
    value_unit::String
    temperature_unit::String = "C"    
"""
Base.@kwdef struct ReaderCurve
    readerplate_well::String = "well"
    kinetic_time::Array{Real}
    reader_value::Array{Real} ## Use Inf, -Inf, NaN rather than Union{Missing, Real}}
    reader_temperature::Array{Union{Missing, Real}} = [missing]
    time_unit::String
    value_unit::String
    temperature_unit::String = "C"
end

function ReaderCurve(df::DataFrame)
    cols_needed(df, string.(fieldnames(ReaderCurve)), "ReaderCurve(::DataFrame)")
    # non_unique = @where(df, nonunique(df, :kinetic_time))
    # if(nrow(non_unique) > 0)
    #     println(non_unique)
    #     @error("Plate: $(unique(df.readerplate_id)), Well: $(unique(df.readerplate_well)). The above times are repeated")
    # end
#    @assert all(.!nonunique(df, :kinetic_time))
    @assert length(unique(df.readerplate_well)) == 1
    @assert length(unique(df.time_unit)) == 1
    @assert length(unique(df.value_unit)) == 1
    @assert length(unique(df.temperature_unit)) == 1
    ReaderCurve(readerplate_well = first(unique(df.readerplate_well)),
                kinetic_time = df.kinetic_time,
                reader_value = ifelse.(ismissing.(df.reader_value), NaN, df.reader_value), ## df.reader_value , ##
                reader_temperature = df.reader_temperature,
                time_unit = first(unique(df.time_unit)),
                value_unit = first(unique(df.value_unit)),
                temperature_unit = first(unique(df.temperature_unit))
                )
end

function cols_needed(df, cols, caller)
    missed_cols = setdiff(cols, names(df))
    if length(missed_cols) > 0
        error("""$caller is missing $(length(missed_cols)): $(join(missed_cols, ", "))""")
    end
    true
end

"""
    ReaderCurveFit: Datastructure for holding reader curves and corresponding fits
    Fields:
    readercurve::ReaderCurve the input readercurve
    fit_method::String name of method to fit (linreg_trim, )
    fit_input_parameters::NamedTuple parameters given to fit method
    predict::Function fitted function. Can be used to predict new fitted values
    slope::Real max slope
    intercept::Real intercept of max slope curve
    fit_mean_absolute_residual::Real average absolute residuals of fit and read
"""
Base.@kwdef struct ReaderCurveFit
    readercurve::ReaderCurve
    fit_method::String
    fit_input_parameters::NamedTuple
    predict::Function
    slope::Real
    intercept::Real
    inflectionpoint::Array{}
    fit_mean_absolute_residual::Real
end

abstract type AbstractPlate end

"""
    ReaderPlate: Structure representing a readerplate
    readerplate_id::String  globally unique eg from UUIDs.uuid4()
    readerplate_barcode::String  can be ""
    readerfile_name::String
    readerplate_geometry::Int  96, 384
    readercurves::Array{ReaderCurve} array of reader curves
"""
Base.@kwdef struct ReaderPlate <: AbstractPlate
    readerplate_id::String ## globally unique eg from UUIDs.uuid4()
    readerplate_barcode::String ## can be ""
    readerfile_name::String
    readerplate_geometry::Int ## 96, 384
    readercurves::Array{ReaderCurve}
end

"""
    ReaderPlateFit: Structure representing a fit of curves on a readerplate
        Very similar to ReaderPlate
    readerplate_id::String   globally unique eg from UUIDs.uuid4()
    readerplate_barcode::String   can be ""
    readerfile_name::String
    readerplate_geometry::Int  96, 384
    readercurves::Array{ReaderCurveFit}
"""
Base.@kwdef struct ReaderPlateFit <: AbstractPlate
    readerplate_id::String ## globally unique eg from UUIDs.uuid4()
    readerplate_barcode::String ## can be ""
    readerfile_name::String
    readerplate_geometry::Int ## 96, 384
    readercurves::Array{ReaderCurveFit}
end

function ReaderPlate(df::DataFrame)
    cols_needed(df, setdiff(string.(fieldnames(ReaderPlate)), ["readercurves"]), "ReaderPlate(::DataFrame)")
    @assert length(unique(df.readerplate_id)) == 1
    @assert length(unique(df.readerplate_barcode)) == 1
    @assert length(unique(df.readerfile_name)) == 1
    @assert length(unique(df.readerplate_geometry)) == 1
    curves = ReaderCurve[]
    for w in unique(df.readerplate_well)
        # push!(curves, ReaderCurve(df[df.readerplate_well .== w,:]))
        push!(curves, ReaderCurve(@where(df, :readerplate_well .== w)))
    end
    ReaderPlate(readerplate_id = first(unique(df.readerplate_id)),
                readerplate_barcode = first(unique(df.readerplate_barcode)),
                readerfile_name = first(unique(df.readerfile_name)),
                readerplate_geometry = first(unique(df.readerplate_geometry)),
                readercurves= curves)
end

Base.@kwdef struct ReaderFile
    equipment::String
    software::String
    run_starttime::Union{DateTime, Missing}
    readerplates::Array{ReaderPlate} 
end

abstract type AbstractRun end

Base.@kwdef struct ReaderRun <: AbstractRun
    equipment::String
    software::String
    run_starttime::Union{DateTime, Missing}
    readerplate_geometry::Int
    readerplates::Array{ReaderPlate} 
end

Base.@kwdef struct ReaderRunFit <: AbstractRun
    equipment::String
    software::String
    run_starttime::Union{DateTime, Missing}
    readerplate_geometry::Int
    readerplates::Array{ReaderPlateFit} 
end


function ReaderPlates(df::DataFrame)::Array{ReaderPlate}
    mycols = setdiff(string.(fieldnames(ReaderPlate)), ["readercurves"])
    cols_needed(df, mycols, "ReaderPlates(::DataFrame)")
    plates = ReaderPlate[]
    for p in unique(df.readerplate_id)
        push!(plates, ReaderPlate(@where(df, :readerplate_id .== p)))
    end
    plates
end

function ReaderRun(df::DataFrame)
    cols_needed(df, setdiff(string.(fieldnames(ReaderRun)), ["readerplates"]), "ReaderRun(::DataFrame)")
    @assert length(unique(df.equipment)) == 1
    @assert length(unique(df.software)) == 1
    @assert length(unique(df.run_starttime)) == 1
    @assert length(unique(df.readerplate_geometry)) == 1
    plates = ReaderPlates(df)
    ReaderRun(equipment = first(unique(df.equipment)),
              software = first(unique(df.software)),
              run_starttime = first(unique(df.run_starttime)),
              readerplate_geometry = first(unique(df.readerplate_geometry)),
              readerplates = plates)
end

function DataFrame(rcf::ReaderCurveFit)
    fits = DataFrame(fit_method = rcf.fit_method, slope = rcf.slope, intercept = rcf.intercept, fit_mean_absolute_residual = rcf.fit_mean_absolute_residual)
end

function DataFrame(rpf::ReaderPlateFit)
    wells = vcat([DataFrame(w) for w in rpf.readercurves]...)
    @transform(wells, readerplate_id = rpf.readerplate_id, readerplate_barcode = rpf.readerplate_barcode, readerfile_name = rpf.readerfile_name)
end

function DataFrame(rrf::ReaderRunFit)
    plates = vcat([DataFrame(p) for p in rrf.readerplates]...)
    @transform(plates, equipment = rrf.equipment, software = rrf.software, run_starttime = rrf.run_starttime, readerplate_geometry = rrf.readerplate_geometry)
end


function Base.length(p::AbstractPlate)
    length(p.readercurves)
end

Base.length(r::ReaderRun) = length(r.readerplates)

geometry(r::ReaderRun) = r.readerplate_geometry
geometry(p::ReaderPlate) = p.readerplate_geometry

well_names(p::ReaderPlateFit) =  map(x -> x.readercurve.readerplate_well, p.readercurves)
well_names(p::ReaderPlate) =  map(x -> x.readerplate_well, p.readercurves)

"""
    Q(::ReaderPlate, q; well96=false)
    Q(::ReaderPlateFit, q; well96=false)
    subset a readerplate to a quadrant
"""
function Q(p::ReaderPlate, q; well96=false)
    if (p.readerplate_geometry == 96) && (q == "Q0")
        return(p)
    end
    @assert occursin(r"^Q[1-4]$",q)
    @assert p.readerplate_geometry == 384
    sub_curves = filter(p.readercurves) do c
        MTP.Q(c.readerplate_well) == q
    end
    if well96
        sub_curves = map(sub_curves) do w
            Setfield.@set w.readerplate_well = MTP.well96(w.readerplate_well)
        end
    end
    geometry = well96 ? 96 : p.readerplate_geometry
    ReaderPlate(
        readerplate_id = p.readerplate_id,
        readerplate_barcode = p.readerplate_barcode,
        readerfile_name = p.readerfile_name,
        readerplate_geometry = geometry,
        readercurves = sub_curves
    )
end
function Q(p::ReaderPlateFit, q; well96=false)
    if (p.readerplate_geometry == 96) && (q == "Q0")
        return(p)
    end
    @assert occursin(r"^Q[1-4]$",q)
    @assert p.readerplate_geometry == 384
    sub_curves = filter(p.readercurves) do c
        MTP.Q(c.readercurve.readerplate_well) == q
    end
    if well96
        sub_curves = map(sub_curves) do w
            Setfield.@set w.readercurve.readerplate_well = MTP.well96(w.readercurve.readerplate_well)
        end
    end
    geometry = well96 ? 96 : p.readerplate_geometry
    ReaderPlateFit(
        readerplate_id = p.readerplate_id,
        readerplate_barcode = p.readerplate_barcode,
        readerfile_name = p.readerfile_name,
        readerplate_geometry = geometry,
        readercurves = sub_curves
    )
end
"""
    well(::ReaderPlate, well)::ReaderCurve
    well(::ReaderPlateFit, well)::ReaderCurveFit
    well(::ReaderPlate, well::Array{String})::Array{ReaderCurve}
    well(::ReaderPlateFit, well::Array{String})::Array{ReaderCurveFit}
    select one or more well(s) from a curve or a fit
"""
function well(p::ReaderPlate, wells::Array{String}) ## Selec wells
    filter(p.readercurves) do w
        w.readerplate_well ∈ wells
    end
end
function well(p::ReaderPlateFit, wells::Array{String}) ## Selec wells
    filter(p.readercurves) do w
        w.readercurve.readerplate_well ∈ wells
    end
end
function well(p::ReaderPlate, well::String) ## Selec a well
    filter(p.readercurves) do w
        w.readerplate_well == well
    end |> first
end
function well(p::ReaderPlateFit, well::String) ## Selec a well
    filter(p.readercurves) do w
        w.readercurve.readerplate_well == well
    end |> first
end

"""
    plate: get plate from a run based on number in run, id or barcode
    plate(r::AbstractRun, n::Int): get plate number n (starting from 1)
    plate(r::AbstractRun, name::String; search=[:readerplate_id, :readerplate_barcode]) search one or more of the fields readerplate_id, readerplate_barcode. Return first match. If partial==true, do a partial match.
"""
function plate(r::AbstractRun, n::Int)
    r.readerplates[n]
end
function plate(r::AbstractRun,name::String; search=[:readerplate_id, :readerplate_barcode], partial=false)
    for field in seach
        for p in r.readerplates
            if partial
                occursin(name, getproperty(p,field)) && return(p)
            end
            (name == getproperty(p,field)) && return(p)
        end
    end
    @info("Did not find plate matching $name")
    nothing
end

"""
    xrange: range of kinetic_time over curve, plate or run
    Returns [xmin, xmax]
"""
xrange(c::ReaderCurve) = [minimum(c.kinetic_time), maximum(c.kinetic_time)]    
function xrange(p::AbstractPlate)
    xmin = Inf
    xmax = -Inf
    for c in p.readercurves
        xmin = minimum([xmin, xrange(c)[1]])
        xmax = maximum([xmax, xrange(c)[2]])
    end
    [xmin, xmax]
end
function xrange(r::AbstractRun)
    xmin = Inf
    xmax = -Inf
    for p in r.readerplates
        xmin = minimum([xmin, xrange(p)[1]])
        xmax = maximum([xmax, xrange(p)[2]])
    end
    [xmin, xmax]
end

## io
"""
    read_table: read data from spread-sheet like format to DatraFrame
    Supports: .xlsx, .csv, .csv.gz, .csv2, .csv2.gz
"""
function read_table(file::String; sheet = 1)
    if endswith(file,".xlsx")
        return(xlsx(file; sheet = sheet))
    elseif (endswith(file, ".csv") || endswith(file, ".csv2"))
        return( DataFrame(CSV.File(file)) )
    elseif (endswith(file, ".csv.gz") || endswith(file, ".csv2.gz"))
        df = GZip.open(file, "r") do io
            CSV.File(io) |> DataFrame
        end
        return(df)
    else
        error("read_table dose not recognize filetype of $file")
    end
end
xlsx(file::String; sheet = 1) = DataFrame(XLSX.readtable(file, sheet)...)

xlsx_write(file::String, df::DataFrame) = XLSX.writetable(file, collect(DataFrames.eachcol(df)), DataFrames.names(df))
xlsx_write(file::String, r::ReaderRun) = xlsx_write(file, DataFrame(r))

## Relative Activity

"""
    struct RelativeActivity
    relative_activity_id::String          Some name
    relative_activity_value::Real         The relative activityvalue
    test_activity::ReaderCurveFit         Input data
    reference_activity::ReaderCurveFit    Input data
    test_activity_x::Real                 where it is measured
    test_activity_y::Real                 where it is measured
    reference_activity_x::Real            where it is measured
    reference_activity_y::Real            where it is measured
    relative_activity_method::String      how it was computed
"""
Base.@kwdef struct RelativeActivity
    relative_activity_id::String 
    relative_activity_value::Real
    test_activity::ReaderCurveFit
    reference_activity::ReaderCurveFit
#    test_activity_x::Real
#    test_activity_y::Real
#    reference_activity_x::Real
#    reference_activity_y::Real
    relative_activity_method::String
end

function RelativeActivity(t::ReaderCurveFit, r::ReaderCurveFit, method = "slope"; id=missing)
    @assert method ∈ ["slope"] ## comon_y
    ra_id = ismissing(id) ? "$t.readerplate_well/$r.readerplate_well" : string(id)
    if method == "slope"
        return(
            RelativeActivity(
                relative_activity_id = ra_id,
                relative_activity_value = t.slope / r.slope,
                test_activity = t,
                reference_activity = r,
                realtive_activity_method = method,
            )
        )
    else
        error("RelativeActivity: This should never happen!")
    end
end
