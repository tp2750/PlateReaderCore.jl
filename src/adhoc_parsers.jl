using DataFramesMeta, Dates, Printf
## Parsers of local fileformats

Def_vals = Dict("reader_temperature" => missing,
                "time_unit" => "time",
                "value_unit" => "value",
                "temperature_unit" => "temperature",
                "readerplate_barcode" => "",
                "equipment" => "unknown",
                "software" => "unknown",
                "run_starttime" => missing,
                "experiment_id" => "",
                )


function tapo_v1(xlsx_file_path, value_unit, time_unit = "sec", temperature_unit = "C", sheet = 1,run_starttime = Dates.now(), equipment = "Neo2", software = "Gen5")
    df =  PlateReaderCore.xlsx(xlsx_file_path; sheet = sheet);
    ReaderRun(@transform(df,
                         equipment = equipment, software= software, run_starttime = run_starttime,
                         time_unit = time_unit, value_unit = value_unit, temperature_unit = temperature_unit,
                         readerplate_geometry = :geometry,
                         experiment_id = string.(:readerfile_name),
                         readerplate_id = :readerfile_name, readerplate_well = :well_name,
                         kinetic_time = :kinetic_sec,
                         reader_value = :absorbance_value, reader_temperature = :reader_temperature_C, 
    ))
end


function srlx_gp(xlsx_file_path, value_unit, time_unit = "sec", temperature_unit = "C", sheet = 1,run_starttime = Dates.now(), equipment = "Growth Profiler", software = "Growth Profiler")
    df = PlateReaderCore.xlsx(xlsx_file_path; sheet=sheet);
    ReaderRun(@transform(df,
                         equipment = equipment, software= software, run_starttime = run_starttime,
                         time_unit = time_unit, value_unit = value_unit, temperature_unit = temperature_unit,
                         readerplate_geometry = Int64.(:geometry),
                         experiment_id = string.(:readerfile_name),
                         readerplate_id = string.(:readerfile_name, :readerplate_number),
                         readerplate_well = string.(:well_name),
                         kinetic_time = Float64.(:kinetic_sec),
                         reader_value = Float64.(:absorbance_value), reader_temperature = :reader_temperature_C, 
    ))
end

function tsa2df(df::DataFrame; def_vals=Def_vals)
    df = @transform(df, readerplate_well = :well,
                    readerplate_id = string.(:exposure),
                    kinetic_time = Float64.(:T),
                    reader_value = Float64.(:Fluo)
                    )
    PlateReaderCore.fill_run_df!(df, def_vals)
    sort!(df, [:readerplate_id, :readerplate_well])
    ## ReaderRun(df)
    df
end
tsa2run(df::DataFrame; def_vals=Def_vals) = ReaderRun(tsa2df(df; def_vals = def_vals))

function vipr_trhs(data_file_path::String, platecount=6, timesteps=20, geometry=96; platenames = missing, time_unit = "ms", value_unit = "Pa", read_starttime = now()) 
    df =  PlateReaderCore.read_table(data_file_path)
    vipr_trhs(df, platecount, timesteps, geometry; experiment_id = basename(data_file_path), platenames = platenames, time_unit = time_unit, value_unit = value_unit, datafile = data_file_path,read_starttime = read_starttime)
end
function vipr_trhs(df::DataFrame, platecount=6, timesteps=20, geometry=96; experiment_id = "", platenames = missing, time_unit = "ms", value_unit = "Pa", datafile = "", read_starttime = now())
    @assert ncol(df)-1 == platecount * timesteps * geometry ## first column is curve time-point. Each column is a vipr-curve
    kinetic_points = df.Time
    @debug "Found $(length(kinetic_points)) curve-points."
    rows = Int64(sqrt(geometry/1.5))
    cols = Int64(sqrt(geometry*1.5))
    @assert geometry == rows * cols
    wells = PlateReaderCore.wells(rows,cols)
    if ismissing(platenames)
        platenames = [@sprintf("Plate%.2d", x) for x in 1:platecount]
    end
    @assert length(platenames) == platecount
    local readercurves = ReaderCurve[]
    local readerplates = ReaderPlate[]
    col_number = 1
    for plate in platenames
        for time_step in 1:timesteps
            for well in wells
                col_number += 1
                push!(readercurves, ReaderCurve(readerplate_well = well, kinetic_time = kinetic_points, reader_value = df[:,col_number], reader_temperature = repeat([missing], length(kinetic_points)), time_unit = time_unit, value_unit = value_unit, temperature_unit = "C"))
            end
            push!(readerplates, ReaderPlate(readerplate_id = @sprintf("%s_%.2d",plate, time_step), readerplate_barcode = "", readerfile_name = datafile, readerplate_geometry = geometry, readercurves = readercurves))
            readercurves = ReaderCurve[]
        end
    end
        ReaderRun(experiment_id = experiment_id, equipment="Vipr", software="", run_starttime=read_starttime, readerplate_geometry = geometry,readerplates= readerplates)
end
