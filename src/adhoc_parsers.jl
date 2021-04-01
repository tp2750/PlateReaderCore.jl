using DataFramesMeta, Dates
## Parsers of local fileformats

Def_vals = Dict("reader_temperature" => missing,
                "time_unit" => "time",
                "value_unit" => "value",
                "temperature_unit" => "temperature",
                "readerplate_barcode" => "",
                "equipment" => "unknown",
                "software" => "unknown",
                "run_starttime" => missing,
                )


function tapo_v1(xlsx_file_path, value_unit, time_unit = "sec", temperature_unit = "C", sheet = 1,run_starttime = Dates.now(), equipment = "Neo2", software = "Gen5")
    df =  PlateReaderCore.xlsx(xlsx_file_path; sheet = sheet);
    ReaderRun(@transform(df,
                         equipment = equipment, software= software, run_starttime = run_starttime,
                         time_unit = time_unit, value_unit = value_unit, temperature_unit = temperature_unit,
                         readerplate_geometry = :geometry,
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
