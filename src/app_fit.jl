## Main function for fitting app

if false
    args = Dict(
        "input_xlsxfile" => "dat_ex.xlsx",
        "fit_file" => "/tmp/TAPO_TEST/dat_ex_fit.xlsx",
        "plot_folder" => "/tmp/TAPO_TEST/",
        "smoothing_parameter" => 1E-3,
        "normalization" => "unit_step", ## "unit_step", "none", "unit_range"
    )
end

"""
    Fits based on .xlsx file with the following columns (* are mandatory):
    - * readerplate_id
    - * readerplate_well
    - * kinetic_time
    - * reader_value
    - reader_temperature
    - time_unit
    - value_unit
    - temperature_unit
    - readerplate_barcode
    - readerfile_name
    - readerplate_geometry
    - equipment
    - software
    - run_starttime    
"""
function app_fit(args)
    println("Input parameters:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    
    @info "Read input file"
    dat_df = PlateReaderCurves.xlsx(args["input_xlsxfile"]; sheet = 1)
    @info "Convert from data frame"
    dat = ReaderRun(dat_df)
    @info "Fit smoothing spline"
    ## select rescale
    x_range = xrange(dat)
    y_range = [0,1]
    if args["normalization"] == "unit_range"
        x_range = [0,1]
    elseif args["normalization"] == "none"
        x_range = missing
        y_range = missing
    end
    dat_fit = rc_fit(dat, "smooth_spline", lambda = args["smoothing_parameter"], x_range= x_range, y_range = y_range) 
    @info "Convert back to data frame"
    fit_df = DataFrame(dat_fit)
    @info "Save fit"
    mkpath(dirname(args["fit_file"]))
    XLSX.writetable(args["fit_file"], collect(DataFrames.eachcol(fit_df)), DataFrames.names(fit_df); overwrite=true, sheetname = "Fits")
    @info "Wrote " * args["fit_file"]
    @info "Plot plates"
    mkpath(dirname(args["plot_folder"]))
    fit_filenames = String[]
    phase_filenames = String[]
    Qs = dat.readerplate_geometry == 96 ? ["Q0"] : ["Q1","Q2","Q3","Q4"]
    for nPlate in 1:length(dat)
        for sQ in Qs 
            sPlate = "$(lpad(nPlate,3,'0'))_$(sQ)"
            @info sPlate
            push!(fit_filenames, joinpath(args["plot_folder"], "plate_$(sPlate)_fit.png"))
            plotObj = plateplot(Q(dat_fit.readerplates[nPlate], sQ))
            png(plotObj, fit_filenames[end])
            @info "did $(fit_filenames[end])"
            push!(phase_filenames, joinpath(args["plot_folder"], "plate_$(sPlate)_phase.png"))
            plotObj = plateplot(Q(dat_fit.readerplates[nPlate], sQ); type = "phase")
            png(plotObj, phase_filenames[end])
            @info "did $(phase_filenames[end])"
        end
    end
    htmlfile = joinpath(args["plot_folder"],"index.html")
    @info "Write $htmlfile"
    ## Use Hyperscript
    h_fitplot_names = [m("h2", basename(x)) for x in fit_filenames]
    h_fitplot = [m("img", src = x) for x in fit_filenames]
    h_pahseplot_names = [m("h2", basename(x)) for x in phase_filenames]
    h_phaseplot = [m("img", src = x) for x in phase_filenames]
    h_page = m("html",
               m("h1", "Plate plots")
               )(m("div").(zip(h_fitplot_names,h_fitplot,h_pahseplot_names,h_phaseplot)))
    savehtml(htmlfile, h_page)
end

Def_vals = Dict("reader_temperature" => missing,
                "time_unit" => "time",
                "value_unit" => "value",
                "temperature_unit" => "temperature",
                "readerplate_barcode" => "",
                "equipment" => "unknown",
                "software" => "unknown",
                "run_starttime" => missing,
                )

function fill_run_df!(df::DataFrame, def_vals; filename = "readerfile") ## add the missing columns
    def_vals = push!(copy(def_vals), ("readerfile_name" => filename))
    for (key, val) in def_vals
        if key ∉ names(df)
            df[:,key] .= val
        end
    end
    if "readerplate_geometry" ∉ names(df)
        df[:,"readerplate_geometry"] .= length(unique(df.readerplate_well))
    end
end

function fill_run_df(df::DataFrame, def_vals; filename = "readerfile")
    df2 = copy(df)
    fill_run_df!(df2, def_vals; filename = filename)
    df2
end
