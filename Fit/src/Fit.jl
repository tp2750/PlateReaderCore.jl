module Fit

using PlateReaderCurves
using DataFrames
using ArgParse

function julia_main()
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input_file", "-i"
            help = "Input file with all reads. Table in xlsx, csv or csv.gz format"
            required = true
        "--fit_file", "-f"
            help = "Output file for fitted values"
            required = true
        "--plot_folder", "-p"
            help = "Folder for plots"
            required = true
        "--smoothing_parameter", "-s"
            help = "Smoothing parameter (1E-3)"
            default = 1E-3
        "--normalization", "-n"
            help = "Normalization scheme: unit_step, unit_range, none. Default: unit_step"
            default = "unit_step"
        
    end

    return parse_args(s)
end


function real_main()
    @show ARGS
    @show Base.PROGRAM_FILE
    @show DEPOT_PATH
    @show LOAD_PATH
    @show pwd()
    @show Base.active_project()
    @show Threads.nthreads()
    @show Sys.BINDIR
#    display(Base.loaded_modules)
    @show unsafe_string(Base.JLOptions().image_file)
    println()

    @info "Running app!"

    parsed_args = parse_commandline()
    PlateReaderCurves.app_fit(parsed_args)

    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end


end # module
