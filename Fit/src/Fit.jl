module Fit

using PlateReaderCore
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
        "--output_file", "-o"
            help = "Output file for slopes etc"
            required = true
        "--fit_file", "-f"
            help = "Output file for fitted values (predicted)"
        # "--plot_folder", "-p"
        #     help = "Folder for plots"
        #     required = true
        "--smoothing_parameter", "-s"
            help = "Smoothing parameter (1E-3)"
            default = 1E-3
            arg_type = Float64
        "--normalization", "-n"
            help = "Normalization scheme: unit_step, unit_range, none. Default: unit_step"
            default = "unit_step"
        
    end

    return parse_args(s)
end


function real_main()
    if false
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
    end


    parsed_args = parse_commandline()
    @info "Parsed arguments"
    @show parsed_args
    PlateReaderCore.app_fit(parsed_args)

    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end


end # module
