var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PlateReaderCore","category":"page"},{"location":"#PlateReaderCore","page":"Home","title":"PlateReaderCore","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PlateReaderCore]","category":"page"},{"location":"#PlateReaderCore.ReaderCurve","page":"Home","title":"PlateReaderCore.ReaderCurve","text":"ReaderCurve: Datastructure for holding reader curves\nFields:\nreaderplate_well::String = \"well\"\nkinetic_time::Array\nreader_value::Array{Union{Missing, Real}}\nreader_temperature::Array{Union{Missing, Real}} = [missing]\ntime_unit::String\nvalue_unit::String\ntemperature_unit::String = \"C\"\n\n\n\n\n\n","category":"type"},{"location":"#PlateReaderCore.ReaderCurveFit","page":"Home","title":"PlateReaderCore.ReaderCurveFit","text":"ReaderCurveFit: Datastructure for holding reader curves and corresponding fits\nFields:\nreadercurve::ReaderCurve the input readercurve\nfit_method::String name of method to fit (linreg_trim, )\nfit_input_parameters::NamedTuple parameters given to fit method\npredict::Function fitted function. Can be used to predict new fitted values\nslope::Real max slope\nintercept::Real intercept of max slope curve\nfit_mean_absolute_residual::Real average absolute residuals of fit and read\n\n\n\n\n\n","category":"type"},{"location":"#PlateReaderCore.ReaderPlate","page":"Home","title":"PlateReaderCore.ReaderPlate","text":"ReaderPlate: Structure representing a readerplate\nreaderplate_id::String  globally unique eg from UUIDs.uuid4()\nreaderplate_barcode::String  can be \"\"\nreaderfile_name::String\nreaderplate_geometry::Int  96, 384\nreadercurves::Array{ReaderCurve,1} array of reader curves\n\n\n\n\n\n","category":"type"},{"location":"#PlateReaderCore.ReaderPlateFit","page":"Home","title":"PlateReaderCore.ReaderPlateFit","text":"ReaderPlateFit: Structure representing a fit of curves on a readerplate\n    Very similar to ReaderPlate\nreaderplate_id::String   globally unique eg from UUIDs.uuid4()\nreaderplate_barcode::String   can be \"\"\nreaderfile_name::String\nreaderplate_geometry::Int  96, 384\nreadercurves::Array{ReaderCurveFit}\n\n\n\n\n\n","category":"type"},{"location":"#PlateReaderCore.RelativeActivity","page":"Home","title":"PlateReaderCore.RelativeActivity","text":"struct RelativeActivity\nrelative_activity_id::String          Some name\nrelative_activity_value::Real         The relative activityvalue\ntest_activity::ReaderCurveFit         Input data\nreference_activity::ReaderCurveFit    Input data\ntest_activity_x::Real                 where it is measured\ntest_activity_y::Real                 where it is measured\nreference_activity_x::Real            where it is measured\nreference_activity_y::Real            where it is measured\nrelative_activity_method::String      how it was computed\n\n\n\n\n\n","category":"type"},{"location":"#PlateReaderCore.Q-Tuple{ReaderPlate, Any}","page":"Home","title":"PlateReaderCore.Q","text":"Q(::ReaderPlate, q; well96=false)\nQ(::ReaderPlateFit, q; well96=false)\nsubset a readerplate to a quadrant\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.app_fit-Tuple{Any}","page":"Home","title":"PlateReaderCore.app_fit","text":"Fits based on .xlsx file with the following columns (* are mandatory):\n- * readerplate_id\n- * readerplate_well\n- * kinetic_time\n- * reader_value\n- reader_temperature\n- time_unit\n- value_unit\n- temperature_unit\n- readerplate_barcode\n- readerfile_name\n- readerplate_geometry\n- equipment\n- software\n- run_starttime\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.linreg-Tuple{Any, Any}","page":"Home","title":"PlateReaderCore.linreg","text":"linreg(x, y): Linear regression\nOutput: (intercept, slope)\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.linreg_trim-Tuple{Any, Any}","page":"Home","title":"PlateReaderCore.linreg_trim","text":"linreg_trim(ReaderCurve; y_low_pct=0, y_high_pct): trimmed linear regression\nlinreg_trim(x,y; y_low_pct=0, y_high_pct)\nSkip the y_low_pct %, and y_high_pct % of the y-range.\nEg linreg_trim(x,y; 5,95) will use the central 90% of the y-range.\nNote it is using the range of the y-values. Not the number of values, as a quantile would do.\nOutput: (intercept, slope) or ReaderCurveFit object\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.plate-Tuple{PlateReaderCore.AbstractRun, Int64}","page":"Home","title":"PlateReaderCore.plate","text":"plate: get plate from a run based on number in run, id or barcode\nplate(r::AbstractRun, n::Int): get plate number n (starting from 1)\nplate(r::AbstractRun, name::String; search=[:readerplate_id, :readerplate_barcode]) search one or more of the fields readerplate_id, readerplate_barcode. Return first match. If partial==true, do a partial match.\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.rc_exp-NTuple{4, Any}","page":"Home","title":"PlateReaderCore.rc_exp","text":"Exponentially asymptotic readercurve\nrc_exp(t,A,k,y0) = y0 + A(1 - exp(-t/k))\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.rc_fit-Tuple{ReaderRun, String}","page":"Home","title":"PlateReaderCore.rc_fit","text":"rc_fit(::ReaderCurve, method::String; y_low_pct=10, y_high_pct=90, lambda = 250, l4p_parameter=100)\nrc_fit(::ReaderPlate, method::String; y_low_pct=10, y_high_pct=90, lambda = 250, l4p_parameter=100)\nFit a readercurve, plate of reader curves or a full run\nReturns a ReaderCurveFit containing the original readercurve and a predict function that can be used to predict new values. It also contains Slope, intercept and mean residual.\nSee @ref ReaderCurveFit\nMethods:\n- linreg_trim: linear regression omitting y_low_pct and y_high_pct of y range.\n- max_slope:\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.read_table-Tuple{String}","page":"Home","title":"PlateReaderCore.read_table","text":"read_table: read data from spread-sheet like format to DatraFrame\nSupports: .xlsx, .csv, .csv.gz, .csv2, .csv2.gz\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.well-Tuple{ReaderPlate, Array{String, N} where N}","page":"Home","title":"PlateReaderCore.well","text":"well(::ReaderPlate, well)::ReaderCurve\nwell(::ReaderPlateFit, well)::ReaderCurveFit\nwell(::ReaderPlate, well::Array{String})::Array{ReaderCurve}\nwell(::ReaderPlateFit, well::Array{String})::Array{ReaderCurveFit}\nselect one or more well(s) from a curve or a fit\n\n\n\n\n\n","category":"method"},{"location":"#PlateReaderCore.xrange-Tuple{ReaderCurve}","page":"Home","title":"PlateReaderCore.xrange","text":"xrange: range of kinetic_time over curve, plate or run\nReturns [xmin, xmax]\n\n\n\n\n\n","category":"method"}]
}