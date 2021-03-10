""" 
    rc_fit(::ReaderCurve, method::String; y_low_pct=10, y_high_pct=90, lambda = 250, l4p_parameter=100)
    rc_fit(::ReaderPlate, method::String; y_low_pct=10, y_high_pct=90, lambda = 250, l4p_parameter=100)
    Fit a readercurve, plate of reader curves or a full run
    Returns a ReaderCurveFit containing the original readercurve and a predict function that can be used to predict new values. It also contains Slope, intercept and mean residual.
    See @ref ReaderCurveFit
    Methods:
    - linreg_trim: linear regression omitting y_low_pct and y_high_pct of y range.
    - max_slope:
"""
function rc_fit(r::ReaderRun, method::String; y_low_pct=10, y_high_pct=90, lambda = 1E-6, l4p_parameter=100, x_range=missing, y_range=missing )
    plate_fits = map(r.readerplates) do p
        rc_fit(p, method;y_low_pct, y_high_pct, lambda, l4p_parameter,  x_range=missing, y_range=missing)
    end
    ReaderRunFit(
        equipment = r.equipment,
        software = r.software,
        run_starttime = r.run_starttime,
        readerplate_geometry = r. readerplate_geometry,
        readerplates = plate_fits
    )
end
function rc_fit(p::ReaderPlate, method::String; y_low_pct=10, y_high_pct=90, lambda = 1E-6, l4p_parameter=100,  x_range=missing, y_range=missing)
    curve_fits = map(p.readercurves) do rc
        rc_fit(rc, method; y_low_pct, y_high_pct, lambda, l4p_parameter,  x_range=missing, y_range=missing)
    end
    ReaderPlateFit(
        readerplate_id = p.readerplate_id,
        readerplate_barcode = p.readerplate_barcode,
        readerfile_name = p.readerfile_name,
        readerplate_geometry = p.readerplate_geometry,
        readercurves = curve_fits
    )
end
function rc_fit(rc::ReaderCurve, method::String; y_low_pct=10, y_high_pct=90, lambda = 1E-6, l4p_parameter=100, x_range=missing, y_range=missing)
    ## method dispatch options: https://discourse.julialang.org/t/dispatch-and-symbols/21162/7?u=tp2750
    (X,Y) = get_finite(rc.kinetic_time, rc.reader_value)
    if(length(Y) == 0)
        return(
            ReaderCurveFit(
                readercurve = rc,
                fit_method = method,
                fit_input_parameters = (;y_low_pct, y_high_pct, lambda, x_range, y_range),
                predict = t -> NaN,
                slope = NaN,
                intercept = NaN,
                inflectionpoint = [NaN,NaN],
                fit_mean_absolute_residual = NaN
            )
        )
    end
    if method == "linreg_trim"
        f1 = linreg_trim(X,Y; y_low_pct, y_high_pct=90)
        pred_fun1(t) = f1.intercept + f1.slope * t
        residuals = abs.(pred_fun1.(X) .- Y)
        min_res_idx = findfirst(residuals .== minimum(residuals))
        return(
            ReaderCurveFit(
                readercurve = rc,
                fit_method = method,
                fit_input_parameters = (;y_low_pct, y_high_pct),
                predict = pred_fun1,
                slope = f1.slope,
                intercept = f1.intercept,
                inflectionpoint = [X[min_res_idx], Y[min_res_idx]],
                fit_mean_absolute_residual = mean(abs.(pred_fun1.(X) .- Y))
            )
        )
    elseif (method == "max_slope")
        f1 = max_slope(X,Y)
        pred_fun2(t) = f1.intercept + f1.slope * t
        return(
            ReaderCurveFit(
                readercurve = rc,
                fit_method = method,
                fit_input_parameters = (;),
                predict = pred_fun2,
                slope = f1.slope,
                intercept = f1.intercept,
                inflectionpoint = f1.inflectionpoint,
                fit_mean_absolute_residual = mean(abs.(pred_fun2.(X) .- Y))
            )
        )
    elseif method == "smooth_spline_old"
        l1 = convert(Float64,lambda)
        X_range = ismissing(x_range) ? [minimum(X), maximum(X)] : x_range
        Y_range = ismissing(y_range) ? [minimum(Y), maximum(Y)] : y_range
        X1 = PlateReaderCore.scale_fwd.(X;x_range = X_range)
        Y1 = PlateReaderCore.scale_fwd.(Y;x_range = Y_range)
        # f1 = smooth_spline_fit(X,Y; lambda = l1)
        # pred_fun3(t) = SmoothingSplines.predict(f1,convert(Float64,t))
        f1 = smooth_spline_fit(X1,Y1; lambda = l1)
        pred_fun3(t) = PlateReaderCore.scale_rev.(SmoothingSplines.predict(f1,PlateReaderCore.scale_fwd.(convert(Float64,t);x_range= X_range)), ;x_range= Y_range)
        ms = max_slope(X,pred_fun3.(X))
        return(
            ReaderCurveFit(
                readercurve = rc,
                fit_method = method,
                fit_input_parameters = (;),
                predict = pred_fun3,
                slope = ms.slope,
                intercept = ms.intercept,
                inflectionpoint = ms.inflectionpoint,
                fit_mean_absolute_residual = mean(abs.(pred_fun3.(X) .- Y))
            )
        )
    elseif method == "smooth_spline"
        l1 = convert(Float64,lambda)
        X_range_in = extrema(X)
        Y_range_in = extrema(Y)
        X_range = ismissing(x_range) ? [minimum(X), maximum(X)] : x_range
        Y_range = ismissing(y_range) ? [minimum(Y), maximum(Y)] : y_range
        X1 = PlateReaderCore.scale_fwd(X, X_range_in, X_range)
        Y1 = PlateReaderCore.scale_fwd(Y, Y_range_in, Y_range)
        # f1 = smooth_spline_fit(X,Y; lambda = l1)
        # pred_fun3(t) = SmoothingSplines.predict(f1,convert(Float64,t))
        f2 = smooth_spline_fit(X1,Y1; lambda = l1)
        pred_fun4(t) = PlateReaderCore.scale_fwd(SmoothingSplines.predict(f2,PlateReaderCore.scale_fwd(convert(Float64,t),X_range_in, X_range)), Y_range, Y_range_in) ## TODO update scale(x,r1,r2)
        ms = max_slope(X,pred_fun4.(X))
        return(
            ReaderCurveFit(
                readercurve = rc,
                fit_method = method,
                fit_input_parameters = (;),
                predict = pred_fun4,
                slope = ms.slope,
                intercept = ms.intercept,
                inflectionpoint = ms.inflectionpoint,
                fit_mean_absolute_residual = mean(abs.(pred_fun4.(X) .- Y))
            )
        )
    elseif method == "exp"
        p0 = [0,1,1,1]
        f1 = LsqFit.curve_fit(rc_exp,rc.kinetic_time, rc.reader_value, p0)
    elseif method == "L4P"
        lc4_params = rc_logistic_fit(X,Y; l4p_parameter=100)
        pred_l4p = t -> rc_logistic_fun(t,lc4_params)
        ms = max_slope(X,pred_l4p.(X))
        return(
            ReaderCurveFit(
                readercurve = rc,
                fit_method = method,
                fit_input_parameters = (;l4p_parameter),
                predict = pred_l4p,
                slope = ms.slope,
                intercept = ms.intercept,
                inflectionpoint = ms.inflectionpoint,
                fit_mean_absolute_residual = mean(abs.(pred_l4p.(X) .- Y))
            )
        )        
    else
        error("This should not happen")
    end
end

"""
    Exponentially asymptotic readercurve
    rc_exp(t,A,k,y0) = y0 + A(1 - exp(-t/k))
"""
rc_exp(t,A,k,y0) = y0 .+ A.*(1 .- exp.(-t ./k))

rc_exp(t,p) = p[3] .+ p[1].*(1 .- exp.(-t ./p[2]))

function smooth_spline_fit(x,y; lambda=250.0)
    ## TODO auto select lambda based on GCV. see paper: ../Lukas_deH_RSA_2015.pdf
    ## From https://github.com/nignatiadis/SmoothingSplines.jl
    X = map(Float64,convert(Array,x))
    Y = map(Float64,convert(Array,y))
    spl = SmoothingSplines.fit(SmoothingSpline, X, Y, lambda)
    # Ypred = SmoothingSplines.predict(spl)
    # Ypred
end


"""
    linreg(x, y): Linear regression
    Output: (intercept, slope)
"""
function linreg(x, y)
    (i,s) = hcat(fill!(similar(x), 1), x) \ y ## https://github.com/JuliaStats/StatsBase.jl/issues/398#issuecomment-417875619
    (intercept = i, slope = s)
end

"""
    linreg_trim(ReaderCurve; y_low_pct=0, y_high_pct): trimmed linear regression
    linreg_trim(x,y; y_low_pct=0, y_high_pct)
    Skip the y_low_pct %, and y_high_pct % of the y-range.
    Eg linreg_trim(x,y; 5,95) will use the central 90% of the y-range.
    Note it is using the range of the y-values. Not the number of values, as a quantile would do.
    Output: (intercept, slope) or ReaderCurveFit object
"""
function linreg_trim(x,y; y_low_pct=10, y_high_pct=90)
    X,Y = get_finite(x,y)
    if(length(Y) == 0)
        return(intercept = NaN, slope = NaN)
    end
    y_cover = maximum(Y) - minimum(Y)
    y_1 = minimum(Y) + y_low_pct /100 * y_cover
    y_2 = minimum(Y) + y_high_pct/100 * y_cover
    idx = (Y .>= y_1) .& (Y .<= y_2)
    linreg(X[idx], Y[idx])
end

function max_slope(x,y)
    X,Y = get_finite(x,y)
    if(length(Y) == 0)
        return(intercept = NaN, slope = NaN, inflectionpoint = NaN)
    end
    slopes = diff(Y) ./ diff(X)
    slope = maximum(slopes)
    slope_idx = findfirst(slopes .== slope)    
    b = y[slope_idx] - slope * x[slope_idx]
    (intercept = b,slope = slope, inflectionpoint = [x[slope_idx], y[slope_idx]])
end

function get_finite(x,y)
    y_finite = isfinite.(y)
    X = x[y_finite]
    Y = y[y_finite]
    (X,Y)
end    

function rc_logistic_fun(t,p)
    @. p[1]/(1+exp(4*p[2]/p[1]*(p[3]-t)+2))+p[4] # logistic
end


function rc_logistic_fit(x,y; l4p_parameter=100)
    guess_parameter = l4p_parameter
    p0 = [0,0,0,0] ## A, µ, λ, bg
    fit1 = linreg_trim(x,y; y_low_pct=10, y_high_pct=90)
    A = maximum(y)
    µ = fit1[2]
    λ = -fit1[1]/fit1[2]
    bg = minimum(y)
    p0 = [A, µ, λ, bg]
    pl = minimum.([1/guess_parameter,guess_parameter, -1/guess_parameter, -guess_parameter].*x for x in [abs(A)+abs(bg), µ, λ, abs(A)+abs(bg)])
    pu = maximum.([1/guess_parameter,guess_parameter, -1/guess_parameter, -guess_parameter].*x for x in [abs(A)+abs(bg), µ, λ, abs(A)+abs(bg)])  ## [10 *A, 10 *µ, 10 *λ, maximum([bg*10, -bg*10])] [0.01,100, -0.01, -100]
    @debug p0 - pl
    @debug pu - p0
    fit2 = LsqFit.curve_fit(rc_logistic_fun, x, y ,p0; lower = pl, upper = pu)
    at_upper = abs.(pu .- coef(fit2)) .<= 1E-10
    at_lower = abs.(pl .- coef(fit2)) .<= 1E-10
    if any(at_upper)
        @show p0
        @show pu
        @show @sprintf("rc_logistic_fit Hitting upper bound at index %s", findall(at_upper))
    end
    if any(at_lower)
        @show p0
        @show pl
        @show @sprintf("rc_logistic_fit Hitting lower bound at index %s", findall(at_lower))
    end
    coef(fit2)
end

Δ(x) = first(diff([extrema(x)...])) ## extrema returns tuple

function lin_i2i(x_interval, y_interval)
    a = Δ(y_interval) / Δ(x_interval)
    b = minimum(y_interval) - a*minimum(x_interval)
    x -> a.*x .+ b
end

function scale_fwd(x, range_in, range_out)
    ## map range_in to range_out
    @assert((length(range_in) == 2) & (length(range_out) == 2))
    lin_i2i(range_in, range_out)(x)
end


function scale_fwd(x;x_range)
    ## maps x to itself if x_range i [0,1]
    @assert length(x_range) == 2
    dx = first(diff(x_range))
    (x-x_range[1])/dx
end

function scale_rev(X;x_range)
    @assert length(x_range) == 2
    dx = first(diff(x_range))
    X*dx + x_range[1]
end

## Simulat ea Hill-type function with slight initial concavity
function sim_hill(; points=100, xmin = 0, xmax = 100, ymin = 0, ymax = 4, sd = 0.05, well = "A01", seed=missing)
    if !ismissing(seed)
        Random.seed!(seed)
    end
    dx = xmax - xmin
    dy = ymax - ymin
    xstep = dx/points
    t = collect(xmin:xstep:xmax);
    y_1 = PlateReaderCore.rc_exp.(t .- xmin, sqrt(dy),dx/5,ymin) ;
    y_2 = PlateReaderCore.rc_exp.(t .- xmin, sqrt(dy),dx/10,ymin) ;
    y = y_1.*y_2 .+ rand.(Normal.(0, sd));
    well = ReaderCurve(readerplate_well = well,
                       kinetic_time = t,
                       reader_value = y,
                       time_unit = "t",
                       value_unit = "y",
                       );
    well
end

function area_under_curve(x, y)
    Good = .!(ismissing.(x) .| ismissing.(y))
    if sum(Good) <= 1 ## less than 1 non-missing value: no area
        return(0)
    end
    X = x[Good]
    Y = y[Good]
    dX = diff(X)
    Y = Y .- Y[1] ## remove base-line
    Ymid = (Y[1:end-1] .+ Y[2:end])./2 ## use mid-points for area
    sum(dX .* Ymid) ## Area
end

area_under_curve(rc::ReaderCurve) = area_under_curve(rc.kinetic_time, rc.reader_value)
area_under_curve(rcf::ReaderCurveFit) = area_under_curve(rcf.readercurve.kinetic_time, rcf.predict.(rcf.readercurve.reader_value))
area_under_curve_ratio(rcf::ReaderCurveFit) =  area_under_curve(rcf) / area_under_curve(rcf.readercurve)
