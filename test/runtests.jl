using PlateReaderCore
using Test
using DataFrames, CSV, DataFramesMeta

@testset "PlateReaderCore.jl" begin

end

@testset "ReaderCurves.jl" begin
    s1 = collect(0:10:100)
    y1 = PlateReaderCore.rc_exp(s1, 4, 100, 0.05) ## 0.05 .+ 4 .* (1 .- exp.(-s1 ./100)),
    A01 = ReaderCurve(readerplate_well = "A01",
                      kinetic_time = s1,
                      reader_value = y1,
                      time_unit = "sec",
                      value_unit = "OD405nm",
                      )
    @test A01.reader_value[11] == 0.05 + 4 * (1 - exp(-1))
    @testset "Fit" begin
        A01_fit1 = rc_fit(A01,"linreg_trim")
        A01_fit2 = rc_fit(A01,"max_slope")
        A01_fit3 = rc_fit(A01,"smooth_spline")
        A01_fit4 = rc_fit(A01,"L4P") # ; lambda = 250
        @test A01_fit1.predict.([0,1]) == [A01_fit1.intercept, A01_fit1.intercept + A01_fit1.slope]
        @test A01_fit2.predict.([0,1]) == [A01_fit2.intercept, A01_fit2.intercept + A01_fit2.slope]
        @test isapprox(A01_fit3.predict.([0,1]), [0.05000000004687248, 0.08883842876426919] ) ## [0.050046522574332335, 0.08886823287981804]) ## [0.05655762222793344, 0.09394291903679455])
        @test isapprox(A01_fit4.predict.([0,1]), [0.05138458324982054, 0.09088904053992053])
    end
    @testset "fit with missing" begin
        s1 = collect(0:10:100)
        A02 = ReaderCurve(readerplate_well = "A02",
                          kinetic_time = s1,
                          reader_value = replace( 2 .* s1, 60 => NaN, 20 => Inf, 40=> -Inf),
                          time_unit = "sec",
                          value_unit = "OD405nm",
                          )
        A03 = ReaderCurve(readerplate_well = "A03",
                          kinetic_time = s1,
                          reader_value = repeat([NaN], length(s1)),
                          time_unit = "sec",
                          value_unit = "OD405nm",
                          )
        A04 = ReaderCurve(readerplate_well = "A04",
                          kinetic_time = s1,
                          reader_value = repeat([Inf], length(s1)),
                          time_unit = "sec",
                          value_unit = "OD405nm",
                          )
        A05 = ReaderCurve(readerplate_well = "A05",
                          kinetic_time = s1,
                          reader_value = repeat([-Inf], length(s1)),
                          time_unit = "sec",
                          value_unit = "OD405nm",
                          )
        A02_fit1 = rc_fit(A02,"linreg_trim")
        A03_fit1 = rc_fit(A03,"linreg_trim")
        A04_fit1 = rc_fit(A04,"linreg_trim")
        A05_fit1 = rc_fit(A05,"linreg_trim")

        A02_fit2 = rc_fit(A02,"max_slope")
        A03_fit2 = rc_fit(A03,"max_slope")
        A04_fit2 = rc_fit(A04,"max_slope")
        A05_fit2 = rc_fit(A05,"max_slope")

        A02_fit3 = rc_fit(A02,"smooth_spline")
        A03_fit3 = rc_fit(A03,"smooth_spline")
        A04_fit3 = rc_fit(A04,"smooth_spline")
        A05_fit3 = rc_fit(A05,"smooth_spline")

        A02_fit4 = rc_fit(A02,"L4P") ## Terrible fit!
        A03_fit4 = rc_fit(A03,"L4P")
        A04_fit4 = rc_fit(A04,"L4P")
        A05_fit4 = rc_fit(A05,"L4P")

        Plate_2 = ReaderPlate(
        readerplate_id = "1117389a-9abe-4cf9-9feb-ec7fa1aa0945",
        readerplate_barcode = "",
        readerfile_name = "testFile",
        readerplate_geometry = 96,
        readercurves = [A02,A03,A04,A05]
        )
        
        Plate_3 = ReaderPlateFit(
        readerplate_id = "1117389a-9abe-4cf9-9feb-ec7fa1aa0931",
        readerplate_barcode = "",
        readerfile_name = "testFile",
        readerplate_geometry = 96,
        readercurves = [A02_fit3,A03_fit3,A04_fit3,A05_fit3]
        )
        
    end    
end
@testset "Bubble" begin
    B01_df = CSV.File("b01.csv") |> DataFrame
    B01 = ReaderCurve(readerplate_well = "B01",
                      kinetic_time = B01_df.kinetic_sec,
                      reader_value = B01_df.absorbance_value,
                      time_unit = "sec",
                      value_unit = "OD405"
                      )
    B01_fit = rc_fit(B01,"linreg_trim")
    @test B01_fit.predict.([0,1]) == [B01_fit.intercept, B01_fit.intercept + B01_fit.slope]
    B01_fit2 = rc_fit(B01,"max_slope")
    @test B01_fit2.predict.([0,1]) == [B01_fit2.intercept, B01_fit2.intercept + B01_fit2.slope]
    B01_fit3 = rc_fit(B01,"smooth_spline") ## ; lambda = 250)
    @test isapprox(B01_fit3.predict.([0,1]), [0.08772661189719624, 0.087772176580951] ) ## [0.0876330427649209, 0.08768712387504587]) ## [0.08762157518257128, 0.08767650822148132]
    B01_fit4 = rc_fit(B01,"L4P" )
    @test isapprox(B01_fit4.predict.([0,1]), [0.0757660594695346, 0.0761047195492072])
end
@testset "Plate" begin
    s1 = collect(0:10:100)
    y1 = PlateReaderCore.rc_exp(s1, 4, 100, 0.05) ## 0.05 .+ 4 .* (1 .- exp.(-s1 ./100)),
    A01 = ReaderCurve(readerplate_well = "A01",
                      kinetic_time = s1,
                      reader_value = y1,
                      time_unit = "sec",
                      value_unit = "OD405nm",
                      )
    A02 = ReaderCurve(readerplate_well = "A02",
                      kinetic_time = s1,
                      reader_value = y1 .+ 0.1,
                      time_unit = "sec",
                      value_unit = "OD405nm",
                      )
    A03 = ReaderCurve(readerplate_well = "A03",
                      kinetic_time = s1,
                      reader_value = y1 .- 0.1,
                      time_unit = "sec",
                      value_unit = "OD405nm",
                      )
    Plate_1 = ReaderPlate(
        readerplate_id = "1117389a-9abe-4cf9-9feb-ec7fa1aa0933",
        readerplate_barcode = "Barcode",
        readerfile_name = "testFile",
        readerplate_geometry = 96,
        readercurves = [A01, A02, A03]
    )
    Plate_f1 = ReaderPlateFit(
        readerplate_id = "1117389a-9abe-4cf9-9feb-ec7fa1aa0932",
        readerplate_barcode = "Barcode",
        readerfile_name = "testFile",
        readerplate_geometry = 96,
        readercurves = [rc_fit(A01, "smooth_spline"), rc_fit(A02, "smooth_spline"), rc_fit(A03, "smooth_spline")]
    )
    
    @test length(Plate_1) == 3
end
@testset "DataFrames" begin
    dat1_df = xlsx("dat1.xlsx"; sheet=1)
    dat1 = ReaderRun(
	@transform(
		rename(dat1_df, :well_name => :readerplate_well, 
			:geometry => :readerplate_geometry, 
			:kinetic_sec => :kinetic_time, 
			:absorbance_value => :reader_value, 
			:reader_temperature_C => :reader_temperature), 
	equipment = "Neo2", 
	software="Gen5", 
	run_starttime = missing,
	time_unit = "sec", 
	value_unit = "OD405nm", 
	temperature_unit="C",
	readerplate_id = :readerfile_name, 
        ));
    @test length(dat1) == 1
    @test length(dat1.readerplates[1]) == geometry(dat1)
    dat1_fit1 = rc_fit(dat1, "linreg_trim")
    ## Subset on a quadrant
    dat1_fit1_Q3 = Q(dat1_fit1.readerplates[1], "Q3")
    @test  well_names(dat1_fit1_Q3) == ["B01", "B03", "B05", "B07", "B09", "B11", "B13", "B15", "B17", "B19", "B21", "B23", "D01", "D03", "D05", "D07", "D09", "D11", "D13", "D15", "D17", "D19", "D21", "D23", "F01", "F03", "F05", "F07", "F09", "F11", "F13", "F15", "F17", "F19", "F21", "F23", "H01", "H03", "H05", "H07", "H09", "H11", "H13", "H15", "H17", "H19", "H21", "H23", "J01", "J03", "J05", "J07", "J09", "J11", "J13", "J15", "J17", "J19", "J21", "J23", "L01", "L03", "L05", "L07", "L09", "L11", "L13", "L15", "L17", "L19", "L21", "L23", "N01", "N03", "N05", "N07", "N09", "N11", "N13", "N15", "N17", "N19", "N21", "N23", "P01", "P03", "P05", "P07", "P09", "P11", "P13", "P15", "P17", "P19", "P21", "P23"]
    @test dat1_fit1_Q3.readerplate_geometry == 384
    dat1_fit1_Q3_96 = Q(dat1_fit1.readerplates[1], "Q3"; well96 = true)
    @test well_names(dat1_fit1_Q3_96) == ["A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12", "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B09", "B10", "B11", "B12", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12", "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09", "E10", "E11", "E12", "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12", "G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", "G11", "G12", "H01", "H02", "H03", "H04", "H05", "H06", "H07", "H08", "H09", "H10", "H11", "H12"]
    @test dat1_fit1_Q3_96.readerplate_geometry == 96
    dat1_fit1_C13 = well(dat1_fit1.readerplates[1], "C13")
    @test dat1_fit1_C13.readercurve.readerplate_well == "C13"
    dat1_1_Q2 = Q(dat1.readerplates[1], "Q2")
    @test well_names(dat1_1_Q2) == ["A02", "A04", "A06", "A08", "A10", "A12", "A14", "A16", "A18", "A20", "A22", "A24", "C02", "C04", "C06", "C08", "C10", "C12", "C14", "C16", "C18", "C20", "C22", "C24", "E02", "E04", "E06", "E08", "E10", "E12", "E14", "E16", "E18", "E20", "E22", "E24", "G02", "G04", "G06", "G08", "G10", "G12", "G14", "G16", "G18", "G20", "G22", "G24", "I02", "I04", "I06", "I08", "I10", "I12", "I14", "I16", "I18", "I20", "I22", "I24", "K02", "K04", "K06", "K08", "K10", "K12", "K14", "K16", "K18", "K20", "K22", "K24", "M02", "M04", "M06", "M08", "M10", "M12", "M14", "M16", "M18", "M20", "M22", "M24", "O02", "O04", "O06", "O08", "O10", "O12", "O14", "O16", "O18", "O20", "O22", "O24"]
    dat1_1_Q2_96 = Q(dat1.readerplates[1], "Q2", well96 = true)
    @test well_names(dat1_1_Q2_96) == ["A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12", "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B09", "B10", "B11", "B12", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12", "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09", "E10", "E11", "E12", "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12", "G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", "G11", "G12", "H01", "H02", "H03", "H04", "H05", "H06", "H07", "H08", "H09", "H10", "H11", "H12"]
    dat1_C13 = well(dat1.readerplates[1], "C13")
    @test dat1_C13.readerplate_well == "C13"    
end
@testset "normalize nonlinear" begin
    t1 = collect(0:100:1000)
    y1 = PlateReaderCore.rc_exp.(t1, 4,500,0.05)
    # @test isapprox(PlateReaderCore.scale_rev.(PlateReaderCore.scale_fwd.(t1;x_range= [100,1000]); x_range= [100,1000]), t1)
    # @test isapprox(PlateReaderCore.scale_rev.(PlateReaderCore.scale_fwd.(y1;x_range= [-2,4]); x_range= [-2,4]), y1)
    # @test isapprox(PlateReaderCore.scale_fwd.(t1;x_range= [0,1]) ,t1) ## Not what I want!
    @test isapprox(PlateReaderCore.lin_i2i([0,1], [1,11])(collect(0:.1:1)), collect(1:11))
    @test isapprox(PlateReaderCore.scale_fwd(t1, extrema(t1), extrema(t1)), t1) ## OBS signature!
    @test isapprox(PlateReaderCore.scale_fwd(PlateReaderCore.scale_fwd(t1, extrema(t1), [0,1]), [0,1], extrema(t1)), t1) ## OBS signature!
end
@testset "DataFrame from file" begin
    dat2_df = xlsx("dat_ex.xlsx"; sheet=1)
    dat2 = ReaderRun(dat2_df)
    @test length(dat2) == 1
end
@testset "Area undercurve" begin
    @test PlateReaderCore.area_under_curve([1,2,3,4,5,6,7,8,9], [1,2,3,4,5,6,7,8,9]) == 32
    @test PlateReaderCore.area_under_curve([-1,9], [-1,9]) == PlateReaderCore.area_under_curve([0,10], [0,10])
    A01 = PlateReaderCore.sim_hill(;points = 100, xmax=1, sd=.1, seed=123, well= "A01")
    A01_fit = rc_fit(A01, "smooth_spline";lambda = 1E-3, x_range = [0,1], y_range = [0,1])
    @test isapprox(PlateReaderCore.area_under_curve(A01), 2.956311533459582)
    @test isapprox(PlateReaderCore.area_under_curve(A01_fit), 2.0307829099202612)
    @test isapprox(PlateReaderCore.area_under_curve_ratio(A01_fit),0.6869312949382456)
end
