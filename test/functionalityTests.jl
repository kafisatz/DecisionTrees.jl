#smoketest_tree.jl

@testset "Functionality Tests" begin

import Random: rand,randstring

tolForTheseTests=1e-15
##################################################
#Read the data
##################################################

datafile=joinpath(datadir,"freMTPL2","freMTPL2.csv")
@test isfile(datafile)
if !isfile(datafile)
    @info("Data not found. Trying to unzip the data: $(datafile)")
    try 
        pathToZip=string(splitext(datafile)[1],".zip") 
        unzipLoc=splitdir(pathToZip)[1]        
        if Sys.iswindows()            
            cmd=`powershell.exe -nologo -noprofile -command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('$(pathToZip)', '$(unzipLoc)'); }"`
            cmd2=`7z e $(pathToZip) -o$(unzipLoc)`
            try
                run(cmd2)
            catch
            end
        else 
            cmd=`unzip $(pathToZip) -d $(unzipLoc)`
        end
        if !isfile(datafile)
            run(cmd)
        end
    catch erO
        @warn("Failed to unzip data. Some tests will not run.")
        @show erO        
    end
end
if !isfile(datafile)
    @test "CSV DATA NOT FOUND"=="freMTPL2.csv was not found. Functionality tests not run."
    @warn("You may want to manually unzip the zip file in this location $(splitdir(datafile)[1]) as some of the tests rely on that CSV file.")
else
    @time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

    #add AreaInteger as new variable, i.e. A:F -> 1:6
        areasSorted=sort(unique(fullData[:Area]))
        AreaInteger=map(x->findall((in)([x]),areasSorted)[1],fullData[:Area])
        fullData[:AreaInteger]=AreaInteger

    #correct for unreasonable observations
    for i=1:size(fullData,1)
        fullData[:ClaimNb][i]=min(4,fullData[:ClaimNb][i])
        fullData[:Exposure][i]=min(1.0,fullData[:Exposure][i])
    end
        
    #set independent variables
    selected_explanatory_vars=["Area","AreaInteger","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"]

    ##############################
    #Prepare the data
    ##############################

    dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);
    updateSettingsMod!(sett,minw=-0.03,model_type="build_tree",boolCalculatePoissonError=true)

    ##############################
    #Run single tree model
    ##############################

    resultingFiles,resM=dtm(dtmtable,sett)
    stats1=resM.exceldata.sheets[2].data[63:end,1:2]
    @test isapprox(stats1[1,2],0.9902398642839035,atol=tolForTheseTests)
    @test isapprox(stats1[5,2],7.893600119119564,atol=tolForTheseTests)
    @test isapprox(stats1[6,2],9.585992687936427,atol=tolForTheseTests)
    @test isapprox(stats1[19,2],0.31235075486038566,atol=tolForTheseTests)
    @test isapprox(stats1[20,2],0.322243286637449,atol=tolForTheseTests)
    @warn("todo: add tests on lift, poisson error, etc.")

    ##################################################
    #run boosting
    ##################################################

    sett.niter=15
    sett.mf=0.1
    sett.model_type="boosted_tree"
    strs,resm2=dtm(dtmtable,sett)
    @test typeof(resm2.modelstats)==DataFrame
    headr=resm2.exceldata.sheets[3].data[1,:] 
    headr2=map(x->headr[x][1],1:length(headr))
    hdrForLoop=["Correlation" "Std Dev Trn" "Std Dev Val" "Std Dev Ratio" "Lift Trn" "Lift Val" "Reversals Trn" "Reversals Val" "RelDiff TrnObs vs ValObs" "RelDiff ValObs vs TrnFitted" "Min Rel Trn" "Max Rel Trn" "Min Observed Ratio Trn" "Max Observed Ratio Trn" "Min Rel Val" "Max Rel Val" "Min Observed Ratio Val" "Max Observed Ratio Val" "Gini Num (Two Arg Fn) Unweighted Trn" "Gini Num (Two Arg Fn) Unweighted Val" "Gini Num Trn (Single Agrument Fn)" "Gini Num Val (Single Argument Fn)" "RSquared of Ratio Trn" "RSquared of Ratio Val" "RSquared of Numerator Trn" "RSquared of Numerator Val" "RSS of Ratio Trn" "RSS of Ratio Val" "RSS of Numerator Trn" "RSS of Numerator Val" "Average Poisson Error Trn" "Average Poisson Error Val"][:]
    #note: you need to add digits in Excel when copy pasting (otherwise pasting "as text" will not show all digits)
    expectedVals=[0.994468759282649000000 0.002375403930973750000 0.002352411531960320000 0.990320636118507000000 6.843274521630170000000 6.450831459366510000000 0.000000000000000000000 1.000000000000000000000 0.050588136489102300000 1.182553313076510000000 0.402697426194537000000 2.755769036603120000000 0.040359671328647100000 0.276192310504699000000 0.423331557393917000000 2.730840528179300000000 0.044077650881683200000 0.284337496962536000000 1.000000000000000000000 1.000000000000000000000 1.000000000000000000000 1.000000000000000000000 0.000000000000000000000 0.000000000000000000000 0.000000000000000000000 0.000000000000000000000 1.000000000000000000000 1.000000000000000000000 1.000000000000000000000 1.000000000000000000000 0.311428129785048000000 0.321224667652410000000][:]
    @test size(expectedVals)==size(hdrForLoop)
    for i=1:length(expectedVals)
        location=indexin([hdrForLoop[i]],headr2)[1]
        @test !(location==nothing)
        if location!=nothing
            @test location>0
            value=resm2.exceldata.sheets[3].data[sett.niter+1,location]
            @test isapprox(value,expectedVals[i])
        end
    end
    #liftvalc=indexin(["Lift Val"],headr2)[1]
    #lifttrnc=indexin(["Lift Trn"],headr2)[1]
    #@test !(liftvalc==nothing)
    #@test !(lifttrnc==nothing)
    #if liftvalc!=nothing
    #    @test liftvalc>0
    #    liv=resm2.exceldata.sheets[3].data[sett.niter+1,liftvalc]
    #end
    #if lifttrnc!=nothing
    #    @test lifttrnc>0
    #    lit=resm2.exceldata.sheets[3].data[sett.niter+1,lifttrnc]
    #end    
    #@show lit,liv
    #@test abs(liv-6.45083145936651) <tolForTheseTests
    #@test abs(lit-6.843274521630169) <tolForTheseTests   
end

end #testset
