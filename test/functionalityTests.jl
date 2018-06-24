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
    hdr=resm2.exceldata.sheets[3].data[1,:] 
    hdr2=map(x->hdr[x][1],1:length(hdr))
    liftvalc=indexin(["Lift Val"],hdr2)[1]
    lifttrnc=indexin(["Lift Trn"],hdr2)[1]
    @test !(liftvalc==nothing)
    @test !(lifttrnc==nothing)
    if liftvalc!=nothing
        @test liftvalc>0
        liv=resm2.exceldata.sheets[3].data[sett.niter+1,liftvalc]
    end
    if lifttrnc!=nothing
        @test lifttrnc>0
        lit=resm2.exceldata.sheets[3].data[sett.niter+1,lifttrnc]
    end
    @test abs(liv-6.332896484066747) <tolForTheseTests
    @test abs(lit-6.610306331060652) <tolForTheseTests
end

end #testset
