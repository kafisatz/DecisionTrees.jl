#


function createGridSearchSettings(sett::ModelSettings;args...)
    settingsVector=[sett]
	seen=[]
	for (smember,valuelist) in args
		if in(smember,seen)
			error("You provided multiple values for the field $(smember)")
		end
        push!(seen,smember)
        try 
            sz=size(valuelist)
            @assert length(sz)==1 "DTM: expected a vector for $(smember). You provided $(typeof(valuelist)). \r\n$(valuelist)"
            n=sz[1]
            @assert n>0
            
            settingsVectorTemplate=deepcopy(settingsVector)
            
            for i=1:n
                #set value for this 'batch'  
                v=valuelist[i]        
                thisSettingsBatch=deepcopy(settingsVectorTemplate)      
                for kk=1:length(thisSettingsBatch)
                    s=thisSettingsBatch[kk]
                    #updateSettingsMod!(s,smember=)
                    oldValue=getfield(s,smember)
                    if typeof(v)!=typeof(oldValue)
                        valConverted=convertFromString(oldValue,string(v)) #tbd maybe string(v) could/should be replaced with v here
                        setfield!(s,smember,valConverted)
                    else
                        setfield!(s,smember,v)
                    end		
                    res=DecisionTrees.checkIfSettingsAreValid(s)
                    @assert res==nothing "DTM: Some settings are invalid! Abort."
                end
                if i==1
                    settingsVector=deepcopy(thisSettingsBatch)
                else 
                    append!(settingsVector,thisSettingsBatch)
                end    
            end
        catch e            
            @show smember,size(valuelist)
            @show valuelist
            @show e
            throw("DTM: unable to create settings vector.")
        end
    end
    return settingsVector
end


#=
valuelist=rand(9)
sz=size(valuelist)
length(sz)
sz[1]
size(22)[1]
=#