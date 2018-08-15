
#=
pt=normpath(joinpath(pwd(),"testing","iterator.jl"))
@assert isfile(pt) 
include(pt)
=#

import Base: length,eltype,==,hash,iterate,append!,isless,resize!,convert

struct MyGrayCodeSubsetsHALF #<: DTSubsets
		size_of_set::Int
		stopat::Int #(1 << it.size_of_set)-2
	end

	eltype(it::MyGrayCodeSubsetsHALF) = Int #Array{eltype(1),1}
	length(it::MyGrayCodeSubsetsHALF) = it.stopat+1 #1 << it.size_of_set
	start(it::MyGrayCodeSubsetsHALF) = zero(1) #state = simply count through the subsets
	done(it::MyGrayCodeSubsetsHALF,state) = state>it.stopat  #if we have equality here, then the full set will be shifted / this is not needed for a binary split, hence we stop beforehand

	function next(it::MyGrayCodeSubsetsHALF, state)
	#returns the element which needs to be flipped
	  bit_to_flip=1+trailing_ones(state) #in the function print_graycode bit_to_flip is defined without the plus 1
	  state+=one(1)
	  bit_to_flip,state
	end
    @inline iterate(x::MyGrayCodeSubsetsHALF) = next(x, start(x))
    @inline iterate(x::MyGrayCodeSubsetsHALF, i) = done(x, i) ? nothing : next(x, i)

	bitflip_graycode_subsetsHALF(xs) = MyGrayCodeSubsetsHALF(length(xs),(1 << (length(xs)-1))-2)
    
    subs=bitflip_graycode_subsetsHALF(["ab","cd","mixa"])
    
    for i in subs
        @show i
    end