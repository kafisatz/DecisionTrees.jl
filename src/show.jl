
function Base.show(io::IO, t::MIME"text/plain",x::Union{DecisionTrees.Tree,DecisionTrees.Leaf,DecisionTrees.Node,DecisionTrees.BaggedTree,DecisionTrees.BoostedTree})
    elt=string(typeof(x))
    println(io,elt)
    #println(io,"")
    show(io,t,fieldnames(typeof(x)))
    #println(io,"")
    #println(io,"")
    #println(io,"Use the dump() function to see the whole content of the object.")
end

function Base.show(io::IO,x::Union{DecisionTrees.Tree,DecisionTrees.Leaf,DecisionTrees.Node,DecisionTrees.BaggedTree,DecisionTrees.BoostedTree})
    elt=string(typeof(x))
    #println(io,"")
    println(io,elt)
    #println(io,"")
    show(io,fieldnames(typeof(x)))
    #println(io,"")
    #println(io,"")
    #println(io,"Use the dump() function to see the whole content of the object.")
end
