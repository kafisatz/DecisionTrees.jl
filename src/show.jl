
function Base.show(io::IO, t::MIME"text/plain",x::Union{DecisionTreesv06.Tree,DecisionTreesv06.Leaf,DecisionTreesv06.Node,DecisionTreesv06.BaggedTree,DecisionTreesv06.BoostedTree})
    elt=string(typeof(x))
    println(io,elt)
    println(io,"")
    show(io,t,fieldnames(x))
    println(io,"")
    println(io,"")
    println(io,"Use the dump() function to see the whole content of the object.")
end

function Base.show(io::IO,x::Union{DecisionTreesv06.Tree,DecisionTreesv06.Leaf,DecisionTreesv06.Node,DecisionTreesv06.BaggedTree,DecisionTreesv06.BoostedTree})
    elt=string(typeof(x))
    println(io,"")
    println(io,elt)
    println(io,"")
    show(io,fieldnames(x))
    println(io,"")
    println(io,"")
    println(io,"Use the dump() function to see the whole content of the object.")
end
