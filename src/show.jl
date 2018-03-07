
function Base.show(io::IO, t::MIME"text/plain",x::Union{Tree,Leaf,Node,BaggedTree,BoostedTree})
    elt=string(typeof(x))
    println(io,elt)
    println(io,"")    
    show(io,t,fieldnames(x))
    show(io,t,"Use the dump() function to see the whole content of the object.")
end

