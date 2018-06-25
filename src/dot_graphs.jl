
# import DTM: number_of_nodes,nodesize,maxdepth,create_leaves_array,fittedratio
export graph 

function graph(tree::Tree) #tree::Node,number_of_num_features::Int,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(undef,0))
    #https://en.wikipedia.org/wiki/DOT_(graph_description_language)
    #http://www.webgraphviz.com/
    #creates a digraph string in DOT language

    dot_graph=""
    node=tree.rootnode
    total_weight_of_tree = nodesize(node)
    candMatWOMaxValues=tree.candMatWOMaxValues
    mappings=tree.charMappings
    sett=tree.settings
    df_name_vector=tree.settings.df_name_vector
    if isa(node,Leaf)
        #the tree is only a leaf....
        dot_graph="! the tree is only a single leaf. No DOT graph was created."
    else
    #we have an actual tree
        orig_id=node.featid
        orig_id<0 ? this_id=tree.settings.number_of_num_features-orig_id : this_id=orig_id
        dot_graph  = "// Tree summary:\n"
        dot_graph  *= "// Total weight/exposure of tree: $(nodesize(tree))\n"
        dot_graph *= "// Number of leaves: $(size(create_leaves_array(tree),1))\n"
        dot_graph *= "// Number of Nodes: $(number_of_nodes(tree))\n"
        dot_graph *= "// Maximal depth of the tree: $(maxdepth(tree))\n\n\n"

        edge_to_the_left = if orig_id>0
             "$(sett.df_name_vector[this_id]) <= $(round(candMatWOMaxValues[this_id][node.subset[end]],sigdigits=5))"
        else
            if size(node.subset,1)==1
                "$(sett.df_name_vector[this_id]) = $(mappings[-orig_id][node.subset][1])"
            else
                "$(sett.df_name_vector[this_id]) in [$(join(mappings[-orig_id][[convert(Int,z) for z in node.subset]],","))]"
            end
        end

        label = "Fit $(round(fittedratio(node),sigdigits=3)) Weight $(round(nodesize(node),digits=1))"		

        name_of_this_node="R"
        dot_graph *= """$(name_of_this_node)[label="$(label)"]\n"""

        (n_nodes_and_leafsr,n_nodes_and_leafsl)=total_n_subnodes(node)

        tmpl=graph(node.left,1,name_of_this_node,sett,edge_to_the_left,mappings,candMatWOMaxValues,total_weight_of_tree)
        dot_graph *= "\n"
        dot_graph *= tmpl
        tmpr=graph(node.right,n_nodes_and_leafsl+1,name_of_this_node,sett,"",mappings,candMatWOMaxValues,total_weight_of_tree)
        dot_graph *= "\n"
        dot_graph *= tmpr
end
leading_comment="""

/*
This is a representation of the tree in DOT format. See also
http://www.webgraphviz.com/
https://en.wikipedia.org/wiki/DOT_(graph_description_language)
http://www.graphviz.org/content/dot-language
*/

"""
    return string("digraph tree {",leading_comment,dot_graph,"\n}")
    
end

function total_n_subnodes(node::Node{T}) where T<:Unsigned
    n_nodes_childl=number_of_nodes(node.left)
    n_leafs_childl=number_of_nodes(node.left)
    n_nodes_childr=number_of_nodes(node.right)
    n_leafs_childr=number_of_nodes(node.right)
    n_nodes_and_leafsr=n_nodes_childr+n_leafs_childr+1
    n_nodes_and_leafsl=n_nodes_childl+n_leafs_childl+1
    return n_nodes_and_leafsr,n_nodes_and_leafsl
end

function graph(node::Node{T},counter::Int,parentname::String,sett::ModelSettings,edge_from_parent_description::String,mappings,candMatWOMaxValues,total_weight_of_tree::Float64) where T<:Unsigned
    orig_id=node.featid
	orig_id<0 ? this_id=sett.number_of_num_features-orig_id : this_id=orig_id

    edge_to_the_left = if orig_id>0
         "$(sett.df_name_vector[this_id]) <= $(round(candMatWOMaxValues[this_id][node.subset[end]],sigdigits=5))"
    else
        if size(node.subset,1)==1
            "$(sett.df_name_vector[this_id]) = $(mappings[-orig_id][node.subset][1])"
        else
            "$(sett.df_name_vector[this_id]) in [$(join(mappings[-orig_id][[convert(Int,z) for z in node.subset]],","))]"
        end
    end
    
	label = "Fit $(round(fittedratio(node),sigdigits=3)) Weight $(round(nodesize(node),digits=1))"		
    

    name_of_this_node="N$(counter)"
    dot_graph = """$(name_of_this_node)[label="$(label)"]\n"""
    dot_graph *= """$(parentname)->$(name_of_this_node) [label="$(edge_from_parent_description)"]\n"""

    (n_nodes_and_leafsr,n_nodes_and_leafsl)=total_n_subnodes(node)

    tmpl=graph(node.left,counter+1,name_of_this_node,sett,edge_to_the_left,mappings,candMatWOMaxValues,total_weight_of_tree)
    dot_graph *= "\n"
    dot_graph *= tmpl
    tmpr=graph(node.right,counter+n_nodes_and_leafsl+1,name_of_this_node,sett,"",mappings,candMatWOMaxValues,total_weight_of_tree)
    dot_graph *= "\n"
    dot_graph *= tmpr

    return dot_graph
end

function graph(node::Leaf,counter::Int,parentname::String,sett::ModelSettings,edge_from_parent_description::String,mappings,candMatWOMaxValues,total_weight_of_tree::Float64)
    thisw = nodesize(node)
	label = "Leaf $(node.id)\nFit $(round(node.fitted,sigdigits=3))\n Weight $(round(thisw,digits=1)) ($(round(100*thisw/total_weight_of_tree,digits=1))%)"
    name_of_this_node="N$(counter)"
    dot_graph = """$(name_of_this_node)[label="$(label)"]\n"""
    dot_graph *= """$(parentname)->$(name_of_this_node) [label="$(edge_from_parent_description)"]\n"""
end

#
#    miha=graph(resulting_model);
# print(miha)



function draw_dot_graph(gvloc::String,dot_txt_file::String,outputfilename::String;format="pdf")
#
#This function prints a representaion of the tree based on an input specified in DOT format, see
#http://www.webgraphviz.com/
#https://en.wikipedia.org/wiki/DOT_(graph_description_language)
#http://www.graphviz.org/content/dot-language
#

#http://edutechwiki.unige.ch/en/Graphviz
# this should work
# "C:\Program Files (x86)\Graphviz2.38\bin\dot.exe mytree.dot.txt -Tpng -o c:\temp\mt.png"
# C:\Program Files (x86)\Graphviz2.38\bin\dot.exe mytree.dot.txt -Tpdf -o c:\temp\mt.pdf
 #dot mytree.dot.txt -Tsvg -o c:\temp\mt.svg
# gvloc="C:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe"

#=
gvloc="C:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe"
outfilename = "c:\\temp\\mo ha\\mi.pdf"
isdir(splitdir(outfilename)[1])
ma=draw_dot_graph(gvloc,"c:\\temp\\maroxi.dot.txt",outfilename)
=#

if length(gvloc)<1 #|| length(graph_as_string)<1
    return nothing
else
    if !isfile(gvloc)
        @warn("DTM: GraphViz executable not found at: \n $(gvloc)")
    else
        if !isfile(dot_txt_file)
            @warn("DTM: GraphViz dot file not found at: \n $(dot_txt_file)")
        else
            if format=="pdf"
                command="-Tpdf"
                # "dot mytree.dot.txt -Tpdf -o c:\temp\mt.pdf"
            else
                @warn("DTM: invalid format for GraphViz was specified. Format = $(format)")
            end
        end
    end
    middle_command="$(dot_txt_file)" #" $(command) -o"
    #command *= " -o $(outputfilename)"
    actual_command=`$(gvloc) $(middle_command) $(command) -o $(outputfilename)`
    chosen_dir = splitdir(outputfilename)[1]
    if isdir(chosen_dir)
        try
            #run(`$(gvloc) $(dot_txt_file) $(command)`)
            #run(`$(command2[])`)
            #GOAL: run(`'C:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe' 'c:\\temp\\maroxi.dot.txt' -Tpdf -o 'c:\\temp\\mi.pdf'`)
            run(actual_command)
			println("Graph was visualized in the file: \n $(outputfilename)")				
        catch e
            @show actual_command
            @warn("DTM: GraphViz failed to crate a graph.")
            @show e
        end
    else
        @warn("DTM: Inexistent location specified for graph file: $(chosen_dir)")
    end

return nothing
end

return nothing
end
