using LightGraphs
using Random

mutable struct Nodes
    parent::Array{Int}
    size::Array{Int}
    num_components::Int
    
    function Nodes(size)
        new(Array(collect(1:size)), zeros(Int, size), size)
    end
end

function findroot(cluster::Nodes, i)
    root = i
    while root != cluster.parent[root]
        root = cluster.parent[root]
    end
    
    #
    # Path Compression
    #
    while i != root
        next = cluster.parent[i]
        cluster.parent[i] = root
        i = next
    end
    return root
end

function unionfind!(cluster::Nodes, i, j)
    rooti = findroot(cluster, i)
    rootj = findroot(cluster, j)
    if rootj != rooti
        if cluster.size[rooti] > cluster.size[rootj]
            cluster.parent[rootj] = rooti
            cluster.size[rooti] += cluster.size[rootj]
        else
            cluster.parent[rooti] = rootj
            cluster.size[rootj] += cluster.size[rooti]
        end
        cluster.num_components -= 1
    end
end

function is_connected(cluster::Nodes, i, j)
    return findroot(cluster, i) == findroot(cluster, j)
end

function components(cluster::Nodes)
    return cluster.num_components
end

function site_percolation(g::AbstractGraph)
    N = nv(g)
    result = ones(Int, N)

    # collect and shuffle the nodes id 
    vertice_set = collect(vertices(g));
    shuffle!(vertice_set);

    nodes = Nodes(N)
    
    for (i,v) in enumerate(vertice_set)
        nodes.size[v] = 1 # activate the site v
        for neighbor in outneighbors(g,v)
            if nodes.size[neighbor] > 0 # check if the site has already been activated
                unionfind!(nodes, v, neighbor)
            end
        end
        result[i] += maximum(nodes.size)
    end
    return result
end

function bond_percolation(g::AbstractGraph)
    N = nv(g)
    M = ne(g)
    result = ones(Int, M)
    
    # collect and shuffle the edges 
    edges_set = collect(edges(g));
    shuffle!(edges_set);

    nodes = Nodes(N)
    nodes.size = ones(N)
    
    for (i,e) in enumerate(edges_set)
        unionfind!(nodes, src(e), dst(e))
        result[i] += maximum(nodes.size)
    end
    return result
end