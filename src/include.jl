let
    global Include
    existing_files = Set{AbstractString}([])
    function Include(path::AbstractString; filename::AbstractString = "None file", linenum::Integer = -1)
        if (path ∉ existing_files)
            print("include the julia file of: "); println(path)
            include(string(@__DIR__) * "/" * path)
            push!(existing_files, path)
        else
            println("TODO: avoid re-inclusion of the julia file of: " * string(path) * " at " * string(filename) * ":" * string(linenum));
        end
    end
end

macro Include(path)
    :(Include($path, filename = $(string(__source__.file)), linenum = $(__source__.line)))
end
