@Include "InputParser.jl"
@Include "Diagnostics.jl"

module Result
    using Plots
    using FileIO
    using ..InputParser
    using ..Diagnostics
    using DataStructures

    pyplot()

    function print_map_z_y(smax, output_dir::String, col_annotation::String)
        N, M = size(smax.Z)
        S = zeros(N, M)
        T = zeros(N, M)
        for i in 1:N
            for j in 1:M
                S[i, j] = smax.S[j]
                T[i, j] = smax.T[j]
            end
        end
        println("===================")
        print("N: ")
        println(N)
        print("M: ")
        println(M)
        println("===================")

        _view_matrix_txt(smax.Z .- 1, zero_based = true, file_path = output_dir * "/Z.txt", row_tag = "Muts", col_tag = "sampleID",
                         col_annotation = col_annotation, col_index = 1)
        _view_matrix_txt(smax.Y .- 1, zero_based = true, file_path = output_dir * "/Y.txt", row_tag = "Muts", col_tag = "sampleID",
                         col_annotation = col_annotation, col_index = 1)
        _view_matrix_txt(S .- 1, zero_based = true, file_path = output_dir * "/S.txt", row_tag = "Muts", col_tag = "sampleID",
                         col_annotation = col_annotation, col_index = 1)
        _view_matrix_txt(T .- 1, zero_based = true, file_path = output_dir * "/T.txt", row_tag = "Muts", col_tag = "sampleID",
                         col_annotation = col_annotation, col_index = 1)
    end

    function parse_list_at(tsv_file::String, col_index::I)::Deque{String} where {I <: Integer}
        ans = Deque{String}()
        if !isfile(tsv_file)
            return ans
        end

        open(tsv_file, "r") do fp
            for line::String in eachline(fp)
                if startswith(line, "#")
                    continue
                end

                line = strip(line)
                line_col::Array{String,1} = split(line, '\t')
                if 1 <= col_index <= length(line_col)
                    push!(ans, line_col[col_index])
                else
                    ans = Deque{String}()
                    return ans
                end
            end
        end
        return ans
    end

    function parse_comment_at(tsv_file::String)
        ans = []
        if !isfile(tsv_file)
            return ans
        end

        open(tsv_file, "r") do fp
            for line::String in eachline(fp)
                if !startswith(line, "#")
                    return ans
                end
                push!(ans, strip(line) * '\n')
            end
        end
        return ans
    end

    function _view_matrix_txt(matrix;
                              zero_based = true, file_path = "",
                              row_tag = "Muts", col_tag = "sampleID",
                              col_annotation = "", col_index = 1)
        R, C = size(matrix)
        annotation_deque::Deque{String} = parse_list_at(col_annotation, col_index)
        (C != length(annotation_deque)) && (annotation_deque = Deque{String}())

        open(file_path, "w") do fp
            output_list = parse_comment_at(col_annotation)

            for output_string in output_list
                write(fp, output_string)
            end
            #     push!(output_list, "#" * row_tag * "/" * col_tag)
            #     for r in 1:R
            #         (zero_based) && (push!(output_list, r-1))
            #         (!zero_based) && (push!(output_list, r))
            #     end
            #     output_string = join( map( x -> string(x), output_list ), '\t') * '\n'
            #     write(fp, output_string)

            for c in 1:C
                output_list = []
                if isempty(annotation_deque)
                    (zero_based)  && (push!(output_list, c-1);)
                    (!zero_based) && (push!(output_list, c);)
                else
                    push!(output_list, popfirst!(annotation_deque))
                end
                for r in 1:R
                    push!(output_list, matrix[r, c])
                end
                output_string = join(map(x->string(x), output_list), '\t') * '\n'
                write(fp, output_string)
            end
        end
    end
end
