
module InputParser

    function parse_input_summary(summary_path::String,
                                 transpose::Bool      = false,
                                 INT::Type{<:Integer} = Int32,
                                 REAL::Type{<:Real}   = Float32)::Array{REAL, 2}
        if transpose
            return _parse_input_summary_sample_mutation(summary_path, INT, REAL)
        else
            return _parse_input_summary_mutation_sample(summary_path, INT, REAL)
        end
    end

    function _parse_input_summary_sample_mutation(summary_path::String,
                                                  INT::Type{<:Integer} = Int32,
                                                  REAL::Type{<:Real}   = Float32)::Array{REAL, 2}
        ans = []
        open(summary_path, "r") do f
            for line::String in eachline(f)
                if startswith(line, "#")
                    continue
                end
                line = strip(line)
                lineCol::Array{String,1} = split(line, '\t')
                nodeID   = lineCol[1]
                sampleID = lineCol[2]
                push!(ans, map(x -> parse(REAL, x), lineCol[3:length(lineCol)]) )
            end
        end
        S = length(ans)
        M = 0
        if S > 0
            M = length(ans[1])
        end
        mat = zeros(REAL, S, M)
        for (s,m) in Iterators.product(1:S,1:M)
            mat[s,m] = ans[s][m]
        end
        return mat
    end

    function _parse_input_summary_mutation_sample(summary_path::String,
                                                 INT::Type{<:Integer} = Int32,
                                                 REAL::Type{<:Real}   = Float32)::Array{REAL, 2}
        ans = []
        open(summary_path, "r") do f
            for line::String in eachline(f)
                if startswith(line, "#")
                    continue
                end
                line = strip(line)
                lineCol::Array{String,1} = split(line, '\t')
                mutationID   = lineCol[1]
                push!(ans, map(x -> parse(REAL, x), lineCol[2:length(lineCol)]) )
            end
        end
        M = length(ans)
        S = 0
        if M > 0
            S = length(ans[1])
        end
        mat = zeros(REAL, S, M)
        for (s,m) in Iterators.product(1:S,1:M)
            mat[s,m] = ans[m][s]
        end
        return mat
    end

    # using ConfParser
    # export Parameters
    # export Annealer
    #
    # mutable struct Annealer{I <: Integer, R <: Real}
    #     num_temperature::I
    #     ln_p_ladders_main::Array{R, 2}
    #     ln_p_ladders_err::Array{R, 2}
    #     p_temperature::Array{R, 1}
    #     q_temperature::Array{R, 2}
    #     coolest_state::I
    # end
    #
    # mutable struct Parameters{I <: Integer, R <: Real}
    #     α_s::R
    #     α_v::R # dirichlet process hyper parameter for sample/variant
    #
    #     δ_s::Array{R, 1} # freq vector for founder existence for B[c,m] = 1, 2, 3. (shared/merge/unique)
    #
    #     α_full::Array{R, 1}
    #     α_back::Array{R, 1}
    #     α_uni::Array{R, 1}
    #     α_err::Array{R, 1}
    #
    #     α_t::Array{R, 1}
    #
    #     ln_p_main_penalty::Array{R, 1} # main block state penalty, [ln_penalty, ln (1-penalty)
    #     ln_p_err_penalty::Array{R, 1} # error block state penalty, [ln_penalty, ln (1-penalty)
    #
    #     p_hap::R
    # end
    #
    # function parse_config_file(config_file::String,
    #                            INT::Type{<:Integer} = Int32,
    #                            REAL::Type{<:Real}   = Float32)::Tuple{Parameters{INT, REAL}, Annealer{INT, REAL}}
    #     conf = ConfParse(config_file)
    #     parse_conf!(conf)
    #     println(conf)
    #     α_s   = parse(REAL, String(retrieve(conf, "model", "alpha_s")) )
    #     α_v   = parse(REAL, String(retrieve(conf, "model", "alpha_v")) )
    #
    #     δ_f   = parse(REAL, String(retrieve(conf, "model", "delta_f")) )
    #     δ_b   = parse(REAL, String(retrieve(conf, "model", "delta_b")) )
    #
    #     α_full_0   = parse(REAL, String(retrieve(conf, "model", "alpha_full_0")) )
    #     α_full_1   = parse(REAL, String(retrieve(conf, "model", "alpha_full_1")) )
    #
    #     α_back_0   = parse(REAL, String(retrieve(conf, "model", "alpha_back_0")) )
    #     α_back_1   = parse(REAL, String(retrieve(conf, "model", "alpha_back_1")) )
    #
    #     α_uni_0   = parse(REAL, String(retrieve(conf, "model", "alpha_uni_0")) )
    #     α_uni_1   = parse(REAL, String(retrieve(conf, "model", "alpha_uni_1")) )
    #
    #     α_err_0   = parse(REAL, String(retrieve(conf, "model", "alpha_err_0")) )
    #     α_err_1   = parse(REAL, String(retrieve(conf, "model", "alpha_err_1")) )
    #
    #     α_t_0   = parse(REAL, String(retrieve(conf, "model", "alpha_t_0")) )
    #     α_t_1   = parse(REAL, String(retrieve(conf, "model", "alpha_t_1")) )
    #     α_t_2   = parse(REAL, String(retrieve(conf, "model", "alpha_t_2")) )
    #
    #     p_hap     = parse(REAL, String(retrieve(conf, "model", "p_hap")) )
    #
    #     ln_p_v_main    = parse(REAL, String(retrieve(conf, "model", "ln_main_penalty")) )
    #     ln_1m_p_v_main = log(ℯ, (REAL)(1.0) - exp(ln_p_v_main))
    #
    #     ln_p_v_err     = parse(REAL, String(retrieve(conf, "model", "ln_err_penalty")) )
    #     ln_1m_p_v_err  = log(ℯ, (REAL)(1.0) - exp(ln_p_v_err))
    #
    #     period    = parse(REAL, String(retrieve(conf, "annealer", "period")) )
    #
    #     ln_p_ladders_main_strings = (retrieve(conf, "annealer", "ln_p_ladders_main"))
    #     ln_p_ladders_main = zeros(REAL, length(ln_p_ladders_main_strings), 2)
    #     for i in (INT)(1):(INT)(length(ln_p_ladders_main_strings))
    #         step = String(ln_p_ladders_main_strings[i])
    #         print("step: "); println(step)
    #         ln_1m_p_step = parse(REAL, step)
    #         ln_p_step    = log(ℯ, (REAL)(1.0) - exp(ln_1m_p_step))
    #         ln_p_ladders_main[i, 1] = ln_1m_p_step
    #         ln_p_ladders_main[i, 2] = ln_p_step
    #     end
    #
    #     ln_p_ladders_err_strings = (retrieve(conf, "annealer", "ln_p_ladders_err"))
    #     ln_p_ladders_err = zeros(REAL, length(ln_p_ladders_err_strings), 2)
    #     for i in (INT)(1):(INT)(length(ln_p_ladders_err_strings))
    #         step = String(ln_p_ladders_err_strings[i])
    #         print("step: "); println(step)
    #         ln_1m_p_step = parse(REAL, step)
    #         ln_p_step    = log(ℯ, (REAL)(1.0) - exp(ln_1m_p_step))
    #         ln_p_ladders_err[i, 1] = ln_1m_p_step
    #         ln_p_ladders_err[i, 2] = ln_p_step
    #     end
    #
    #     # @assert (length(ln_p_ladders_main_strings) = length(ln_p_ladders_err_strings))
    #     T::INT = (INT)(length(ln_p_ladders_main_strings))
    #
    #     p_temperature::Array{REAL, 1} = zeros(REAL, T)
    #     p_temperature .= (REAL)( 1.0 / T )
    #
    #     q_temperature::Array{REAL, 2} = zeros(REAL, T, 3)
    #     for i in 2:(T-1)
    #         for j in 1:3
    #             q_temperature[i, j] = (REAL)(1.0/3.0)
    #         end
    #     end
    #     q_temperature[1,1] = 0.0; q_temperature[T,1] = 0.5
    #     q_temperature[1,2] = 0.5; q_temperature[T,2] = 0.5
    #     q_temperature[1,3] = 0.5; q_temperature[T,3] = 0.0
    #
    #     return (Parameters{INT,REAL}(α_s, α_v,
    #                                 [δ_f, δ_b],
    #                                 [α_full_0, α_full_1],
    #                                 [α_back_0, α_back_1],
    #                                 [α_uni_0,  α_uni_1],
    #                                 [α_err_0,  α_err_1],
    #                                 [α_t_0,  α_t_1, α_t_2],
    #                                 [ln_p_v_main, ln_1m_p_v_main],
    #                                 [ln_p_v_err,  ln_1m_p_v_err],
    #                                 p_hap),
    #            Annealer{INT, REAL}(T, ln_p_ladders_main, ln_p_ladders_err, p_temperature, q_temperature, 1))
    # end


end
