@Include "InputParser.jl"
@Include "RandomUtil.jl"

module SamplerType
    using ..InputParser
    using ..RandomUtil
    using Random
    using Distributions
    using SpecialFunctions

    # TODO: re-define parameters...
    mutable struct Sampler{I,R}
        T::Array{I, 1}          # Somatic state, S[j] == 1 error and S[j] > 1 somatic
        Z::Array{I, 2}          # Z == 1 err, Z == 2 mutation

        S::Array{I, 1}          # Somatic state, S[j] == 1 error and S[j] > 1 somatic
        Y::Array{I, 2}          # Z == 1 err, Z == 2 mutation

        ln_p_evidence::Array{R, 3}
        ln_p_data::Array{R, 3}

        p_consistent::R

        p_fp_data::R
        p_fp_evidence::R

        ln_θ::R
        ln_rho::R
        ln_prob::R

        count_Y::Array{I, 3}
    end
    export Sampler

    mutable struct StateSummary{I,R}
        ln_prob::R
        sum_t::I
        sum_z::I
        sum_s::I
        sum_y::I
    end
    export StateSummary

    function convert_to_summary(samp::Sampler{I, R})::StateSummary{I,R} where {I <: Integer, R <: Real}
        return StateSummary(samp.ln_prob, sum(samp.T .- (I)(1)), sum(samp.Z .- (I)(1)), sum(samp.S .- (I)(1)), sum(samp.Y .- (I)(1)) )
    end

    function _ln_p_evidence_z(samp::Sampler{I,R})::Array{R, 1} where {I<:Integer, R<:Real}
        N::I, M::I = size(samp.Z)
        ans_ln::Array{R, 1} = convert.(R, zeros(M))
        for j in (I)(1):M
            for i in (I)(1):N
                ans_ln[j] += samp.ln_p_data[i, j, samp.Z[i, j]]
            end
        end
        return ans_ln
    end

    function _ln_p_data_y_s_t(samp::Sampler{I,R})::Array{R, 1} where {I<:Integer, R<:Real}
        N::I, M::I = size(samp.Z)
        ans_ln::Array{R, 1} = convert.(R, zeros(M))
        for j in (I)(1):M
            for i in (I)(1):N
                ans_ln[j] += samp.ln_p_data[i, j, samp.Y[i, j]]
                ans_ln[j] += (samp.Y[i, j] == 2) * (samp.S[j] == 2) * samp.ln_rho
                ans_ln[j] += (samp.Y[i, j] == 2) * (samp.T[j] == 2) * samp.ln_θ
            end
        end
        return ans_ln
    end

    function _ln_p_y(samp::Sampler{I,R})::Array{R, 1} where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ans::R = (R)(0.0)
        ans_ln::Array{R, 1} = convert.(R, zeros(M))
        for j in (I)(1):M
            p::R = (R)( (samp.S[j] == 1) * samp.p_fp_data + (samp.S[j] == 2) * 0.5)
            for i in (I)(1):N
                ans_ln[j] += _ln_p_ber(samp.Z[i, j]-(I)(1), p)
            end
        end
        return ans_ln
    end

    function _ln_p_z(samp::Sampler{I,R})::Array{R, 1} where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ans::R = (R)(0.0)
        ans_ln::Array{R, 1} = convert.(R, zeros(M))
        for j in (I)(1):M
            p::R = (R)( (samp.T[j] == 1) * samp.p_fp_evidence + (samp.T[j] == 2) * 0.5)
            for i in (I)(1):N
                ans_ln[j] += _ln_p_ber(samp.Z[i, j]-(I)(1), p)
            end
        end
        return ans_ln
    end

    function _ln_p_s_t(samp::Sampler{I, R})::Array{R, 1} where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ans_ln::Array{R, 1} = convert.(R, zeros(M))
        for j in (I)(1):M
            p::R = (R)( (samp.T[j] == 1) * 0.5 + (samp.T[j] == 2) * samp.p_consistent)
            ans_ln[j] += _ln_p_ber( (I)(samp.S[j] > (I)(1)), p)
        end
        return ans_ln
    end

    function _ln_p_t(samp::Sampler{I, R})::Array{R, 1} where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ans_ln::Array{R, 1} = convert.(R, zeros(M))
        for j in (I)(1):M
            ans_ln[j] += _ln_p_ber( (I)(samp.S[j] > (I)(1)), (R)(0.5))
        end
        return ans_ln
    end


    function ln_p_all(samp::Sampler{I,R}, debug::Bool = false)::R where {I<:Integer, R<:Real}
        ans::R = (R)(0.0)
        S::I, M::I = size(samp.Z)
        ans += sum(_ln_p_evidence_z(samp))
        ans += sum(_ln_p_data_y_s_t(samp))
        ans += sum(_ln_p_y(samp))
        ans += sum(_ln_p_z(samp))

        ans += sum(_ln_p_s_t(samp))
        ans += sum(_ln_p_t(samp))
        return ans
    end

    function ln_p_all_array(samp::Sampler{I,R}, debug::Bool = false)::Array{R, 1} where {I<:Integer, R<:Real}
        return _ln_p_evidence_z(samp) + _ln_p_data_y_s_t(samp) + _ln_p_y(samp) + _ln_p_z(samp) + _ln_p_s_t(samp) + _ln_p_t(samp)
    end

    function sample_t!(samp::Sampler{I,R})::Nothing where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ln_p_fp_evidence::Array{R, 1} = [ (R)(log(ℯ, 1.0-samp.p_fp_evidence)), (R)(log(ℯ, samp.p_fp_evidence))]
        ln_p_consistent::Array{R, 1} = [ (R)(log(ℯ, 1.0-samp.p_consistent)), (R)(log(ℯ, samp.p_consistent))]
        ln_p_half::R = (R)(log(ℯ, 0.5))

        ln_p::Array{R, 1} = [(R)(0.0), (R)(0.0)]
        for j in (I)(1):M
            ln_p    .= (R)(0.0)
            ln_p[1] += ln_p_half
            ln_p[2] += (R)((samp.S[j] == 1) * ln_p_consistent[1] + (samp.S[j] == 2) * ln_p_consistent[2])
            for i in (I)(1):N
                ln_p[1] += (R)((samp.Z[i, j] == 1) * ln_p_fp_evidence[1] + (samp.Z[i, j] == 2) * ln_p_fp_evidence[2])
                ln_p[2] += ln_p_half + (samp.Y[i,j] == (I)(2)) * samp.ln_θ
            end
            _exp_normalize!(ln_p)
            samp.T[j] = (I)(argmax( RandomUtil.sample_multinomial((I)(1), ln_p) ))
        end
        return nothing
    end

    function sample_z!(samp::Sampler{I,R})::Nothing where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ln_p::Array{R, 1} = [(R)(0.0), (R)(0.0)]
        ln_p_fp_evidence::Array{R, 1} = [ (R)(log(ℯ, 1.0-samp.p_fp_evidence)), (R)(log(ℯ, samp.p_fp_evidence))]
        # ln_p_half::R = (R)(log(ℯ, 0.5))
        for j in (I)(1):M
            for i in (I)(1):N
                ln_p .= (R)(0.0)
                ln_p[1] += (samp.T[j] == 1) * ln_p_fp_evidence[1]
                ln_p[2] += (samp.T[j] == 1) * ln_p_fp_evidence[2]
                ln_p[2] += samp.ln_p_evidence[i, j, 2]
                _exp_normalize!(ln_p)
                samp.Z[i, j] = (I)(argmax( RandomUtil.sample_multinomial((I)(1), ln_p) ))
            end
        end
        return nothing
    end

    function sample_s!(samp::Sampler{I,R})::Nothing where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)

        ln_p_fp_data::Array{R, 1} = [ (R)(log(ℯ, 1.0-samp.p_fp_data)), (R)(log(ℯ, samp.p_fp_data))]
        ln_p_consistent::Array{R, 1} = [ (R)(log(ℯ, 1.0-samp.p_consistent)), (R)(log(ℯ, samp.p_consistent))]
        ln_p_half::R = (R)(log(ℯ, 0.5))

        ln_p::Array{R, 1} = [(R)(0.0), (R)(0.0)]
        for j in (I)(1):M
            ln_p    .= (R)(0.0)
            ln_p[1] += (R)( (samp.T[j] == 2) * ln_p_consistent[1] )
            ln_p[2] += (R)( (samp.T[j] == 2) * ln_p_consistent[2] )
            for i in (I)(1):N
                ln_p[1] += (R)( (samp.Y[i,j] == (I)(1)) * ln_p_fp_data[1] + (samp.Y[i,j] == (I)(2)) * ln_p_fp_data[2] )
                ln_p[2] += ln_p_half + (R)((samp.Y[i,j] == (I)(2)) * samp.ln_rho)
            end
            _exp_normalize!(ln_p)
            samp.S[j] = (I)(argmax( RandomUtil.sample_multinomial((I)(1), ln_p) ))
        end
        return nothing
    end

    function sample_y!(samp::Sampler{I,R})::Nothing where {I <: Integer, R <: Real}
        N::I, M::I = size(samp.Z)
        ln_p::Array{R, 1} = [(R)(0.0), (R)(0.0)]
        ln_p_fp_data::Array{R, 1} = [ (R)(log(ℯ, 1.0-samp.p_fp_data)), (R)(log(ℯ, samp.p_fp_data))]
        ln_p_half::R = (R)(log(ℯ, 0.5))

        for j in (I)(1):M
            for i in (I)(1):N
                ln_p .= (R)(0.0)
                ln_p[1] += (R)((samp.S[j] == 1) * ln_p_fp_data[1])
                ln_p[2] += (R)((samp.S[j] == 1) * ln_p_fp_data[2])
                ln_p[2] += (R)(samp.ln_p_data[i, j, 2] + (samp.T[j] == 2) * samp.ln_θ + (samp.S[j] == 2) * samp.ln_rho)
                _exp_normalize!(ln_p)
                samp.Y[i, j] = (I)(argmax( RandomUtil.sample_multinomial((I)(1), ln_p) ))
            end
        end
        return nothing
    end

    function sample_map_y!(samp::Sampler{I, R},
                           seed::I = (I)(0),
                           iter::I = (I)(100000),
                           thin::I = (I)(1),
                           burnin::I = (I)(0))::Tuple{Sampler{I, R}, Array{StateSummary{I, R}, 1}} where {I <:Integer, R <: Real}
        Random.seed!(seed)
        S::I, M::I = size(samp.Z)
        traced_variables::Array{StateSummary{I, R}, 1} = []

        map_state::Sampler{I, R} = deepcopy(samp)
        max_ln_prob::Array{R, 1} = convert.(R, zeros(M))
        max_ln_prob .= (R)(-Inf)
        ln_p_now::Array{R, 1} = convert.(R, zeros(M))

        for count in (I)(1):(iter+burnin)
            sample_t!(samp)
            sample_z!(samp)
            sample_s!(samp)
            sample_y!(samp)

            samp.ln_prob = ln_p_all(samp)
            if (count <= burnin) || (count % thin == 0)
                for i in (I)(1):S
                    for j in (I)(1):M
                        samp.count_Y[i,j,1] += (I)( samp.Y[i,j] == 1 )
                        samp.count_Y[i,j,2] += (I)( samp.Y[i,j] == 2 )
                    end
                end
                push!(traced_variables, convert_to_summary(samp))

                ln_p_now .= ln_p_all_array(samp)
                for j in (I)(1):M
                    if ln_p_now[j] > max_ln_prob[j]
                        map_state.T[j] = samp.T[j]
                        map_state.Z[:, j] .= samp.Z[:, j]
                        map_state.S[j] = samp.S[j]
                        map_state.Y[:, j] .= samp.Y[:, j]
                        max_ln_prob[j] = ln_p_now[j]
                    end
                end
            end
            if count % 100 == 0
                println(count)
                println(samp.ln_prob)
            end
        end
        print("samp.p_consistent: ")
        println(samp.p_consistent)
        print("samp.p_fp_evidence: ")
        println(samp.p_fp_evidence)
        print("samp.p_fp_data: ")
        println(samp.p_fp_data)
        print("samp.ln_θ: ")
        println(samp.ln_θ)
        print("samp.ln_rho: ")
        println(samp.ln_rho)
        print("S: ")
        println(S)
        print("M: ")
        println(M)
        # for j in (I)(1):M
        #     for i in (I)(1):S
        #         map_state.Y[i,j] = (I)(samp.count_Y[i,j,1] < samp.count_Y[i,j,2]) + (I)(1)
        #     end
        # end
        return map_state, traced_variables
    end

    # # TODO: re-define parameters...
    # mutable struct Sampler{I,R}
    #     T::Array{I, 1}          # Somatic state, S[j] == 1 error and S[j] > 1 somatic
    #     Z::Array{I, 2}          # Z == 1 err, Z == 2 mutation
    #
    #     S::Array{I, 1}          # Somatic state, S[j] == 1 error and S[j] > 1 somatic
    #     Y::Array{I, 2}          # Z == 1 err, Z == 2 mutation
    #
    #     ln_p_evidence::Array{R, 3}
    #     ln_p_data::Array{R, 3}
    #
    #     p_consistent::R
    #
    #     p_fp_evidence::R
    #     p_fp_data::R
    #
    #     ln_θ::R
    #     ln_rho::R
    #     ln_prob::R
    # end

    # TODO
    function init(data_score_path::String, evidence_score_path::String, p_fp_data, p_fp_evidence, p_consistent, log_10_θ, log_10_rho,
                  transpose::Bool = true,
                  INT::Type{<:Integer} = Int32, REAL::Type{<:Real} = Float32)
        ln_p_e::Array{REAL, 3} = _parseData(evidence_score_path, transpose, INT, REAL)
        ln_p_d::Array{REAL, 3} = _parseData(data_score_path, transpose, INT, REAL)

        N::INT = size(ln_p_d, 1)
        M::INT = size(ln_p_d, 2)

        T::Array{INT, 1}        = convert.(INT, fill(1,M))
        Z::Array{INT, 2}        = convert.(INT, fill(1,N,M)) # init ℤ, Z[i,j] ∈ {1,2}, 1: error, 2: tumor
        S::Array{INT, 1}        = convert.(INT, fill(1,M))
        Y::Array{INT, 2}        = convert.(INT, fill(1,N,M)) # init ℤ, Z[i,j] ∈ {1,2}, 1: error, 2: tumor
        count_Y::Array{INT, 3}  = convert.(INT, fill(0,N,M,2))

        ln_θ::REAL   = log_10_θ * log(ℯ, 10)
        ln_rho::REAL = log_10_rho * log(ℯ, 10)
        samp::Sampler{INT, REAL} = Sampler{INT, REAL}(T, Z, S, Y,  ln_p_e, ln_p_d,
                                                      p_consistent, p_fp_data, p_fp_evidence,
                                                      (REAL)(ln_θ), (REAL)(ln_rho),
                                                      (REAL)(0.0), count_Y)
        samp.ln_prob = ln_p_all(samp)
        print("samp.p_consistent: ")
        println(samp.p_consistent)
        print("samp.p_fp_evidence: ")
        println(samp.p_fp_evidence)
        print("samp.p_fp_data: ")
        println(samp.p_fp_data)
        print("samp.ln_θ: ")
        println(samp.ln_θ)
        print("samp.ln_rho: ")
        println(samp.ln_rho)
        print("N: ")
        println(N)
        print("M: ")
        println(M)
        return samp
    end

    function evaluate_softmax_loss(samp::Sampler{I, R}, summary_score::Array{R, 2}) where {I <: Integer, R <: Real}
        N::I, M::I = size(summary_score)
        loss::R = (R)(0.0)

        ln_posterior_y::Array{R, 3} = convert.(R, zeros(N, M, 2))
        for j in (I)(1):M
            for i in (I)(1):N
                ln_posterior_y[i, j, 1] = (R)((samp.count_Y[i, j, 1] + 1)/ (samp.count_Y[i, j, 1] + samp.count_Y[i, j, 2] + 2))
                ln_posterior_y[i, j, 2] = (R)((samp.count_Y[i, j, 2] + 1)/ (samp.count_Y[i, j, 1] + samp.count_Y[i, j, 2] + 2))
            end
        end
        ln_posterior_y .= log.(ℯ, ln_posterior_y)
        for j in (I)(1):M
            for i in (I)(1):N
                loss -= (R)( (1.0 - summary_score[i, j]) * ln_posterior_y[i, j, 1] + (summary_score[i, j]) * ln_posterior_y[i, j, 2] )
            end
        end
        return loss
    end

    function exec_train(param_path::String, ans_path::String, data_score_path::String, evidence_score_path::String,
                        p_fp_data, p_fp_evidence, p_consistent, log_10_θ, log_10_rho,
                        transpose::Bool = false,
                        seed::INT = (I)(0), iter::INT = (I)(100000), thin::INT = (I)(10), burnin::INT = (I)(10),
                        REAL::Type{<:Real} = Float32) where {INT <: Integer}

        summary_score::Array{REAL, 2} = InputParser.parse_input_summary(ans_path, transpose, INT, REAL)

        min_loss::REAL       = (REAL)(Inf)
        min_log_10_rho::REAL = (REAL)(0.0)
        min_log_10_θ::REAL   = (REAL)(0.0)

        open(param_path, "r") do param_file
            for line::String in eachline(param_file)
                if startswith(line, "#")
                    continue
                end
                line = strip(line)
                line_col::Array{String,1} = split(line, '\t')

                log_10_rho = parse(REAL, line_col[1])
                log_10_θ   = parse(REAL, line_col[2])
                samp = SamplerType.init(data_score_path, evidence_score_path, p_fp_data, p_fp_evidence, p_consistent, log_10_θ, log_10_rho, transpose, INT, REAL)

                map, traced = SamplerType.sample_map_y!(samp, seed, iter, thin, burnin)

                loss::REAL = evaluate_softmax_loss(samp, summary_score)

                println("########################################")
                print("loss: "); println(loss)
                print("log_10_rho: "); println(log_10_rho)
                print("log_10_θ: "); println(log_10_θ)
                println("########################################")

                if loss < min_loss
                    min_loss = loss
                    min_log_10_rho = log_10_rho
                    min_log_10_θ   = log_10_θ
                end
            end
        end

        println("########################################")
        print("min_loss: "); println(min_loss)
        print("min_log_10_rho: "); println(min_log_10_rho)
        print("min_log_10_θ: "); println(min_log_10_θ)
        println("########################################")

        return (min_log_10_rho, min_log_10_θ, min_loss)
    end

    function exec_map(data_score_path::String, evidence_score_path::String, p_fp_data, p_fp_evidence, p_consistent, log_10_θ, log_10_rho,
                      transpose::Bool = false,
                      seed::INT = (I)(0), iter::INT = (I)(100000), thin::INT = (I)(10), burnin::INT = (I)(10),
                      REAL::Type{<:Real} = Float32) where {INT <: Integer}
        samp = SamplerType.init(data_score_path, evidence_score_path, p_fp_data, p_fp_evidence, p_consistent, log_10_θ, log_10_rho, transpose, INT, REAL)
        map, traced = SamplerType.sample_map_y!(samp, seed, iter, thin, burnin)
        return (map, traced)
    end

    function _ping_sampler(err_scores::String, mat_scores::String, pat_scores::String, ini_file::String)
        samp = SamplerType.init(err_scores, mat_scores, pat_scores, ini_file)
        sampled = SamplerType.sample_all!(samp)
        return sampled
    end

    function __ping_sampler()
        samp = SamplerType.init("../../simulationTree/err.score.txt", "../../simulationTree/mat.score.txt",
                            "../../simulationTree/pat.score.txt", "./simpleModel.ini")
        sampled = SamplerType.sample_all!(samp)
        return sampled
    end

    #
    # private funcitons
    #
    function _parseData(mut_score_path::String,
                        transpose::Bool = true,
                        INT::Type{<:Integer} = Int32, REAL::Type{<:Real}=Float32)::Array{REAL, 3}
        pat_score::Array{REAL, 2} = InputParser.parse_input_summary(mut_score_path, transpose, INT, REAL)
        S::INT, M::INT = size(pat_score)
        err_score::Array{REAL, 2} = convert.(REAL, zeros(S,M))

        ln_p::Array{REAL, 3} = convert.(REAL, zeros(S,M,2))
        for (s,m) in Iterators.product((INT)(1):S,(INT)(1):M)
            ln_p[s,m,1] = err_score[s,m]
            ln_p[s,m,2] = pat_score[s,m]
        end
        return ln_p
    end

    function _exp_normalize!(ln_p::Array{R,1})::Nothing where {R <: Real}
        ln_p .= (ln_p .- maximum(ln_p))
        ln_p .= exp.(ln_p)
        ln_p .= ln_p ./ sum(ln_p)
        return nothing
    end

    function _exp_normalize!(ln_p::AbstractArray{R,1})::Nothing where {R <: Real}
        ln_p .= (ln_p .- maximum(ln_p))
        ln_p .= exp.(ln_p)
        ln_p .= ln_p ./ sum(ln_p)
        return nothing
    end

    function _log_sum_exp(ln_p::Array{R,1})::R where {R <: Real}
        maxVal::R = maximum(ln_p)
        return log(ℯ, sum(exp.( (ln_p .- maxVal) ))) + maxVal
    end

    # x ∈ {0,1}, x = 1 w.p. p
    function _ln_p_ber(x::I, p::R)::R where {R <: Real, I <: Integer}
        return x * log(ℯ, p) + (1-x) * log(ℯ, (1.0-p))
    end

    function _ln_p_beta(p::R, α::R, β::R)::R where {R <: Real}
        return (α-(R)(1.0))*log(ℯ, p) + (β-(R)(1.0))*log(ℯ,((R)(1.0)-p)) - SpecialFunctions.lbeta(α, β)
    end

    function _ln_p_dir(p::Array{R, 1}, α::Array{R, 1})::R where {R <: Real}
        ans::R = 0.0
        @assert length(p) == length(α)
        for i in 1:length(p)
            ans += (α[i]-(R)(1.0))*log(ℯ, p[i])
            ans -= SpecialFunctions.lgamma(α[i])
        end
        ans += SpecialFunctions.lgamma(sum(α))
        return ans
    end

end
