
module Diagnostics
    using Distributions
    using Statistics
    import StatsBase: autocov
    using PyPlot
    # autocor, countmap, counts, describe, predict, quantile, sample, sem, summarystats

    function mcse(x::Vector{T}, method::Symbol=:imse; args...) where {T<:Real}
      method == :bm ? mcse_bm(x; args...) :
      method == :imse ? mcse_imse(x) :
      method == :ipse ? mcse_ipse(x) :
        throw(ArgumentError("unsupported mcse method $method"))
    end

    function mcse_bm(x::Vector{T}; size::Integer=100) where {T<:Real}
      n = length(x)
      m = div(n, size)
      m >= 2 ||
        throw(ArgumentError(
          "iterations are < $(2 * size) and batch size is > $(div(n, 2))"
        ))
      mbar = [mean(x[i * size .+ (1:size)]) for i in 0:(m - 1)]
      sem(mbar)
    end

    function mcse_imse(x::Vector{T}) where {T<:Real}
      n = length(x)
      m = div(n - 2, 2)
      ghat = autocov(x, [0, 1])
      Ghat = sum(ghat)
      value = -ghat[1] + 2 * Ghat
      for i in 1:m
        Ghat = min(Ghat, sum(autocov(x, [2 * i, 2 * i + 1])))
        Ghat > 0 || break
        value += 2 * Ghat
      end
      sqrt(value / n)
    end

    function mcse_ipse(x::Vector{T}) where {T<:Real}
      n = length(x)
      m = div(n - 2, 2)
      ghat = autocov(x, [0, 1])
      value = ghat[1] + 2 * ghat[2]
      for i in 1:m
        Ghat = sum(autocov(x, [2 * i, 2 * i + 1]))
        Ghat > 0 || break
        value += 2 * Ghat
      end
      sqrt(value / n)
    end

    function geweke(value_list::Vector{T}, burnin; former::Float64 = 0.1, latter::Float64 = 0.5, method::Symbol=:imse) where {T<:Real}
        @assert length(value_list) > burnin

        chain::Array{Float64, 1} = convert.(Float64, value_list)

        former_avg::Float64 = (Float64)(0.0)
        former_se::Float64  = (Float64)(0.0)

        # former_start = (Int64)(length(chain) - burnin)
        # former_end   = former_start + (Int64)(round((length(chain) - burnin) * former))
        # latter_start = (Int64)(round((length(chain) - burnin) * latter ))
        # latter_end   = (Int64)(round((length(chain) - burnin)))

        former_start = (Int64)(burnin + 1)
        former_end   = (Int64)(burnin + (round((length(chain) - burnin) * former)))
        latter_start = (Int64)(burnin + (round((length(chain) - burnin) * latter )))
        latter_end   = (Int64)(length(chain))

        println("former_start, former_end, latter_start, latter_end")
        println(former_start)
        println(former_end)
        println(latter_start)
        println(latter_end)

        println("chain[former_start:former_end]")
        println(chain[former_start:former_end])
        former_avg   = (Float64)(mean(chain[former_start:former_end]))
        former_mcse  = (Float64)(mcse(chain[former_start:former_end], method = method))

        println("chain[latter_start:latter_end]")
        println(chain[latter_start:latter_end])
        latter_avg   = (Float64)(mean(chain[latter_start:latter_end]))
        latter_mcse  = (Float64)(mcse(chain[latter_start:latter_end], method = method))

        # TODO: _var function is not appropriate. Use monte calro variance estimeters.
        Z::Float64   = (former_avg - latter_avg) / (Float64)(sqrt(former_mcse^2 + latter_mcse^2))
        println("Z:" * string(Z))
        Norm = Distributions.Normal(0.0, 1.0)

        return 1.0 - cdf(Norm, abs(Z))
    end


    function viewtraced(trace_list, output_dir::String, variable_symbol::String, burnin, acf_length, font_size = 30)
        # Plots.PyPlot.rc("font", size = 30.0)
        PyPlot.rc("agg.path", chunksize = 100000 )
        PyPlot.rc("font", size = font_size )

        if typeof(variable_symbol) != Symbol
            variable_symbol = Symbol(variable_symbol)
        end

        @assert length(trace_list) > 0
        @assert typeof(getproperty(trace_list[1], variable_symbol)) <: Real
        @assert acf_length > 0  && typeof(acf_length) <: Integer

        acf_length = min(acf_length, length(trace_list))
        if acf_length <= 0
          println("acf_length <= 0 @ function viewtraced(trace_list, output_dir::String, variable_symbol::String, burnin, acf_length) in Diagnostics.jl")
          return
        end

        item_list::Array{Real, 1} = []
        for i in 1:length(trace_list)
            @assert typeof(getproperty(trace_list[i], variable_symbol)) <: Real
            push!(item_list, getproperty(trace_list[i], variable_symbol))
        end

        iters = collect(1:length(item_list))
        PyPlot.clf()
        PyPlot.plot(iters, item_list)
        # PyPlot.title("Geweke:" * String(variable_symbol)* ":" * string(geweke_ans))
        PyPlot.title(String(variable_symbol))
        PyPlot.savefig(output_dir * "/trace_" * String(variable_symbol) * ".png")

        geweke_ans = geweke(item_list, burnin, former = 0.1, latter = 0.5)

        # iters = collect(1:length(item_list))
        # PyPlot.clf()
        # PyPlot.plot(iters, item_list)
        # # PyPlot.title("Geweke:" * String(variable_symbol)* ":" * string(geweke_ans))
        # PyPlot.title(String(variable_symbol))
        # PyPlot.savefig(output_dir * "/trace_" * String(variable_symbol) * ".png")

        # plot(iters, item_list, size=(2000,2000), title = "Geweke:" * String(variable_symbol)* ":" * string(geweke_ans))
        # savefig(output_dir * "/trace_" * String(variable_symbol) * ".png")

        PyPlot.clf()
        acf = _acf(item_list, acf_length)
        ess = _ess(acf, length(item_list))
        PyPlot.bar((1):(length(acf)), acf)
        PyPlot.title("Geweke:" * String(variable_symbol)* ":" * string(geweke_ans)  * ", ESS: " * string(ess))
        PyPlot.savefig(output_dir * "/acf_" * String(variable_symbol) * ".png")
    end

    function _acf(s::Array{R, 1}, K::I)::Array{R, 1} where { R <: Real, I <:Integer }
        S = length(s)
        ans::Array{R, 1} = zeros(R, min((S-1), K))
        s_avg::R = sum(s) / S
        s_var::R = sum( (s .- s_avg) .* (s .- s_avg) ) / (R)(S)

        for t in 1:min((S-1), K)
            corr::R = 0.0
            for i in 1:(S-t)
                corr += (s[i] - s_avg) * (s[i+t] - s_avg)
            end
            corr = corr / (R)(S)
            ans[t] = corr / s_var
        end
        return ans
    end

    function _ess(acf::Array{R, 1}, N::I)::R where { I <: Integer, R <: Real }
        println("== ess ==")
        print("N: ")
        println(N)
        print("1 + 2 * sum(acf): ")
        println(1 + 2 * sum(acf))
        print("1 + 2 * sum(acf): ")
        println(1 + 2 * sum(acf))
        print("N /( 1 + 2 * sum(acf) ): ")
        println((R)( N / ( 1 + 2 * sum(acf))))
        return (R)( N / ( 1 + 2 * sum(acf)))
    end

end
