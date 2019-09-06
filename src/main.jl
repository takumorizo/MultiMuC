using Pkg
Pkg.activate(@__DIR__)
module MultiMuC
    include(string(@__DIR__) * "/include.jl")
    @Include "Result.jl"
    @Include "SamplerType.jl"
    @Include "Diagnostics.jl"

    using .SamplerType
    using .Diagnostics
    using .Result
    using DocOpt
    using Dates

    export pingResult

    function run(INT::Type{<:Integer} = Int32, REAL::Type{<:Real} = Float32)
        nowTime = (Dates.value(Dates.now())) #Int(Dates.now())
        cwd = pwd()
        doc = """majority_voting

        Usage:
            main_simple.jl MAP <data> <evidence> [options]
            main_simple.jl TRAIN <param> <ans> <result> <data> <evidence> [options]

        Options:
          -o <DIR> --outDir=<DIR>  Output directory. We use current working directory if unspecified. [default: $cwd]
          -s <SEED> --seed=<SEED>  Seed of randomness. We use current time (msec) if unspecified. [default: $nowTime]
          -n <NUM> --number=<NUM>  Number of iteration after burnin. [default: 1000]
          -b <BURNIN> --burnin=<BURNIN>  Number of iteration during burnin. [default: 0]
          -t <THIN> --thin=<THIN>  Duration between sampling. [default: 1]
          -m <THRES> --threshold=<THRES>  Threshold value for call mutation from file. [default: 0.0]
          --pdata=<FPDATA>  Prob for false positive in data. [default: 0.2]
          --pevidence=<FPEVIDENCE>  Prob for false positive in evidence. [default: 0.02]
          --pcon=<CONSISTENCE>  Prob for consistency between evidence and data. [default: 0.999]
          --theta=<THETA>  Down threshold value from evidence in log_10 form. [default: 0.5]
          --rho=<RHO>  Down threshold value from data in log_10 form. [default: 0.1]
          --transpose  Transpose input matrix. If transpose, input col is sample id and node id, otherwise col is mutation.
          -h --help  Show this screen.

        """
        # println(doc)
        args = docopt(doc)
        println(args)

        if args["MAP"] == "MAP"
            @time smax, traced = SamplerType.exec_map(args["<data>"], args["<evidence>"],
                                                      parse(REAL, args["--pdata"]),
                                                      parse(REAL, args["--pevidence"]),
                                                      parse(REAL, args["--pcon"]),
                                                      parse(REAL, args["--theta"]),
                                                      parse(REAL, args["--rho"]),
                                                      args["--transpose"],
                                                      parse(INT, args["--seed"]),
                                                      parse(INT, args["--number"]),
                                                      parse(INT, args["--thin"]),
                                                      parse(INT, args["--burnin"]),
                                                      REAL)
            Result.print_map_z_y(smax, args["--outDir"], args["<data>"])
            Diagnostics.viewtraced(traced, args["--outDir"], "sum_s",   parse(INT, args["--burnin"]), min(100, length(traced)), 10)
            Diagnostics.viewtraced(traced, args["--outDir"], "sum_y",   parse(INT, args["--burnin"]), min(100, length(traced)), 10)

            Diagnostics.viewtraced(traced, args["--outDir"], "sum_t",   parse(INT, args["--burnin"]), min(100, length(traced)), 10)
            Diagnostics.viewtraced(traced, args["--outDir"], "sum_z",   parse(INT, args["--burnin"]), min(100, length(traced)), 10)
            Diagnostics.viewtraced(traced, args["--outDir"], "ln_prob", parse(INT, args["--burnin"]), min(100, length(traced)), 10)
        elseif args["TRAIN"] == "TRAIN"
            @time min_log_10_rho, min_log_10_θ, min_loss = SamplerType.exec_train(args["<param>"], args["<ans>"],
                                                                                  args["<data>"], args["<evidence>"],
                                                                                  parse(REAL, args["--pdata"]),
                                                                                  parse(REAL, args["--pevidence"]),
                                                                                  parse(REAL, args["--pcon"]),
                                                                                  parse(REAL, args["--theta"]),
                                                                                  parse(REAL, args["--rho"]),
                                                                                  args["--transpose"],
                                                                                  parse(INT, args["--seed"]),
                                                                                  parse(INT, args["--number"]),
                                                                                  parse(INT, args["--thin"]),
                                                                                  parse(INT, args["--burnin"]),
                                                                                  REAL)

            open( args["<result>"], "w" ) do fp
                write(fp, "#min_log_10_rho" * "\t" * "min_log_10_θ" * "\t" * "min_loss" * "\n")
                write(fp, string(min_log_10_rho) * "\t" * string(min_log_10_θ) * "\t" * string(min_loss) * "\n")
            end
        end
    end
end

using ..MultiMuC
MultiMuC.run(Int64, Float64)
