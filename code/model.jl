#=
macro_model_climate:
- Julia version: 
- Author: mreguant
- Date: 2021-04-20
=#

using CSV, DataFrames, DataFramesMeta
using JuMP, PiecewiseLinearOpt # used for mathematical programming
using Ipopt # solver
using Cbc, Gurobi
using Printf
using LinearAlgebra
using Statistics


if ENV["USER"]=="mreguant"
    path = "/Users/mreguant/Documents/git/macro-annual-comment/"
end


# this function takes correlations in uncertainty and generates transition paths
function uncertainty_model(df_clim::DataFrame, T::Int64, Y::Int64, D::Int64; damage="none", future="random")

	df = DataFrame(p = Int64[], s = Int64[], s_1 = Int64[], ind_tr = Array{Int64,1}[], tr = Array{Float64,1}[], 
		zeta = Float64[], g3 = Float64[], Y_last = Bool[]);

    # year 0
    g = [0.0, 0.032, 0.138];
    if (damage=="Nordhaus")
        g = [g[1]]
    end
    if (damage=="Medium")
        g = [g[2]]
    end
    if (damage=="Extreme")
        g = [g[3]]
    end

    # start year
    push!(df, [1, 1, 0, [Int(i) for i=2:1+D], [1.0/D for i=1:D], Y*mean(df_clim.Column1), 0.0, false])

    # populating the future
    if (future=="random")
        ct = 1;
        for p = 2:T
            Dp = D^(p-1);   # draws accumulated this period
            Dp_1 = D^(p-2); # draws previous period
            for d = 1:D
                push!(df, [p, ct+d, max(1+D*(p-3)+d,1), [Int(i) for i=2+D*(p-1):1+D*p], [1.0/D for i=1:D], Y*sample(df_clim.Column1), sample(g), (p==T)])
                #push!(df, [p, 1+(p-2)*D+d, max(1+D*(p-3)+d,1), [Int(1+D*(p-1)+d)], [1.0], Y*sample(df_clim.Column1), sample(g), (p==T)])
            end
            ct = ct + D;
        end
    end

	return df;

end

# this function solves for the dynamic game at T0 (year 2020)
function solve_model(pm::Dict{String,Float64}, df::DataFrame)

    model = Model(
        optimizer_with_attributes(Ipopt.Optimizer,  "print_level" => 1,
            "acceptable_tol" => 1e-12)
        );

    # model = Model(
    #     optimizer_with_attributes(NLopt.Optimizer,  "algorithm" => :LN_COBYLA)
    #     );


    S = Int(nrow(df));
    S_1 = Int(maximum(df.s .* (1 .-df.Y_last)));
    S_E = Int(maximum(df.s .* (df.p.<=10)));

    # Convert d

    # main variables to solve for
    @variable(model, K[1:S]);  # capital
    @variable(model, C[1:S] >= 0.0);  # consumption
    @variable(model, I[1:S] >= 0.0);  # investment
    @variable(model, N[1:S] >= 1.0);  # damages
    @variable(model, 0.0 <= E[1:S] <= 60.0);  # energy (allow at most 60 Gt per year)
    @variable(model, Y[1:S] >= 1.0);  # temperature
    @variable(model, U[1:S]);  # utility function
    @variable(model, V[1:S]);  # value function


    # objective function
    @objective(model, Max, V[1]);

    # initialization constraints
    @constraint(model, K[1]==1000.0);
    @constraint(model, Y[1]==1.2);

    # no borrowing at capital T
    @constraint(model, [s in S_1+1:S], I[s] >= 0.0);

    # definition utility function
    @NLconstraint(model, [s=1:S], U[s] == (C[s]/N[s])^(1-pm["eta"]) * E[s]^pm["eta"]);

    # value function - Bellman
    @constraint(model, [s in S_1+1:S], V[s] == U[s]);
    @constraint(model, [s in 1:S_1], V[s] == U[s] + pm["beta"] * sum(df.tr[s][i] * V[df.ind_tr[s][i]] for i in 1:length(df.tr[s])));

    # law of motion and definition constraints
    @constraint(model, [s=2:S], K[s] == I[df.s_1[s]] + (1.0-pm["delta"]) * K[df.s_1[s]]);
    @constraint(model, [s=1:S], C[s] + I[s] == pm["alpha"] * K[s]);
    @constraint(model, [s=2:S], Y[s] == Y[df.s_1[s]] + df.zeta[s] * E[df.s_1[s]]/3700.0);  # temperature equation
    @NLconstraint(model, [s=1:S], log(N[s]) == pm["g1"] * Y[s] + pm["g2"] * Y[s]^2 + df.g3[s] * (Y[s]-Y[1])^3);  # damage equation

    # constraint on last periods
    @constraint(model, [s=S_1+1:S], E[s] <= sum(E[spast] for spast in 1:S_1)/S_1);

    optimize!(model);

    f = open(string(path,"code/model.lp"), "w")
    print(f, model)
    close(f)

    status = @sprintf("%s", JuMP.termination_status(model));
    if (status=="LOCALLY_SOLVED") | (status=="ALMOST_LOCALLY_SOLVED")

        results = DataFrame(p = Int64[], s = Int64[], U = Float64[], V = Float64[], K = Float64[],
            I = Float64[], C = Float64[], E = Float64[], Y = Float64[], N = Float64[]);
        for s=1:S
            push!(results, [df.p[s], s, JuMP.value(U[s]), JuMP.value(V[s]), JuMP.value(K[s]),
                JuMP.value(I[s]), JuMP.value(C[s]), JuMP.value(E[s]), JuMP.value(Y[s]), JuMP.value(N[s])])
        end
        return results;

    else

        return status;

    end

end


function run_model(T::Int64, Y::Int64, D::Int64, beta::Float64,  pm::Dict{String,Float64}, dmg::String)
    
    params["beta"] = 1/((2.0-beta)^Y);

    # Reading in distribution of climate sensitivities and computing uncertainty
    df_clim = CSV.read(string(path,"input/model144.csv"), header=false, DataFrame);

    # run several draws
    results = DataFrame(p = Int64[], s = Int64[], U = Float64[], V = Float64[], K = Float64[],
            I = Float64[], C = Float64[], E = Float64[], Y = Float64[], N = Float64[]);
    for d in 1:50
        df = uncertainty_model(df_clim, T, Y, D, damage=dmg);
        res = solve_model(params, df);
        append!(results, res);
    end

    # keep only up to 2100
    results.p = (results.p .- 1.0)*Y .+ 2020;
    @linq df_plt = results |> where(results.p .<= 2100);
    df_plt = combine(groupby(df_plt,[:p, :s]),:U=>mean=>:U, :V=>mean=>:V, :K=>mean=>:K,
        :I=>mean=>:I, :C=>mean=>:C, :E=>mean=>:E, :Y=>mean=>:Y, :N=>mean=>:N);

    return df_plt;

end


function plot_model(df_plt::DataFrame, ref_case::DataFrame; scat=true)
    
    if (scat==true)
        p1 = scatter(df_plt.p, df_plt.E/3.7, xlabel = "Emissions", ylabel = "Gt C") # Make a line plot
        p2 = scatter(df_plt.p, df_plt.Y, xlabel = "Temperature", ylabel = "degree C") # Make a scatter plot
        p3 = scatter(df_plt.p, 1.0 ./df_plt.N, xlabel = "1/N")
        p4 = scatter(df_plt.p, df_plt.I, xlabel = "Investment")
    else
        df_E = combine(groupby(df_plt,:p), :E=>mean, :Y=>mean, :N=>mean, :I=>mean, :C=>mean, :U=>mean);
        df_ref = combine(groupby(ref_case,:p), :E=>mean, :Y=>mean, :C=>mean, :U=>mean);
        p1 = plot(df_E.p, df_E.E_mean/3.7, xlabel = "Emissions", ylabel = "Gt C") # Make a line plot
        p1 = plot!(df_ref.p, df_ref.E_mean/3.7, xlabel = "Emissions", ylabel = "Gt C") # Make a line plot
        p2 = plot(df_E.p, df_E.Y_mean, xlabel = "Temperature", ylabel = "degree C") # Make a scatter plot
        p2 = plot!(df_ref.p, df_ref.Y_mean, xlabel = "Temperature", ylabel = "degree C") # Make a scatter plot
        p3 = plot(df_E.p, 1.0 ./df_E.N_mean, xlabel = "1/N")
        p4 = plot(df_E.p, df_E.U_mean./df_ref.U_mean, xlabel = "Rel. U to Case 1")
    end
        
    plot(p1, p2, p3, p4, layout = (2, 2), legend = false)

end

