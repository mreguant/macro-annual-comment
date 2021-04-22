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
using Statistics, StatsBase
using Plots


if ENV["USER"]=="mreguant"
    path = "/Users/mreguant/Documents/git/macro-annual-comment/"
end
include(string(path,"code/model.jl"));



## RUN Test Model  -----------------------------


# Setting model parameters
params_simple = Dict("beta" => .97, "eta" => 0.032, "alpha" => .115, "delta" => .01, 
	"g1" => .000177, "g2" =>.0044, "g31" => .000, "g32" => .039,  "g33" => .138);

# Reading test data frame
df_simple = CSV.read(string(path,"input/simple_test.csv"), DataFrame);
df_simple[!,:tr] = repeat([[.4,.3,.3]],nrow(df_simple));
df_simple[!,:ind_tr] = repeat([[2,3,4]],nrow(df_simple));

# Solve
simple_test = solve_model(params_simple, df_simple);



## RUN actual model  --------------------------

# Reading in distribution of climate sensitivities and computing uncertainty
df_clim = CSV.read(string(path,"input/model144.csv"), header=false, DataFrame);

function run_test(T::Int64, Y::Int64, D::Int64, beta::Float64,  pm::Dict{String,Float64}, dmg::String)
	
	params["beta"] = 1/((2.0-beta)^Y);
	df = uncertainty_model(df_clim, T, T1, Y, D, damage=dmg);
	results = solve_model(pm, df);
	results.p = (results.p .- 1.0)*Y .+ 2020;
	
	@linq df_plt = results |> where(results.p .<= 2100);

	return df_plt;

end

function plot_test(df_plt::DataFrame)
	df_E = combine(groupby(results,:p),:E=>mean);
	p1 = plot(df_E.p, df_E.E_mean/3.7, xlabel = "Emissions", ylabel = "Gt C") # Make a line plot
	p2 = scatter(df_plt.p, df_plt.Y, xlabel = "Temperature", ylabel = "degree C") # Make a scatter plot
	p3 = scatter(df_plt.p, 1.0 ./df_plt.N, xlabel = "1/N")
	p4 = scatter(df_plt.p, df_plt.I, xlabel = "Investment")
	plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
end

include(string(path,"code/model.jl"));
params = Dict("beta" => .99, "eta" => 0.032, "alpha" => .115, "delta" => .05, 
	"g1" => .000177, "g2" =>.0044, "g31" => .000, "g32" => .039,  "g33" => .771);

results = run_test(20, 8, 6, .99, params, "Nordhaus");
plot_test(results)

results = run_test(20, 8, 6, .99, params, "Extreme");
plot_test(results)

results = run_test(20, 8, 6, .99, params, "Medium");
plot_test(results)

results = run_test(20, 8, 6, .99, params, "none");
plot_test(results)

