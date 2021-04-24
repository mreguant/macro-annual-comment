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


## RUN model  --------------------------

# include file
include(string(path,"code/model.jl"));

# Settings
params = Dict("beta" => .99, "eta" => 0.032, "alpha" => .115, "delta" => .05, 
	"g1" => .000177, "g2" =>.0044, "g31" => .000, "g32" => .039,  "g33" => .138);

# Uncertainty modeling
T = 18;
T1 = 3;
Y = 10;
D = 3;

# Running and plotting

for c in ["none", "average"]
	for cu in ["none", "persistent"]

		res1 = run_model(T, T1, Y, D, .99, params, dmg="Nordhaus", clim=c, climunc=cu);

		res2 = run_model(T, T1, Y, D, .99, params, dmg="Medium", clim=c, climunc=cu);
		plot_model(res2, res1, scat=false)
		Plots.savefig(string(path,"output/outcome_Medium_",c,"_",cu,".pdf"));
		plot_model(res2, res1, scat=true)
		Plots.savefig(string(path,"output/scatter_Medium_",c,"_",cu,".pdf"));

		res3 = run_model(T, T1, Y, D, .99, params, dmg="Extreme", clim=c, climunc=cu);
		plot_model(res3, res1, scat=false)
		Plots.savefig(string(path,"output/outcome_Extreme_",c,"_",cu,".pdf"));
		plot_model(res3, res1, scat=true)
		Plots.savefig(string(path,"output/scatter_Extreme_",c,"_",cu,".pdf"));

		for f in ["unknown", "gradual"]	
			res = run_model(T, T1, Y, D, .99, params, dmgunc=f, clim=c, climunc=cu);
			plot_model(res, res1, scat=false)
			Plots.savefig(string(path,"output/outcome_",f,"_",c,"_",cu,".pdf"));
			plot_model(res, res1, scat=true)
			Plots.savefig(string(path,"output/scatter_",f,"_",c,"_",cu,".pdf"));
		end

	end
end
