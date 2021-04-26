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
	"g1" => .000177, "g2" =>.0044, "g31" => .000, "g32" => .039,  "g33" => .06,
	"Emin" => 0.0, "Emax" => 50.0);

# Uncertainty modeling
T = 18;
T1 = 3;
Y = 10;
D = 3;
nummc = 10;
beta = .99;

# Running and plotting

for c in ["average"]
	for cu in ["random", "persistent"]

		res1 = run_model(T, T1, Y, D, beta, params, dmg="Nordhaus", clim=c, climunc=cu, mc=nummc);

		res2 = run_model(T, T1, Y, D, beta, params, dmg="Medium", clim=c, climunc=cu, mc=nummc);
		plot_model(res2, res1, scat=false)
		Plots.savefig(string(path,"output/outcome_Medium_",c,"_",cu,".pdf"));
		plot_model(res2, res1, scat=true)
		Plots.savefig(string(path,"output/scatter_Medium_",c,"_",cu,".pdf"));

		res3 = run_model(T, T1, Y, D, beta, params, dmg="Extreme", clim=c, climunc=cu, mc=nummc);
		plot_model(res3, res1, scat=false)
		Plots.savefig(string(path,"output/outcome_Extreme_",c,"_",cu,".pdf"));
		plot_model(res3, res1, scat=true)
		Plots.savefig(string(path,"output/scatter_Extreme_",c,"_",cu,".pdf"));

		plot_model(res3, res2, scat=false)
		Plots.savefig(string(path,"output/outcome_comparison_",c,"_",cu,".pdf"));

		for f in ["unknown", "gradual"]	
			res = run_model(T, T1, Y, D, beta, params, dmgunc=f, clim=c, climunc=cu, mc=nummc);
			plot_model(res, res1, scat=false)
			Plots.savefig(string(path,"output/outcome_",c,"_",cu,"_",f,".pdf"));
			plot_model(res, res1, scat=true)
			Plots.savefig(string(path,"output/scatter_",c,"_",cu,"_",f,".pdf"));
		end

	end
end
