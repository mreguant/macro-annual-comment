#=
macro_model_climate:
- Julia version: 
- Author: mreguant
- Date: 2021-04-20
=#

using CSV, DataFrames, DataFramesMeta
using JuMP, PiecewiseLinearOpt # used for mathematical programming
using Ipopt # solver
using Printf
using LinearAlgebra
using Statistics, StatsBase
using Plots


if ENV["USER"]=="mreguant"
    path = "/Users/mreguant/Documents/git/macro-annual-comment/"
end


## RUN model  --------------------------

# include file
include(string(path,"code/model_cake.jl"));

# Settings
params = Dict("beta" => .99, "eta" => 0.032, 
	"g1" => .000177, "g2" =>.0044, "g31" => .000, "g32" => .016,  "g33" => .09,  
	"pr1" => 1.0/3.0, "pr2" => 1.0/3.0,  "pr3" => 1.0/3.0,
	"Y0" => 1.2, "Emin" => 1.0, "Emax" => 45.0, "E0" => 0.0, "Leontieff" => 0.0);

# Uncertainty modeling
T = 20;
T1 = 3;
Y = 10;
D = 3;  # needs to be a multiple of 3
nummc = 20;
beta = .99;
c = "average";
cu = "random";
f  = "unknown";

# Running and plotting baseline compared to BAU
res1 = run_model(T, T1, Y, D, beta, params, dmg="Nordhaus", clim=c, climunc=cu, mc=nummc);
resbase = run_model(T, T1, Y, D, beta, params, dmgunc=f, clim=c, climunc=cu, mc=nummc);
plot_model(resbase, res1, only_emissions=true, name=string(path,"output/outcome_baseline"));
plot_model(resbase, res1, scat=true)
Plots.savefig(string(path,"output/scatter_baseline.pdf"));

# adding ambiguity beliefs (more conservative)
params["pr1"] = 0.1;
params["pr2"] = 0.35;
params["pr3"] = 0.55;
res = run_model(T, T1, Y, D, beta, params, dmgunc=f, clim=c, climunc=cu, mc=nummc);
plot_model(res, resbase, only_emissions=true, name=string(path,"output/outcome_ambiguity"));
plot_model(res, resbase, scat=true)
Plots.savefig(string(path,"output/scatter_ambiguity.pdf"));

# adding gradual learning about scenarios
params["pr1"] = 1.0/3.0;
params["pr2"] = 1.0/3.0;
params["pr3"] = 1.0/3.0;
res = run_model(T, T1, Y, D, beta, params, dmgunc="gradual", dmgalpha=0.1, clim=c, climunc=cu, mc=nummc);
plot_model(res, resbase, only_emissions=true, name=string(path,"output/outcome_gradual_learning"));
plot_model(res, resbase, scat=true)
Plots.savefig(string(path,"output/scatter_gradual_learning.pdf"));

# adding emissions constraints
params["Emin"] = 15.0;
res2 = run_model(T, T1, Y, D, beta, params, dmg="Medium", clim=c, climunc=cu, mc=nummc);
res3 = run_model(T, T1, Y, D, beta, params, dmg="Extreme", clim=c, climunc=cu, mc=nummc);
plot_model(res3, res2, only_emissions=true, name=string(path,"output/outcome_comparison_Emin"));

# adding gradual constraints
params["Emin"] = 1.0;
params["E0"] = 45.0;
res2 = run_model(T, T1, Y, D, beta, params, dmg="Medium", clim=c, climunc=cu, mc=nummc);
res3 = run_model(T, T1, Y, D, beta, params, dmg="Extreme", clim=c, climunc=cu, mc=nummc);
plot_model(res3, res2, only_emissions=true, name=string(path,"output/outcome_comparison_Decrease"));

# adding gradual constraints + min 
params["Emin"] = 15.0;
params["E0"] = 45.0;
res2 = run_model(T, T1, Y, D, beta, params, dmg="Medium", clim=c, climunc=cu, mc=nummc);
res3 = run_model(T, T1, Y, D, beta, params, dmg="Extreme", clim=c, climunc=cu, mc=nummc);
plot_model(res3, res2, only_emissions=true, name=string(path,"output/outcome_comparison_Constraints"));

# adding gradual reductions + Leontieff constraints
params["Emin"] = 1.0;
params["Leontieff"] = 0.02;
res2 = run_model(T, T1, Y, D, beta, params, dmg="Medium", clim=c, climunc=cu, mc=nummc);
res3 = run_model(T, T1, Y, D, beta, params, dmg="Extreme", clim=c, climunc=cu, mc=nummc);
plot_model(res3, res2, only_emissions=true, name=string(path,"output/outcome_comparison_Leontieff"));
