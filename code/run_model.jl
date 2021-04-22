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


## RUN model  --------------------------

# Settings

params = Dict("beta" => .99, "eta" => 0.032, "alpha" => .115, "delta" => .05, 
	"g1" => .000177, "g2" =>.0044, "g31" => .000, "g32" => .039,  "g33" => .138);

T = 20;
Y = 10;
D = 10;

# Running and plotting

res1 = run_model(T, Y, D, .99, params, "Nordhaus");

res2 = run_model(T, Y, D, .99, params, "Medium");
plot_model(res2, res1, scat=false)

res3 = run_model(T, Y, D, .99, params, "Extreme");
plot_model(res3, res1, scat=false)

res = run_model(T, Y, D, .99, params, "none");
plot_model(res, res1, scat=false)
