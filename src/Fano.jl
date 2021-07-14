module Fano


using LinearAlgebra
using DelimitedFiles
using Plots
using Optim

global temps = Float64.([10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 70; 80; 90])

include("Utility.jl")
include("Load.jl")

#load data
VH, VV = load()

include("Ham.jl")
include("MainFit.jl")


end #end module
