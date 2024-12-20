using DrWatson, Test
@quickactivate "CoEvolution"

# Here you include files using `srcdir`
include(srcdir("ERC_functions.jl"))

# Run test suite
println("Starting tests for ERC_simulations.jl")
ti = time()
trees = filter(contains("test_tree"),readdir(projectdir("test","data"),join=true))
species_tree = open(parsenewick,projectdir("test","data","species_tree.newick"))
t1 = read_tree(trees[1])
t2 = read_tree(trees[2])
t3 = read_tree(trees[3])
t4 = read_tree(trees[4])
@testset "ERC functions tests" begin
    @test 1 == 1
    @test calculate_ERC(1,2,t1,t2,species_tree,5)[:r2] == 1
    @test calculate_ERC(1,3,t1,t3,species_tree,5)[:r2] == 1
    @test calculate_ERC(1,3,t1,t3,species_tree,5)[:r2] == 1 
    @test calculate_ERC(3,4,t3,t4,species_tree,5)[:r2] == 1 
    @test runERC_files(trees)[:,:r2] == [1,1,1,1,1,1]
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

using UnicodePlots
include(srcdir("ERC_simulations.jl"))
println("Testing ERC_simulations.jl")
ti = time()
@testset "ERC_simulations tests" begin
    trees = run_simulations(species_tree,10000,500)
    ERC = runERC_collection(trees)
    histogram(ERC[:,:r2])
end
ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60),digits=3," minutes to generate 1000 trees and calculate ERC")