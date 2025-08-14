using DrWatson, Test
@quickactivate "CoEvolution"

# Here you include files using `srcdir`
include(srcdir("Phylo_utilities.jl"))
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
    @test calculate_ERC(1,2,t1,t2,species_tree;cutoff=5)[:r] == 1
    @test calculate_ERC(1,3,t1,t3,species_tree;cutoff=5)[:r] == 1
    @test calculate_ERC(1,3,t2,t3,species_tree;cutoff=5)[:r] == 1 
    @test calculate_ERC(3,4,t3,t4,species_tree;cutoff=5)[:r] == 1 
    @test runERC_files(trees,species_tree)[:,:r] == [1,1,1,1,1,1]
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

include(srcdir("ERC_simulations.jl"))
println("Testing ERC_simulations.jl")
ti = time()
trees_simulated = simulate_tree_set(species_tree,0.01,100,length(getleaves(species_tree)),length(getleaves(species_tree)))
ERC_simulation = runERC_collection(trees_simulated,species_tree)

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60,digits=3)," minutes to generate 100 trees and calculate ERC")

