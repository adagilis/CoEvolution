using DrWatson, Test
@quickactivate "CoEvolution"

# Here you include files using `srcdir`
include(srcdir("ERC_functions.jl"))

# Run test suite
println("Starting tests")
ti = time()
trees = filter(contains("test_tree"),readdir(projectdir("test","data"),join=true))
species_tree = open(parsenewick,projectdir("test","data","species_tree.newick"))
@testset "CoEvolution tests" begin
    @test 1 == 1
    @test calculate_ERC(1,2,species_tree,5)[:r2] == 1
    @test calculate_ERC(1,3,species_tree,5)[:r2] == 1
    @test calculate_ERC(2,3,species_tree,5)[:r2] == 1 #TODO: fix issues when branch order not identical
    @test calculate_ERC(3,4,species_tree,5)[:r2] == 1 #TODO: fix issues when branch order not identical
    @test runERC()[:,:r2] == [1,1,1,1,1,1]
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
