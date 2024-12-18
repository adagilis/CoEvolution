using DrWatson, Test
@quickactivate "CoEvolution"

# Here you include files using `srcdir`
include(srcdir("ERC_functions.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "CoEvolution tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
