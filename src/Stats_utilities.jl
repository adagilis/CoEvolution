using Statistics
using Distributions
using FLoops
using DataFrames
using HypothesisTests
using Loess

"""
    cor_test_matrix(m,dim) -> DataFrame(i,j,r,pval)
    performs correlation tests for all rows or columns of a matrix, and returns a data frame listing the row vs column compared, the r value and the pvalue of a t-test
"""

function cor_test_matrix(m,dim)
    dims=size(m)
    num_comp = binomial(dims[dim],2)
    cor_res = DataFrame(zeros(num_comp,4),[:i,:j,:r,:pval])
    @floop ThreadedEx() for (i,j) in Iterators.product(1:(dims[dim]-1),1:dims[dim])
        if(j>i)
            index = index_func(i,j,dims[dim])
            if dim==1
                cor_test = CorrelationTest(m[i,:],m[j,:])
                cor_res[index,:] = Dict(:i=>i,:j=>j,:r=>cor_test.r,:pval=>pvalue(cor_test))
            elseif dim==2
                cor_test = CorrelationTest(m[:,i],m[:,j])
                cor_res[index,:] = Dict(:i=>i,:j=>j,:r=>cor_test.r,:pval=>pvalue(cor_test))
            end
        end
    end
    return(cor_res)
end

"""
    index_func(i,j,n) -> Int
    returns the index of the pair i and j among n choose 2 elements.
"""
function index_func(i,j,n)
    k = 0
    if(i>1)
        k = (i-1)*n-sum(1:i .- 1)
    end
    return(k+j-i)
end



"""
    zscores([x]) -> [z]
    returns the z-scores of a vector
"""
function zscores(x)
    return((x .- mean(x))./ std(x))
end


"""
    loess!(plot,x,y)
    adds a loess curve to a plot
"""

function loess!(p,x,y)
    model = loess(x,y)
    us = range(extrema(x)...;step=(extrema(x)[2]-extrema(x)[1])/100)
    vs = predict(model,us)
    plot!(p,us,vs,legend=false)
end