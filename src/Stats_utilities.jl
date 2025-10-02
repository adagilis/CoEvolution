using Statistics
using Distributions
using FLoops
using DataFrames
using HypothesisTests
using Loess
using Missings
using CategoricalArrays

"""
    cor_test_matrix(m,dim) -> DataFrame(i,j,r,pval)
performs correlation tests for all rows or columns of a matrix, and returns a data frame listing the row vs column compared, the r value and the pvalue of a t-test
"""
function cor_test_matrix(m,dim;cutoff=nothing)
    dims=size(m)
    num_comp = binomial(dims[dim],2)
    cor_res = DataFrame(missings(Float64,num_comp,4),[:i,:j,:r,:pval])
    @floop ThreadedEx() for (i,j) in Iterators.product(1:(dims[dim]-1),1:dims[dim])
        if(j>i)
            index = index_func(i,j,dims[dim])
            if dim==1
                sx, sy = collect.(skipmissings(m[i,:],m[j,:]))
            elseif dim==2
                sx, sy = collect.(skipmissings(m[:,i],m[:,j])) 
            end
            if isfinite(cutoff)
                idx = intersect(findall(sx .< cutoff),findall(sy .< cutoff))
                sx = sx[idx]
                sy = sy[idx]
            end
            if length(sx)>4

                cor_test = CorrelationTest(sx,sy)
                cor_res[index,:] = Dict(:i=>i,:j=>j,:r=>cor_test.r,:pval=>pvalue(cor_test))
            else 
                cor_res[index,:] = Dict(:i=>i,:j=>j,:r=>missing,:pval=>missing)
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

"""
    nonanmissing(x)
quick convenience function, similaar to skipmissing, but also skips nans
"""
function nonanmissing(x)
    idnan = findall((!isnan).(x))
    idmiss = findall((!ismissing).(x[idnan]))
    return(x[idnan[idmiss]])
end

"""
    nonanmissing(x,y)
quick convenience function, similar to skipmissings, but here checking for missings in either of two vectors. 
"""
function nonanmissing(x,y)
    idnan = intersect(findall((!isnan).(x)),findall((!isnan).(y)))
    idmiss = intersect(findall((!ismissing).(x[idnan])),findall((!ismissing).(y[idnan])))
    return([x[idnan[idmiss]],y[idnan[idmiss]]])
end

"""
    range_dist(x,y)
Returns the distance between two ranges, returning 0 if they overlap.
"""
function range_dist(x,y)
    if sort(x)[1] < sort(y)[1]
        return(maximum([0,sort(y)[1]-sort(x)[2]]))
    else
        return(maximum([0,sort(x)[1]-sort(y)[2]]))
    end
end

"""
    binned_mean(df::DataFrame,x_wnd::Int,y_wnd::Int)
creates a 2d set of bins, with specified window sizes
Assumed that the data frames first two columns are the x and y values.
Currently only returns the mean, but should be re-writable to return any function desired.
"""
function binned_mean(df::DataFrame,x_wnd::Int,y_wnd::Int)
    x_edge = 0:x_wnd:(maximum(df[:,2])+1)
    y_edge = 0:y_wnd:(maximum(df[:,2])+1)
    fmt(from, to, i; leftclosed, rightclosed) = (from + to)/2
    df.x_bin = cut(df[:,1], x_edge,labels=fmt;extend=true)
    df.y_bin = cut(df[:,2], y_edge,labels=fmt;extend=true)
    df.x_bin = unwrap.(df.x_bin)
    df.y_bin = unwrap.(df.y_bin)
    return(combine(groupby(df,[:x_bin,:y_bin]),nrow, 3=>mean))
end

"""
    binned_mean(x,y,x_wnd)
A one dimensional binning approach, where x indicates the values to be binned along, while y are the values to return the mean of.
e.g. chromosomal position versus some statistic of interest. Currently only supports returning the mean, and does not support sliding windows.
"""
function binned_mean(x::Vector,y::Vector,x_wnd::Int)
    x_edge = 0:x_wnd:(maximum(df[:,2])+1)
    fmt(from, to, i; leftclosed, rightclosed) = (from + to)/2
    x_bin = cut(x, x_edge,labels=fmt;extend=true)
    x_bin = unwrap.(x_bin)
    df = DataFrame(x_bin=x_bin,y=y)
    return(combine(groupby(df,:x_bin),nrow, 2=>mean))
end
