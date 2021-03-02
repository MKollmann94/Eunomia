__precompile__()

module Eunomia

##########################################################################################

#### 	REQUIREMENTS

##########################################################################################

using Eirene
using DelimitedFiles
using CSV
using Suppressor
using SparseArrays

##########################################################################################

#### 	EXPORTS

##########################################################################################

export eunomia

##########################################################################################

#### 	GLOBAL CONSTANTS

##########################################################################################

const FILEPATH = @__DIR__


function eunomia(M)
    n, m = size(M)
    ev = numOfSimplices(M)
    #fv = ones(sum(ev))
    fv, simp = dowkerActivationTimes(M)

    #it is important to load these arrays as Int64, otherwise Eirene will not work properly
    rv = readdlm(FILEPATH * "\\distanceFiles\\rv_$m.txt", Int64)[:,1]
    cp = readdlm(FILEPATH * "\\distanceFiles\\cp_$m.txt", Int64)[:,1]

    C = eirene(model = "complex", rv = rv, cp = cp, ev = ev, fv = fv)

    return C
end


#=
    returns the number of simplices for each dimension (up to 3)
=#
function numOfSimplices(M)
    n, m = size(M)
    ev = []

    for k in 1 : 4
        push!(ev, binomial(m, k))
    end
    push!(ev, 0)

    return ev
end

#=
    calculate the activation time of the simplices and
    add them to fv
    it is necessary to use different for-loops
    since Eirene needs the simplices ordered by dimension
=#
function dowkerActivationTimes(M)
    n, m = size(M)
    fv = []
    simp = []

    for i in 1:m
        push!(fv, minimum(M[:,i]))
        push!(simp, (i))
    end

    for i in 1:m
        for j in i+1:m
            push!(fv, minimum(max(M[col, i], M[col, j]) for col in 1:n))
            push!(simp, (i, j))
        end
    end

    for i in 1:m
        for j in i+1:m
            for k in j+1:m
                push!(fv, minimum(max(M[col, i], M[col, j], M[col, k]) for col in 1:n))
                push!(simp, (i, j, k))
            end
        end
    end

    for i in 1:m
        for j in i+1:m
            for k in j+1:m
                for l in k+1:m
                    push!(fv, minimum(max(M[col, i], M[col, j], M[col, k], M[col, l]) for col in 1:n))
                    push!(simp, (i, j, k, l))
                end
            end
        end
    end

    return fv, simp
end

#=
    calculates the distance matrix for a matrix with m columns
    and stores the vectors needed for Eirene in the given filepath
=#
function distanceD(m)
    println(m)
    D = zeros(Int8, sumnchoosek(m, 4), sumnchoosek(m, 4))
    simp = []

    for i in 1:m
        push!(simp, (i))
    end
    for i in 1:m
        for j in i+1:m
            push!(simp, (i,j))
        end
    end

    for i in 1:m
        for j in i+1:m
            for k in j+1:m
                push!(simp, (i,j,k))
            end
        end
    end

    for i in 1:m
        for j in i+1:m
            for k in j+1:m
                for l in k+1:m
                    push!(simp, (i,j,k,l))
                end
            end
        end
    end

    #D[i,j] = 1 iff i is a codimension-1 face of cell j
    for k in 1:3
        for idx_lower in sumnchoosek(m, k - 1) + 1 : sumnchoosek(m, k)
            for idx_higher in sumnchoosek(m, k) + 1 : sumnchoosek(m, k + 1)
                contained = true
                for element in simp[idx_lower]
                    contained = contained && element in simp[idx_higher]
                end
                if contained
                    D[idx_lower, idx_higher] = 1
                end
            end
        end
    end

    S = sparse(D)
    rv = S.rowval
    cp = S.colptr

    writedlm(FILEPATH * "\\distanceFiles\\rv_$m.txt", rv)
    writedlm(FILEPATH * "\\distanceFiles\\cp_$m.txt", cp)
end


function sumnchoosek(m,k)
    sum = 0
    for i in 1:k
        sum += binomial(m, i)
    end
    return sum
end


function saveDistances(dfInput, k, folderUrl)
    println("Start with k = $k")
    M = Array(dfInput)
    n, m = size(M)
    distancesDim0 = []
    distancesDim1 = []
    distancesDim2 = []
    C = callEirene(M[:, 1 : k])
    for i in 2 : (m - k) + 1
        if i % 1000 == 0
            println(i)
        end
        D = callEirene(M[:, i : i + k - 1])
        c0 = barcode(C, dim = 0)
        d0 = barcode(D, dim = 0)
        c1 = barcode(C, dim = 1)
        d1 = barcode(D, dim = 1)
        c2 = barcode(C, dim = 2)
        d2 = barcode(D, dim = 2)
        @suppress begin
            push!(distancesDim0, wasserstein_distance(c0, d0, p = 1))
            push!(distancesDim1, wasserstein_distance(c1, d1, p = 1))
            push!(distancesDim2, wasserstein_distance(c2, d2, p = 1))
        end
        C = D
    end
    distances = Array{Any, 2}(undef, size(distancesDim0, 1) + 1, 4)
    distances[1, 1] = "Dates"
    distances[1, 2] = "Distances Dim 0"
    distances[1, 3] = "Distances Dim 1"
    distances[1, 4] = "Distances Dim 2"
    distances[2 : end, 1] = names(dfInput)[1 : size(distancesDim0, 1)]
    distances[2 : end, 2] = distancesDim0
    distances[2 : end, 3] = distancesDim1
    distances[2 : end, 4] = distancesDim2
    writedlm("$folderUrl\\distances_$k.csv", distances, ";")
end

end
