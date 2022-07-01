using Random, StatsBase, Statistics

Random.seed!(0)

function groupFormation(N::Int, η::Float64)
    c = 2 #the number of founders
    n = 2 #the number of admitted members
    m = c #m evaluating members at each time step
    while n < N
        w = (1 - 2η) * c / n + η
        f = 1 / 2 #f∈[0,1]
        tp = f * w^m / (f * w^m + (1 - f) * (1 - w)^m)
        if tp > rand()
            c += 1
            n += 1
        else
            n += 1
        end
    end
    return c / N
end

function measure(N::Int, η::Float64; iterations=500)
    cohesion_vec = zeros(iterations)
    for iter in 1:iterations
        cohesion = groupFormation(N, η)
        cohesion_vec[iter] = cohesion
    end
    return mean(cohesion_vec), std(cohesion_vec)
end













