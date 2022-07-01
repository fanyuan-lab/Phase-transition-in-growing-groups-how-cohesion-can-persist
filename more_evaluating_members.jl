using Random, StatsBase, Statistics

Random.seed!(0)

function groupFormation(N::Int, η::Float64, f::Float64, m::Int)
    c = m #initial positive founders
    n = m #number of admitted nodes
    while n < N
        w = (1 - 2η) * c / n + η
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

function calc_cohesion(N::Int, η::Float64, f::Float64, m::Int; iterations=500)
    cohesion_vec = zeros(iterations)
    for iter in 1:iterations
        cohesion = groupFormation(N, η, f, m)
        cohesion_vec[iter] = cohesion
    end
    return mean(cohesion_vec), std(cohesion_vec)
end
