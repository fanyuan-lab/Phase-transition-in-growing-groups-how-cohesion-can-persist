using Random, Statistics, StatsBase

Random.seed!(0)

function groupFormation(N::Int, η::Float64, N0::Int, rule::String, f::Float64)
    R = zeros(Int, N, N) #relation matrix
    founders = [1:Int(N0);] #founder members
    types = zeros(Int, N) #types of all group members 
    types[founders] .= 1 #founder nodes are of positive types 
    if rule == "UC" #UC case
        for t in N0+1:N #group size N
            while true
                types[t] = sample([1,-1],pweights([f,1-f])) #candidate t is of type +1/-1 with probability f/(1-f)
                evaluator = rand(1:t-1) #one group member is chosen uniformly, i.e. UC case
                types[t] == types[evaluator] ? ((η < rand()) && (R[t, evaluator] = 1)) : ((η > rand()) && (R[t, evaluator] = 1))
                if R[t,evaluator] == 1
                    break 
                end
            end
        end
    elseif rule == "PA" #PA case
        counters = zeros(Int, N)
        N0 == 1 ? counters[1] = 1 : counters[founders] = repeat([N0 - 1], N0)
        for t in N0+1:N
            evaluator = 0
            while true
                types[t] = sample([1,-1],pweights([f,1-f])) #candidate t is of type +1/-1 with probability f/(1-f)
                evaluator = sample([1:t-1;], pweights(counters ./ sum(counters))) #preferential attachment, i.e. PA case
                types[t] == types[evaluator] ? ((η < rand()) ? (R[t, evaluator] = 1) : (R[t, evaluator] = -1)) : ((η > rand()) ? (R[t, evaluator] = 1) : (R[t, evaluator] = -1))
                R[t, evaluator] == -1 ? ((R[t, evaluator] = 0); continue) : break
            end
            (t == 2) && (counters[evaluator] -= 1)
            counters[evaluator] += 1
            counters[t] += 1
        end
    elseif rule == "DS" #Benchmark: Dictatorship 
        for t in N0+1:N 
            while true
                types[t] = sample([1,-1],pweights([f,1-f])) #candidate t is of type +1/-1 with probability f/(1-f)
                types[t] == 1 ? ((rand() < 1-η) ? break : nothing) : ((rand() < η) ? break : nothing)
            end
        end
    end
    return types
end


function calc_cohesion(N::Int, η::Float64, N0::Int, rule::String, f::Float64; iterations=10000)
    cohesion_vec = zeros(iterations)
    for iter in 1:iterations
        types = groupFormation(N, η, N0, rule, f)
        cohesion = sum(types .== 1) / N
        cohesion_vec[iter] = cohesion
    end
    return mean(cohesion_vec), std(cohesion_vec)
end




