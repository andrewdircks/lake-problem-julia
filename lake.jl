using Roots
using Distributions
include("borg.jl")

# lake parameters
b = 0.42
q = 2.0

# natural phosphorus inflow parameters
mu = 0.03
sigma = (10 ^ -5) ^ 0.5

# economic parameters
alpha = 0.4
delta  =0.98

# thresholds
reliability_threshold = 0.85
inertia_threshold = -0.02

# decision variables, objectives and constraints 
nvars = 100
# nobjs = 4
nobjs = 2
nYears = 100
nSamples = 100
nSeeds = 1
# nconstrs = 1
nconstrs = 0

function pCrit(x)
    return (x^q) / (1+x^q) - b*x
end

critical_threshold = find_zero(pCrit, 0.5)

function lake(vars, objs, constrs)::Cvoid 
    # load vars
    vars = [unsafe_load(vars, i) for i = 1:nvars]

    # initialize arrays
    average_annual_P = zeros(Float64, nYears)
    discounted_benefit = zeros(Float64, nSamples)
    yrs_inertia_met = zeros(Float64, nSamples)
    yrs_Pcrit_met = zeros(Float64, nSamples)
    lake_state = zeros(Float64, nYears+1)
    _objs = zeros(Float64, nobjs)
    _constrs = zeros(Float64, nconstrs)

    # generate Monte Carlo sampoles of natural phosphorus inflows
    natFlow = zeros(Float64, nSamples, nYears)
    _mean =log((mu^2) / ((mu^2+sigma^2) ^ 0.5))
    _sigma = log((sigma^2+mu^2)/mu^2) ^ 0.5
    for i in 1:nSamples
        inflow_distribution = LogNormal(_mean, _sigma)
        natFlow[i,:] = rand(inflow_distribution, nYears)
    end

    # perform simulation
    for s in 1:nSamples
        for i in 1:nYears
            lake_state[i+1] = lake_state[i]*(1-b) + (lake_state[i]^q)/(1+(lake_state[i]^q)) + vars[i] + natFlow[s,i]
            average_annual_P[i] = average_annual_P[i] + lake_state[i+1]/nSamples
            discounted_benefit[s] = discounted_benefit[s] + alpha*vars[i]*(delta^i)
            
            if i>=2 && ((vars[i] - vars[i-1]) > inertia_threshold)
                yrs_inertia_met[s] = yrs_inertia_met[s] + 1
            end
            if lake_state[i+1] < critical_threshold
                yrs_Pcrit_met[s] = yrs_Pcrit_met[s] + 1
            end
        end
    end

    # objectives
    _objs[1] = -1*mean(discounted_benefit)                   #average economic benefit
    _objs[2] = maximum(average_annual_P)                     #minimize the max average annual P concentration
    # _objs[3] = -1*sum(yrs_inertia_met)/((nYears-1)*nSamples) #average percent of transitions meeting inertia thershold
    # _objs[4] = -1*sum(yrs_Pcrit_met)/(nYears*nSamples)       #average reliability

    # constraint
    # _constrs[1] = max(0.0, reliability_threshold - (-1*_objs[3]))

    # store objs
    for i = 1:nobjs
      unsafe_store!(objs, _objs[i], i)
    end

    # store constrs
    # for i = 1:nconstrs
    #   unsafe_store!(constrs, _constrs[i], i)
    # end
    return Cvoid()
end


bounds = repeat([(0.01, 0.1)], nvars)
epsilons = [0.01, 0.01, 0.0001, 0.0001]
test = borg(nvars, nobjs, 1, bounds, epsilons, 1000, lake)
settings = Dict("frequency" => 100, "runtimefile" => "test.txt")
result = run(test, settings)
