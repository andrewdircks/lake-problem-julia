# https://github.com/andrewdircks/borg-julia

using Libdl, Dates, Printf

# TODO: make this variable
const temp = "borg/libborg.so"

mutable struct Configuration
    so_dir::String
    lib::Ptr{Nothing}
end

function Configure()::Configuration
    so_dir = "borg/libborg.so"
    lib = dlopen(so_dir)
    return Configuration(so_dir, lib)
end

# configure and set finalizer
config = Configure()
finalizer(config) do x
    dlclose(x.lib)
end

struct Borg
    nVars::Int32
    nObjs::Int32
    nCons::Int32
    bounds::Vector{Tuple{Float64,Float64}}
    epsilons::Vector{Float64}
    nfe::Int32
    problem::Base.CFunction
end

struct Solution 
    vars::Vector{Float64}
    objs::Vector{Float64}
    constrs::Vector{Float64}
end

function borg(nVars, nObjs, nCons, bounds, epsilons, nfe, julia_problem)::Borg
    c_problem = @cfunction($julia_problem, Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}))
    return Borg(nVars, nObjs, nCons, bounds, epsilons, nfe, c_problem)
end


# set the bounds to this BORG_Problem reference
function setBounds(bounds, nvars, ref)
    for i = 0:(nvars - 1)
        ccall((:BORG_Problem_set_bounds, temp), Cvoid, (Ptr{Cvoid}, Int32, Cdouble, Cdouble), ref, i, bounds[i + 1][1], bounds[i + 1][2])
    end
end


# set the epsilons to this BORG_Problem reference
function setEpsilons(epsilons, nobjs, ref)
    for i = 0:nobjs - 1
        ccall((:BORG_Problem_set_epsilon, temp), Cvoid, (Ptr{Cvoid}, Int32, Cdouble), ref, i, epsilons[i + 1])
    end
end


function getObjectives(sol, nObjs)
    return [ccall((:BORG_Solution_get_objective), Cdouble, (Ptr{Cvoid}, Int32), sol, i) for i = 0:(nObjs - 1)]
end

function getVariables(sol, nVars)
    return [ccall((:BORG_Solution_get_variable), Cdouble, (Ptr{Cvoid}, Int32), sol, i) for i = 0:(nVars - 1)]
end

function getConstraints(sol, nCons)
    return [ccall((:BORG_Solution_get_constraint), Cdouble, (Ptr{Cvoid}, Int32), sol, i) for i = 0:(nCons - 1)]
end


# process the result (pointer to archive)
function process(result, borg) 
    archive_size = ccall((:BORG_Archive_get_size, temp), Int32, (Ptr{Cvoid},), result)
    solutions = Vector{Solution}(undef, archive_size)

    for i = 0:(archive_size - 1)
        sol = ccall((:BORG_Archive_get, temp), Ptr{Cvoid}, (Ptr{Cvoid}, Int32), result, i)
        vars = getVariables(sol, borg.nVars)
        objs = getObjectives(sol, borg.nObjs)
        constrs = getConstraints(sol, borg.nCons)
        solutions[i + 1] = Solution(vars, objs, constrs)
    end

    return solutions

end

function problem_setup(borg::Borg)
    # create problem
    ref = ccall((:BORG_Problem_create), Ptr{Cvoid}, (Int64, Int64, Int64, Ptr{Cvoid}), borg.nVars, borg.nObjs, borg.nCons, borg.problem)

    # set bounds
    setBounds(borg.bounds, borg.nVars, ref)

    # set epsilons
    setEpsilons(borg.epsilons, borg.nObjs, ref)

    return ref
end


function run(borg::Borg, settings::Dict=Dict())
    # create problem
    ref = problem_setup(borg)

    start = now()

    # PM
    pm_f = dlsym(config.lib, :BORG_Operator_PM)
    pm_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "PM", 1, 1, 2, pm_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), pm_op, 0, get(settings, "pm.rate", 1 / borg.nVars))
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), pm_op, 1, get(settings, "pm.distributionIndex", 20.0))

    # SBX
    sbx_f = dlsym(config.lib, :BORG_Operator_SBX)
    sbx_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "SBX", 2, 2, 2, sbx_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), sbx_op, 0, get(settings, "sbx.rate", 1.0))
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), sbx_op, 1, get(settings, "sbx.distributionIndex", 15.0))
    ccall((:BORG_Operator_set_mutation), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), sbx_op, pm_op)

    # DE 
    de_f = dlsym(config.lib, :BORG_Operator_DE)
    de_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "DE", 4, 1, 2, de_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), de_op, 0, get(settings, "de.crossoverRate", 0.1))
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), de_op, 1, get(settings, "de.stepSize", 0.5))
    ccall((:BORG_Operator_set_mutation), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), de_op, pm_op)

    # UM
    um_f = dlsym(config.lib, :BORG_Operator_UM)
    um_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "UM", 1, 1, 1, um_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), um_op, 0, get(settings, "um.rate", 1 / borg.nVars))

    # SPX 
    spx_f = dlsym(config.lib, :BORG_Operator_SPX)
    spx_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "SPX", get(settings, "spx.parents", 10), get(settings, "spx.offspring", 2), 1, spx_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), spx_op, 0, get(settings, "um.epsilon", 3.0))

    # PCX 
    pcx_f = dlsym(config.lib, :BORG_Operator_PCX)
    pcx_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "PCX", get(settings, "pcx.parents", 10), get(settings, "pcx.offspring", 2), 2, pcx_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), pcx_op, 0, get(settings, "pcx.eta", 0.1))
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), pcx_op, 1, get(settings, "pcx.zeta", 0.1))

    # UNDX
    undx_f = dlsym(config.lib, :BORG_Operator_UNDX)
    undx_op = ccall((:BORG_Operator_create), Ptr{Cvoid}, (Cstring, Int32, Int32, Int32, Ptr{Nothing}), "UNDX", get(settings, "undx.parents", 10), get(settings, "undx.offspring", 2), 2, undx_f)
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), undx_op, 0, get(settings, "undx.zeta", 0.5))
    ccall((:BORG_Operator_set_parameter), Cvoid, (Ptr{Cvoid}, Int32, Float64), undx_op, 1, get(settings, "undx.eta", 0.35))

    ## set up algorithm 
    algorithm = ccall((:BORG_Algorithm_create), Ptr{Cvoid}, (Ptr{Cvoid}, Int32), ref, 6)
    ccall((:BORG_Algorithm_set_operator), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Cvoid}), algorithm, 0, sbx_op)
    ccall((:BORG_Algorithm_set_operator), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Cvoid}), algorithm, 1, de_op)
    ccall((:BORG_Algorithm_set_operator), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Cvoid}), algorithm, 2, pcx_op)
    ccall((:BORG_Algorithm_set_operator), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Cvoid}), algorithm, 3, spx_op)
    ccall((:BORG_Algorithm_set_operator), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Cvoid}), algorithm, 4, undx_op)
    ccall((:BORG_Algorithm_set_operator), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Cvoid}), algorithm, 5, um_op)

    fp = nothing
    frequency = nothing
    last_snapshot = nothing
    write_runtime = false
    
    if (haskey(settings, "runtimefile"))
        fp = open(get(settings, "runtimefile", nothing), "w")
        frequency = get(settings, "frequency", borg.nVars / 100)
        last_snapshot = 0
        write_runtime = true
    end


    # algorithm run
    while (local evals = ccall((:BORG_Algorithm_get_nfe), Int32, (Ptr{Cvoid},), algorithm)) < borg.nfe
        ccall((:BORG_Algorithm_step), Cvoid, (Ptr{Cvoid},), algorithm)

        if write_runtime & (evals - last_snapshot >= frequency)
            entry = Dict()
            entry["NFE"] = evals
            entry["ElapsedTime"] = now() - start
            entry["SBX"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), sbx_op)
            entry["DE"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), de_op)
            entry["PCX"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), pcx_op)
            entry["SPX"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), spx_op)
            entry["UNDX"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), undx_op)
            entry["UM"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), um_op)
            entry["SBX"] = ccall((:BORG_Operator_get_probability), Cdouble, (Ptr{Cvoid},), sbx_op)
            entry["Improvements"] = ccall((:BORG_Algorithm_get_number_improvements), Int32, (Ptr{Cvoid},), algorithm)
            entry["Restarts"] = ccall((:BORG_Algorithm_get_number_restarts), Int32, (Ptr{Cvoid},), algorithm)
            entry["PopulationSize"] = ccall((:BORG_Algorithm_get_population_size), Int32, (Ptr{Cvoid},), algorithm)
            entry["ArchiveSize"] = ccall((:BORG_Algorithm_get_archive_size), Int32, (Ptr{Cvoid},), algorithm)

            archive = process(ccall((:BORG_Algorithm_get_result), Ptr{Cvoid}, (Ptr{Cvoid},), algorithm), borg)
            metrics = ["NFE", "ElapsedTime", "SBX", "DE", "PCX", "SPX", "UNDX", "UM", "Improvements", "Restarts", "PopulationSize", "ArchiveSize"]

            for metric in metrics
                write(fp, "//$metric=$(entry[metric])\n")
            end

            for sol in archive
                for var in sol.vars
                    write(fp, "$var ")
                end
                for obj in sol.objs
                    write(fp, "$obj ")
                end
                for constr in sol.constrs
                    write(fp, "$constr ")
                end
                write(fp, "\n")
            end


            last_snapshot = evals
        end
    end
    
    # process results
    result = ccall((:BORG_Algorithm_get_result), Ptr{Cvoid}, (Ptr{Cvoid},), algorithm)
    solutions = process(result, borg)

    # close file 
    if (fp != nothing)
        close(fp)
    end

    # free memory
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), sbx_op)
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), de_op)
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), pm_op)
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), um_op)
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), spx_op)
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), pcx_op)
    ccall((:BORG_Operator_destroy), Cvoid, (Ptr{Cvoid},), undx_op)
    ccall((:BORG_Algorithm_destroy), Cvoid, (Ptr{Cvoid},), algorithm)

    return solutions
end

const mpi_dylib_path = "borg/libborg.so"
const borgms = "borg/libborgms.so"
    
function runMPI(borg::Borg, settings::Dict=Dict())
    mpi_lib = dlopen(mpi_dylib_path, RTLD_GLOBAL)
    borg_lib = dlopen(borgms, RTLD_GLOBAL)

    ccall((:BORG_Algorithm_ms_startup, borgms), Cvoid, (Ptr{Cint}, Ptr{Cint}), C_NULL, C_NULL)
    ccall((:BORG_Algorithm_ms_max_evaluations, borgms), Cvoid, (Int32,), 1000)
    ref = problem_setup(borg)

    result = ccall((:BORG_Algorithm_ms_run, borgms), Ptr{Cvoid}, (Ptr{Cvoid},), ref)
    solutions = process(result, borg)
    println(solutions)

    ccall((:BORG_Algorithm_ms_shutdown, borgms), Cvoid, ())
end