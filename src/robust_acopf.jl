
function solveRobustACOPF(Scen::Array{Scenario},
                          pertLoads::Set{String},
                          pert::Float64
                          )
    println("> COMPUTE ROBUST SOLUTION")
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    alpha=DefAlpha(pertLoads,pert)
    omegaMax=sum(pertPUp[i]-pd[i] for i=pertLoads)
    omegaMin=sum(pd[i]-pertPDown[i] for i=pertLoads)
    p0min=Dict{String,Float64}(pmin)
    p0max=Dict{String,Float64}(pmax)
    for i=NonSlackGen
        p0min[i]+=alpha[i]*omegaMin
        p0max[i]-=alpha[i]*omegaMax
    end

    m=JuMP.Model()
    ScID=1:length(Scen)
    @variable(m, pmin[i] <= varpg[i=keys(netgen),ScID] <= pmax[i])
    @variable(m, p0min[i] <= p0[i=keys(netgen)] <= p0max[i])
    @variable(m, qmin[i] <= varqg[i=keys(netgen),ScID] <= qmax[i])
    @variable(m, varva[keys(netbus),ScID])
    @variable(m, vmin[i]<= varvm[i=keys(netbus),ScID] <= vmax[i])
    @variable(m, vgen[i=keys(netgen)])
    @variable(m, -rate_a[i] <= varfp_fr[i=keys(netbranch),ScID] <= rate_a[i])
    @variable(m, -rate_a[i] <= varfp_to[i=keys(netbranch),ScID] <= rate_a[i])
    @variable(m, -rate_a[i] <= varfq_fr[i=keys(netbranch),ScID] <= rate_a[i])
    @variable(m, -rate_a[i] <= varfq_to[i=keys(netbranch),ScID] <= rate_a[i])

    @constraint(m, DefpgRef[i=RefGen],
                p0[i]==sum(varpg[i,s] for s=ScID)/length(Scen))

    @constraint(m, Defvgen[i=keys(netgen),s=ScID],
                vgen[i]==varvm[string(netgen[i]["gen_bus"]),s])

    @constraint(
        m, Defpg[i=NonSlackGen, s=ScID],
        varpg[i,s]
        == p0[i] + alpha[i]*sum(values(Scen[s].loadP))
    )

    @constraint(m, RefAngle[i=ref_bus,s=ScID],
                varva[i,s]==network["bus"][i]["va"])

    @NLconstraint(
        m, Deffp_fr[i=keys(netbranch),s=ScID],
        varfp_fr[i,s] == (g[i]+g_fr[i])*(varvm[f_bus[i],s]/tap[i])^2
        - g[i]*varvm[f_bus[i],s]/tap[i]*varvm[t_bus[i],s]*(
            cos(varva[f_bus[i],s]-varva[t_bus[i],s]-shift[i]))
        - b[i]*varvm[f_bus[i],s]/tap[i]*varvm[t_bus[i],s]*(
            sin(varva[f_bus[i],s]-varva[t_bus[i],s]-shift[i]))
    )
    @NLconstraint(
        m, Deffq_fr[i=keys(netbranch),s=ScID],
        varfq_fr[i,s] == -(b[i]+b_fr[i])*(varvm[f_bus[i],s]/tap[i])^2
        + b[i]*varvm[f_bus[i],s]/tap[i]*varvm[t_bus[i],s]*(
            cos(varva[f_bus[i],s]-varva[t_bus[i],s]-shift[i]))
        - g[i]*varvm[f_bus[i],s]/tap[i]*varvm[t_bus[i],s]*(
            sin(varva[f_bus[i],s]-varva[t_bus[i],s]-shift[i]))
    )
    @NLconstraint(
        m, Deffp_to[i=keys(netbranch),s=ScID],
        varfp_to[i,s] == (g[i]+g_to[i])*(varvm[t_bus[i],s])^2
        - g[i]*varvm[t_bus[i],s]*varvm[f_bus[i],s]/tap[i]*(
            cos(varva[t_bus[i],s]-varva[f_bus[i],s]+shift[i]))
        - b[i]*varvm[t_bus[i],s]*varvm[f_bus[i],s]/tap[i]*(
            sin(varva[t_bus[i],s]-varva[f_bus[i],s]+shift[i]))
    )
    @NLconstraint(
        m, Deffq_to[i=keys(netbranch),s=ScID],
        varfq_to[i,s] == -(b[i]+b_to[i])*(varvm[t_bus[i],s])^2
        + b[i]*varvm[t_bus[i],s]*varvm[f_bus[i],s]/tap[i]*(
            cos(varva[t_bus[i],s]-varva[f_bus[i],s]+shift[i]))
        - g[i]*varvm[t_bus[i],s]*varvm[f_bus[i],s]/tap[i]*(
            sin(varva[t_bus[i],s]-varva[f_bus[i],s]+shift[i]))
    )
    @constraint(
        m, LimFlowFr[i=keys(netbranch),s=ScID],
        varfp_fr[i,s]^2+varfq_fr[i,s]^2 <= rate_a[i]^2
    )
    @constraint(
        m, LimFlowTo[i=keys(netbranch),s=ScID],
        varfp_to[i,s]^2+varfq_to[i,s]^2 <= rate_a[i]^2
    )

    @constraint(
        m, RealPowerBalance[i=keys(netbus),s=ScID],
        + sum(varfp_fr[b,s] for b=branch_f[i])
        + sum(varfp_to[b,s] for b=branch_t[i])
        ==
        + sum(varpg[g,s]
              for g=gen_i[i])
        - sum(pd[d]+Scen[s].loadP[d] for d=load_i[i])
        - sum(netshunt[s]["gs"] for s=shunt_i[i])*varvm[i,s]^2
    )

    @constraint(
        m, ReactivePowerBalance[i=keys(netbus),s=ScID],
        + sum(varfq_fr[b,s] for b=branch_f[i])
        + sum(varfq_to[b,s] for b=branch_t[i])
        ==
        + sum(varqg[g,s] for g=gen_i[i])
        - sum(qd[d]+Scen[s].loadQ[d] for d=load_i[i])
        + sum(netshunt[s]["bs"] for s=shunt_i[i])*varvm[i,s]^2
    )

    @objective(
        m, Min,
        sum(netgen[g]["cost"][1]*p0[g] + netgen[g]["cost"][2]
            for g=keys(netgen) if length(netgen[g]["cost"]) == 2)
        + sum(netgen[g]["cost"][1]*p0[g]^2 + netgen[g]["cost"][2]*p0[g]
              + netgen[g]["cost"][3]
              for g=keys(netgen) if length(netgen[g]["cost"]) == 3)
    )

    optimize!(m, with_optimizer(Ipopt.Optimizer))

    PG0=Dict{String, Float64}(i=>value(p0[i]) for i=NonSlackGen)
    V0=Dict{String, Float64}(i=>value(varvm[i,1])
                             for i=keys(netbus) if netbus[i]["bus_type"]!=1)

    return PG0,V0
end

function solveDeterministic(pertLoads::Set{String},
                            pert::Float64
                            )
    return solveRobustACOPF([CentralScenario()],pertLoads,pert)
end

function solveOPF()
    res=run_ac_opf(network,with_optimizer(Ipopt.Optimizer))
    sol=res["solution"]
    PG0=Dict{String, Float64}(i=>sol["gen"][i]["pg"] for i=NonSlackGen)
    V0=Dict{String, Float64}(i=>sol["bus"][i]["vm"]
                             for i=keys(netbus) if netbus[i]["bus_type"]!=1)
    return PG0, V0
end

function distance_to_opf(PG0,V0)
    PG00,V00=solveOPF()
    return (sqrt(sum(abs2,PG0[i]-PG00[i] for i=keys(PG0))),
            sqrt(sum(abs2,V0[i]-V00[i] for i=keys(V0))))
end

function GenCost(res)
    sol=res["solution"]["gen"]
    return  (mapreduce(sum,+,netgen[g]["cost"][1]*sol[g]["pg"]
                       + netgen[g]["cost"][2]
                       for g=keys(netgen) if length(netgen[g]["cost"]) == 2
                       ;init=0)
             + mapreduce(sum,+,netgen[g]["cost"][1]*sol[g]["pg"]^2
                         + netgen[g]["cost"][2]*sol[g]["pg"]
                         + netgen[g]["cost"][3]
                         for g=keys(netgen) if length(netgen[g]["cost"]) == 3
                         ;init=0)
             + mapreduce(sum,+, netgen[g]["cost"][1]
                         for g=keys(netgen) if length(netgen[g]["cost"]) == 1
                         ;init=0))
end

function SampleScenario(Scen::Scenario,
                        PG0::Dict{String,Float64},V0::Dict{String,Float64},
                        pertLoads::Set{String}, pert::Float64)
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    alpha=DefAlpha(pertLoads,pert)
    for i=NonSlackGen
        network["gen"][i]["pg"]=PG0[i]
    end
    for i=keys(netgen)
        network["gen"][i]["vg"]=V0[string(netgen[i]["gen_bus"])]
    end
    for i=keys(netbus)
        if netbus[i]["bus_type"]!=1
            network["bus"][i]["vm"]=V0[i]
        end
    end
    for i=pertLoads
        network["load"][i]["pd"]+=Scen.loadP[i]
        network["load"][i]["qd"]+=Scen.loadQ[i]
    end
    w=sum(network["load"][i]["pd"] for i=keys(netload)) - sum(values(pd))
    for i=NonSlackGen
        network["gen"][i]["pg"]=PG0[i]+alpha[i]*w
    end
    results=run_ac_pf(network,with_optimizer(Ipopt.Optimizer,
                                             print_level=0,
                                             max_iter=250))
    if results["termination_status"] != LOCALLY_SOLVED
        println(">> Issue in power flow recourse")
    end
    ConstrViol,ViolConstr=checkSolution(results,true)
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    println("Cost of Sample: $(GenCost(results))")
    return results
end


function SampleScenarios(Scens::Array{Scenario},
                  PG0::Dict{String,Float64},V0::Dict{String,Float64},
                  pertLoads::Set{String}, pert::Float64)
    println("> SAMPLING : $(length(Scens)) samples")
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    alpha=DefAlpha(pertLoads,pert)
    passedScen=Array{Scenario}(undef,0)
    failedScen=Array{Scenario}(undef,0)
    ConstrViolFailedScen=Array{Array{String}}(undef,0)
    MaxViolFailedScen=Array{Float64}(undef,0)
    ConstrViolToScen=Dict{String,Array{Tuple{Float64,Scenario,Scenario}}}()
    ConstrMaxViol=Array{String}(undef,0)
    results=Any
    nbInfeasible=0
    nbPFInfeasible=0
    maxViol=0.0
    for i=NonSlackGen
        network["gen"][i]["pg"]=PG0[i]
    end
    for i=keys(netgen)
        network["gen"][i]["vg"]=V0[string(netgen[i]["gen_bus"])]
    end
    for i=keys(netbus)
        if netbus[i]["bus_type"]!=1
            network["bus"][i]["vm"]=V0[i]
        end
    end
    for s=1:length(Scens)
        if s % 100 == 0 println("s=",s) end
        for i=pertLoads
            network["load"][i]["pd"]=Scens[s].loadP[i]
            network["load"][i]["qd"]=Scens[s].loadQ[i]
        end
        Scen=Scenario(
            Dict{String,Float64}(i=>network["load"][i]["pd"]-pd[i]
                                 for i=keys(netload)),
            Dict{String,Float64}(i=>network["load"][i]["qd"]-qd[i]
                                 for i=keys(netload))
        )
        w=sum(network["load"][i]["pd"] for i=keys(netload)) - sum(values(pd))
        for i=NonSlackGen
            network["gen"][i]["pg"]=PG0[i]+alpha[i]*w
        end
        results=run_ac_pf(network,with_optimizer(Ipopt.Optimizer,
                                                 print_level=0,
                                                 max_iter=250))
        if results["termination_status"] != LOCALLY_SOLVED
            nbInfeasible+=1
            continue
        end
        ConstrViol,ViolConstr=checkSolution(results,false)
        if  length(ConstrViol) > 0
            maxViol=max(maxViol,maximum(ViolConstr))
            #println("omega:",w)
            push!(failedScen,Scen)
            push!(ConstrViolFailedScen,ConstrViol)
            push!(MaxViolFailedScen,maximum(ViolConstr))
            push!(ConstrMaxViol,ConstrViol[findfirst((j->j==maximum(ViolConstr)),
                                          ViolConstr)])
            for c=1:length(ConstrViol)
                push!(get!(ConstrViolToScen,ConstrViol[c],
                           Array{Tuple{Float64,Scenario,Scenario}}(undef,0)),
                      (ViolConstr[c],Scen,CentralScenario())
                      )
            end
            nbPFInfeasible+=1
        else
            push!(passedScen,Scen)
        end
    end
    # println("There are ", nbPFInfeasible,
    #         " infeasible PF solutions violating OPF constraints out of ",
    #         NbSamples, ".")
    println("Max violation of constraints: ", maxViol*100, "%")
    if !isempty(failedScen)
        println("Expected max violation when failing: ",
                Statistics.mean(MaxViolFailedScen)*100, "%" )
        println("Expected number of violated constraints when failing: ",
                ceil(Int64,Statistics.mean(length(i)
                                           for i=ConstrViolFailedScen)))
    end
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    return (nbPFInfeasible,maxViol*100,passedScen,failedScen,
            ConstrViolFailedScen,MaxViolFailedScen,ConstrViolToScen,
            ConstrMaxViol)
end

function SamplePF(NbSamples::Int64,
                  PG0::Dict{String,Float64},V0::Dict{String,Float64},
                  pertLoads::Set{String}, pert::Float64)
    println("> SAMPLING : $NbSamples samples")
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    alpha=DefAlpha(pertLoads,pert)
    passedScen=Array{Scenario}(undef,0)
    failedScen=Array{Scenario}(undef,0)
    ConstrViolFailedScen=Array{Array{String}}(undef,0)
    MaxViolFailedScen=Array{Float64}(undef,0)
    ConstrViolToScen=Dict{String,Array{Tuple{Float64,Scenario,Scenario}}}()
    ConstrMaxViol=Array{String}(undef,0)
    results=Any
    nbInfeasible=0
    nbPFInfeasible=0
    maxViol=0.0
    for i=NonSlackGen
        network["gen"][i]["pg"]=PG0[i]
    end
    for i=keys(netgen)
        network["gen"][i]["vg"]=V0[string(netgen[i]["gen_bus"])]
    end
    for i=keys(netbus)
        if netbus[i]["bus_type"]!=1
            network["bus"][i]["vm"]=V0[i]
        end
    end
    for s=1:NbSamples
        if s % 100 == 0 println("s=",s) end
        for i=pertLoads
            network["load"][i]["pd"]=(pertPDown[i]
                                      +rand()*(pertPUp[i]-pertPDown[i]))
            network["load"][i]["qd"]=(pertQDown[i]
                                      +rand()*(pertQUp[i]-pertQDown[i]))
        end
        Scen=Scenario(
            Dict{String,Float64}(i=>network["load"][i]["pd"]-pd[i]
                                 for i=keys(netload)),
            Dict{String,Float64}(i=>network["load"][i]["qd"]-qd[i]
                                 for i=keys(netload))
        )
        w=sum(network["load"][i]["pd"] for i=keys(netload)) - sum(values(pd))
        for i=NonSlackGen
            network["gen"][i]["pg"]=PG0[i]+alpha[i]*w
        end
        results=run_ac_pf(network,with_optimizer(Ipopt.Optimizer,
                                                 print_level=0,
                                                 max_iter=250))
        if results["termination_status"] != LOCALLY_SOLVED
            nbInfeasible+=1
            continue
        end
        ConstrViol,ViolConstr=checkSolution(results,false)
        if  length(ConstrViol) > 0
            maxViol=max(maxViol,maximum(ViolConstr))
            #println("omega:",w)
            push!(failedScen,Scen)
            push!(ConstrViolFailedScen,ConstrViol)
            push!(MaxViolFailedScen,maximum(ViolConstr))
            push!(ConstrMaxViol,ConstrViol[findfirst((j->j==maximum(ViolConstr)),
                                          ViolConstr)])
            for c=1:length(ConstrViol)
                push!(get!(ConstrViolToScen,ConstrViol[c],
                           Array{Tuple{Float64,Scenario}}(undef,0)),
                      (ViolConstr[c],Scen,CentralScenario())
                      )
            end
            nbPFInfeasible+=1
        else
            push!(passedScen,Scen)
        end
    end
    println("There are ", nbPFInfeasible,
            " infeasible PF solutions violating OPF constraints out of ",
            NbSamples, ".")
    println("Max violation of constraints: ", maxViol*100, "%")
    if !isempty(failedScen)
        println("Expected max violation when failing: ",
                Statistics.mean(MaxViolFailedScen)*100, "%" )
        println("Expected number of violated constraints when failing: ",
                ceil(Int64,Statistics.mean(length(i)
                                           for i=ConstrViolFailedScen)))
    end
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    return (nbPFInfeasible,maxViol*100,passedScen,failedScen,
            ConstrViolFailedScen,MaxViolFailedScen,ConstrViolToScen,
            ConstrMaxViol)
end

function checkSolution(res,verbosity::Bool)
    gen=res["solution"]["gen"]
    bus=res["solution"]["bus"]

    pg=Dict{String,Float64}(i => gen[i]["pg"] for i=keys(gen))
    qg=Dict{String,Float64}(i => gen[i]["qg"] for i=keys(gen))
    vm=Dict{String,Float64}(i => bus[i]["vm"] for i=keys(bus))
    va=Dict{String,Float64}(i => bus[i]["va"] for i=keys(bus))

    fp_fr=Dict{String,Float64}(
        i=> (g[i]+g_fr[i])*(vm[f_bus[i]]/tap[i])^2
        - g[i]*vm[f_bus[i]]/tap[i]*vm[t_bus[i]]*(
            cos(va[f_bus[i]]-va[t_bus[i]]-shift[i]))
        - b[i]*vm[f_bus[i]]/tap[i]*vm[t_bus[i]]*(
            sin(va[f_bus[i]]-va[t_bus[i]]-shift[i]))
        for i=keys(netbranch))
    fq_fr=Dict{String,Float64}(
        i=> -(b[i]+b_fr[i])*(vm[f_bus[i]]/tap[i])^2
        + b[i]*vm[f_bus[i]]/tap[i]*vm[t_bus[i]]*(
            cos(va[f_bus[i]]-va[t_bus[i]]-shift[i]))
        - g[i]*vm[f_bus[i]]/tap[i]*vm[t_bus[i]]*(
            sin(va[f_bus[i]]-va[t_bus[i]]-shift[i]))
        for i=keys(netbranch))
    fp_to=Dict{String,Float64}(
        i=> (g[i]+g_to[i])*(vm[t_bus[i]])^2
        - g[i]*vm[t_bus[i]]*vm[f_bus[i]]/tap[i]*(
            cos(va[t_bus[i]]-va[f_bus[i]]+shift[i]))
        - b[i]*vm[t_bus[i]]*vm[f_bus[i]]/tap[i]*(
            sin(va[t_bus[i]]-va[f_bus[i]]+shift[i]))
        for i=keys(netbranch))
    fq_to=Dict{String,Float64}(
        i=> -(b[i]+b_to[i])*(vm[t_bus[i]])^2
        + b[i]*vm[t_bus[i]]*vm[f_bus[i]]/tap[i]*(
            cos(va[t_bus[i]]-va[f_bus[i]]+shift[i]))
        - g[i]*vm[t_bus[i]]*vm[f_bus[i]]/tap[i]*(
            sin(va[t_bus[i]]-va[f_bus[i]]+shift[i]))
        for i=keys(netbranch))

    VDown=Set{String}()
    VUp=Set{String}()
    PgDown=Set{String}()
    QgDown=Set{String}()
    PgUp=Set{String}()
    QgUp=Set{String}()
    AngUp=Set{String}()
    AngDown=Set{String}()
    FlowLimFr=Set{Any}()
    FlowLimTo=Set{Any}()
    ConstrViol=Array{String}(undef,0)
    ViolConstr=Array{Float64}(undef,0)
    for i=keys(bus)
        if mapreduce(sum,+,pg[g] - pmax[g] for g=gen_i[i];init=0) > tol
            push!(PgUp, i)
            push!(ConstrViol,"PgUp$i")
            push!(ViolConstr,
                  abs(mapreduce(sum,+,pg[g] - pmax[g] for g=gen_i[i];init=0))/(
                      abs(mapreduce(sum,+,pmax[g]-pmin[g] for g=gen_i[i];init=1)))
                  )
        end
        if mapreduce(sum,+,pmin[g] - pg[g] for g=gen_i[i];init=0) > tol
            push!(PgDown, i)
            push!(ConstrViol,"PgDown$i")
            push!(ViolConstr,
                  abs(mapreduce(sum,+,pg[g] - pmin[g] for g=gen_i[i];init=0))/(
                      abs(mapreduce(sum,+,pmax[g]-pmin[g] for g=gen_i[i];init=1)))
                  )
        end
        if mapreduce(sum,+,qg[g] - qmax[g] for g=gen_i[i];init=0) > tol
            push!(QgUp, i)
            push!(ConstrViol,"QgUp$i")
            push!(ViolConstr,
                  abs(mapreduce(sum,+,qg[g] - qmax[g] for g=gen_i[i];init=0))/(
                      abs(mapreduce(sum,+,qmax[g]-qmin[g] for g=gen_i[i];init=1)))
                  )
        end
        if mapreduce(sum,+,qmin[g] - qg[g] for g=gen_i[i];init=0) > tol
            push!(QgDown, i)
            push!(ConstrViol,"QgDown$i")
            push!(ViolConstr,
                  abs(mapreduce(sum,+,qg[g] - qmin[g] for g=gen_i[i];init=0))/(
                      abs(mapreduce(sum,+,qmax[g]-qmin[g] for g=gen_i[i];init=1)))
                  )
        end
        if vm[i] - vmax[i] > tol
            push!(VUp, i)
            push!(ConstrViol,"VUp$i")
            push!(ViolConstr,abs(vm[i]-vmax[i])/abs(vmax[i]-vmin[i]))
        end
        if vmin[i] - vm[i] > tol
            push!(VDown, i)
            push!(ConstrViol,"VDown$i")
            push!(ViolConstr,abs(vm[i]-vmin[i])/abs(vmax[i]-vmin[i]))
        end
    end
    for b=keys(netbranch)
        if fp_fr[b]^2+fq_fr[b]^2 - rate_a[b]^2 > tol
            push!(FlowLimFr,[f_bus[b],t_bus[b]])
            push!(ConstrViol,"FlowLimFr$b")
            push!(ViolConstr,(fp_fr[b]^2+fq_fr[b]^2-rate_a[b]^2)/rate_a[b]^2)
        end
        if fp_to[b]^2+fq_to[b]^2 - rate_a[b]^2 > tol
            push!(FlowLimTo,[f_bus[b],t_bus[b]])
            push!(ConstrViol,"FlowLimTo$b")
            push!(ViolConstr,(fp_to[b]^2+fq_to[b]^2-rate_a[b]^2)/rate_a[b]^2)
        end
    end
    nbConstrViol=length(ConstrViol)
    if verbosity && nbConstrViol > 0
        maxViol=maximum(ViolConstr)
        println("################################################")
        println("### Checking solution ###")
        println("# Real power generation violations: ",
                length(PgDown)+length(PgUp), ", ", union(PgDown, PgUp))
        println("# Reactive power generation violations: ",
                length(QgUp) + length(QgDown), ", ", union(QgDown, QgUp))
        println("# Voltage violations: ", length(VDown)+length(VUp),
                ", ", union(VDown, VUp))
        println("# Voltage angle violations: ", length(AngUp) + length(AngDown))
        println("# Flow limit violations: ", length(union(FlowLimFr,FlowLimTo)),
                ", ", union(FlowLimFr,FlowLimTo))
        println("Max violation of constraints: ", maxViol*100, "%")
        println("Number of constraints violated: ", nbConstrViol)
    end
    return ConstrViol,ViolConstr
end



function IterateRACOPF(Clusters::Array{Scenario},
                       nbSamples::Int64,nbScenToAdd::Int64,
                       maxNbScen::Int64, pertLoads::Set{String}, pert::Float64,
                       method::String,pushing::String
                       )
    NbIterations=0
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    PG0=Dict{String,Float64}
    V0=Dict{String,Float64}
    while length(Clusters) <= maxNbScen
        PG0,V0=solveRobustACOPF(Clusters,pertLoads,pert)
        println("\n\n$(length(Clusters)) scenarios in scenario set.")
        res=SamplePF(nbSamples,PG0,V0,pertLoads,pert)
        if length(res[4])==0
            println(">> SUCCESS: It's all done!")
            return PG0,V0,NbIterations
        end
        if length(Clusters)==maxNbScen
            println("\n\n$(length(Clusters)) (=Max!) scenarios",
                    " in scenario set.")
            return PG0,V0,NbIterations
        end
        if method=="Random"
            for i=1:min(nbScenToAdd,maxNbScen-length(Clusters))
                push!(Clusters, Scenario(
                    Dict{String,Float64}(i=>pertPDown[i]-pd[i]
                                         +rand()*(pertPUp[i]-pertPDown[i])
                                         for i=keys(netload)),
                    Dict{String,Float64}((i=>pertQDown[i]-qd[i]
                                          +rand()*(pertQUp[i]-pertQDown[i])
                                          for i=keys(netload))))
                      )
            end
            continue
        end
        hashConstrViol=[(hash((res[5][i],res[6][i])),hash(res[5][i]),
                         length(res[5][i]),res[6][i]) for i=1:length(res[5])]
        setHashConstrViol=Set{UInt64}([i[2] for i=hashConstrViol])
        println("There are $(length(setHashConstrViol)) different scenario ",
                "sets out of $(length(res[5])) failing scenarios")
        maxMaxViol=maximum(i[4] for i=hashConstrViol)
        maxNbViol=maximum(i[3] for i=hashConstrViol)
        countHashConstrViol=[(i[4]/maxMaxViol+i[3]/maxNbViol,i[3],
                              count(j->(j[2]==i[2]),hashConstrViol),i[1],i[2])
                             for i=hashConstrViol]
        if method=="NbConstr"
            countHashConstrViol=[(i[3],i[4],
                                  count(j->(j[2]==i[2]),hashConstrViol),
                                  i[1],i[2]) for i=hashConstrViol]
        end
        if method=="MaxViol"
            countHashConstrViol=[(i[4],i[3],
                                  count(j->(j[2]==i[2]),hashConstrViol),
                                  i[1],i[2]) for i=hashConstrViol]
        end
        countHCV=sort([i for i=countHashConstrViol], rev=true)
        ScenAddHashKey=Set{UInt64}()
        i=1
        p=1
        println("> SCENARIO SELECTION")
        while (i <= nbScenToAdd && i <= length(setHashConstrViol))
            if length(Clusters)==maxNbScen
                println("\n\n$(length(Clusters)) (=Max!) scenarios",
                        " in scenario set.")
                PG0,V0=solveRobustACOPF(Clusters,pertLoads,pert)
                res=SamplePF(nbSamples,PG0,V0,pertLoads,pert)
                return PG0,V0,NbIterations
            end
            if countHCV[p][5] in ScenAddHashKey
                p+=1
                continue
            end
            push!(ScenAddHashKey,countHCV[p][5])
            hashKey=countHCV[p][4]
            println(countHCV[p])
            for j=1:length(res[5])
                if hashKey==hash((res[5][j],res[6][j]))
                    if pushing=="ExtremePush"
                        push!(Clusters,PutScenExtreme(res[4][j],pertLoads,pert))
                    elseif pushing=="GradPush"
                        push!(Clusters,
                              GradPush(res[4][j],res[7],PG0,V0,
                                       pertLoads,pert))
                    else
                        push!(Clusters, res[4][j])
                    end
                    continue
                end
            end
            i+=1
            p+=1
        end
        NbIterations+=1
    end
    #println(">> FINISHED: $(length(Clusters)) > $(maxNbScen).")
    return PG0,V0,NbIterations
end
