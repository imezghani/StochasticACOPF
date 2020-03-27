

function PutScenExtreme(Scen::Scenario,pertLoads::Set{String},pert::Float64
                        )::Scenario
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    ExtrScen=Scenario(Dict{String,Float64}(i => 0.0 for i=keys(netload)),
                      Dict{String,Float64}(i => 0.0 for i=keys(netload)))
    for i=keys(netload)
        if abs(Scen.loadP[i]) > 0.1*(pertPUp[i]-pertPDown[i])
            if Scen.loadP[i] > 0
                ExtrScen.loadP[i]=pertPUp[i]-pd[i]
            else
                ExtrScen.loadP[i]=pertPDown[i]-pd[i]
            end
        end
        if abs(Scen.loadQ[i]) > 0.1*(pertQUp[i]-pertQDown[i])
            if Scen.loadQ[i] > 0
                ExtrScen.loadQ[i]=pertQUp[i]-qd[i]
            else
                ExtrScen.loadQ[i]=pertQDown[i]-qd[i]
            end
        end
    end
    return ExtrScen
end


function ComputeWScen(CVTS::Dict{String,Array{Tuple{Float64,Scenario,Scenario}}},
                      Constr::String,pertLoads::Set{String})::Scenario
    if sum(abs,values(CVTS[Constr][1][3].loadP))!=0 return CVTS[Constr][1][3] end
    keysWS=union([(i,"p") for i=pertLoads],[(i,"q") for i=pertLoads])
    #m=Model(with_optimizer(Ipopt.Optimizer,print_level=0))
    m=Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0))
    @variable(m,wscen[keysWS])
    @variable(m,ws0)
    @variable(m,lambda[keysWS]>=0)
    @constraint(m,defLambda1[i=keysWS], wscen[i]  <= lambda[i])
    @constraint(m,defLambda2[i=keysWS], -wscen[i] <= lambda[i])
    @objective(m,Min,
               sum((100*y[1] - (ws0 + sum(y[2].loadP[i]*wscen[(i,"p")]
                                      + y[2].loadQ[i]*wscen[(i,"q")]
                                      for i=pertLoads)))^2
                   for y=CVTS[Constr])/length(CVTS[Constr])
               +0.01*sum(lambda[i] for i=keysWS)/(sqrt(length(CVTS[Constr])))#*length(keysWS))
               )
    optimize!(m)
    CVTS[Constr][1]=(CVTS[Constr][1][1],CVTS[Constr][1][2],
                    Scenario(merge(Dict{String,Float64}(i=> value(wscen[(i,"p")])
                                               for i=pertLoads),
                          Dict{String,Float64}(i=> 0.0
                                               for i=setdiff(keys(netload),
                                                             pertLoads))),
                    merge(Dict{String,Float64}(i=> value(wscen[(i,"q")])
                                               for i=pertLoads),
                          Dict{String,Float64}(i=> 0.0
                                               for i=setdiff(keys(netload),
                                                             pertLoads)))
                    ))
    return CVTS[Constr][1][3]
end

function ComputeCVTS(CriticalScen::Array{Scenario},
                     PG0::Dict{String,Float64},V0::Dict{String,Float64},
                     pertLoads::Set{String},pert::Float64
                     )::Dict{String,Array{Tuple{Float64,Scenario,Scenario}}}
    alpha=DefAlpha(pertLoads,pert)
    CVTS=Dict{String,Array{Tuple{Float64,Scenario,Scenario}}}()
    for i=keys(netgen)
        network["gen"][i]["vg"]=V0[string(netgen[i]["gen_bus"])]
    end
    for i=keys(netbus)
        if netbus[i]["bus_type"]!=1
            network["bus"][i]["vm"]=V0[i]
        end
    end
    for scen=CriticalScen
        for i=pertLoads
            network["load"][i]["pd"]=pd[i] + scen.loadP[i]
            network["load"][i]["qd"]=qd[i] + scen.loadQ[i]
        end
        w=sum(scen.loadP[i] for i=keys(netload))
        for i=NonSlackGen
            network["gen"][i]["pg"]=PG0[i]+alpha[i]*w
        end
        results=run_ac_pf(network,with_optimizer(Ipopt.Optimizer,
                                                 print_level=0,
                                                 max_iter=250))
        if results["termination_status"] != LOCALLY_SOLVED
            println(">>> Issue when solving scenario ", s)
            continue
        end
        ConstrViol,ViolConstr=checkSolution(results, false)
        if  length(ConstrViol) > 0
            for c=1:length(ConstrViol)
                push!(get!(CVTS,ConstrViol[c],
                           Array{Tuple{Float64,Scenario,Scenario}}(undef,0)),
                      (ViolConstr[c],scen,CentralScenario())
                      )
            end
        end
    end
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    return CVTS
end

function GradPush(Scen::Scenario,
                  CVTS::Dict{String,Array{Tuple{Float64,Scenario,Scenario}}},
                  PG0,V0,
                  pertLoads::Set{String},
                  pert::Float64)
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    alpha=DefAlpha(pertLoads,pert)
    for i=keys(netgen)
        network["gen"][i]["vg"]=V0[string(netgen[i]["gen_bus"])]
    end
    for i=keys(netbus)
        if netbus[i]["bus_type"]!=1
            network["bus"][i]["vm"]=V0[i]
        end
    end
    for i=pertLoads
        network["load"][i]["pd"]=pd[i] + Scen.loadP[i]
        network["load"][i]["qd"]=qd[i] + Scen.loadQ[i]
    end
    w=sum(Scen.loadP[i] for i=keys(netload))
    for i=NonSlackGen
        network["gen"][i]["pg"]=PG0[i]+alpha[i]*w
    end
    results=run_ac_pf(network,with_optimizer(Ipopt.Optimizer,
                                             print_level=0,
                                             max_iter=250))
    if results["termination_status"] != LOCALLY_SOLVED
        println(">>> Issue when solving scenario ", Scen)
    end
    ConstrViol,ViolConstr=checkSolution(results,false)
    nbPushP=0
    nbPushQ=0
    PushScen=Scenario(copy(Scen.loadP),copy(Scen.loadQ))
    for c=ConstrViol
        WS=ComputeWScen(CVTS,c,pertLoads)
        maxWSP=max(1e-3,maximum(abs,values(WS.loadP)))
        maxWSQ=max(1e-3,maximum(abs,values(WS.loadQ)))
        for i=pertLoads
            if (abs(WS.loadP[i]) > maxWSP*1e-2
                && PushScen.loadP[i]==Scen.loadP[i])#tol
                nbPushP+=1
                if sign(WS.loadP[i]) > 0
                    PushScen.loadP[i]=pertPUp[i]-pd[i]
                else
                    PushScen.loadP[i]=pertPDown[i]-pd[i]
                end
            end
            if (abs(WS.loadQ[i]) > maxWSQ*1e-2
                && PushScen.loadQ[i]==Scen.loadQ[i])#tol
                nbPushQ+=1
                if sign(WS.loadQ[i]) > 0
                    PushScen.loadQ[i]=pertQUp[i]-qd[i]
                else
                    PushScen.loadQ[i]=pertQDown[i]-qd[i]
                end
            end
        end
    end
    println("Pushed P: $nbPushP / $(length(pertLoads))")
    println("Pushed Q: $nbPushQ / $(length(pertLoads))")
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    return PushScen
end
