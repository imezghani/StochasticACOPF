
function CentralScenario()::Scenario
    return Scenario(Dict{String,Float64}(i => 0.0
                                         for i=keys(netload)),
                    Dict{String,Float64}(i => 0.0
                                         for i=keys(netload)))
end


function ExtremeScenarios(pertLoads::Set{String}, pert::Float64)::Array{Scenario}
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    ExtremeScen=Array{Scenario}(undef,0)
    push!(ExtremeScen,
          Scenario(Dict{String,Float64}(i => pertPUp[i]-pd[i]
                                        for i=keys(netload)),
                   Dict{String,Float64}(i => pertQUp[i]-qd[i] 
                                        for i=keys(netload))))
    push!(ExtremeScen,
          Scenario(Dict{String,Float64}(i => pertPUp[i]-pd[i] 
                                        for i=keys(netload)),
                   Dict{String,Float64}(i => pertQDown[i]-qd[i] 
                                        for i=keys(netload))))
    push!(ExtremeScen,
          Scenario(Dict{String,Float64}(i => pertPDown[i]-pd[i] 
                                        for i=keys(netload)),
                   Dict{String,Float64}(i => pertQUp[i]-qd[i]
                                        for i=keys(netload))))
    push!(ExtremeScen,
          Scenario(Dict{String,Float64}(i => pertPDown[i]-pd[i] 
                                        for i=keys(netload)),
                   Dict{String,Float64}(i => pertQDown[i]-qd[i]
                                        for i=keys(netload))))
    return ExtremeScen
end


function OPF_constraints(PG0::Dict{String,Float64},
                         V0::Dict{String,Float64},
                         pertLoads::Set{String},pert::Float64)::Model
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    alpha=DefAlpha(pertLoads,pert)
    m=JuMP.Model()
    @variable(m, varpg[i=keys(netgen)])
    @variable(m, qmin0[i] <= varqg[i=keys(netgen)] <= qmax0[i])
    @variable(m, varva[keys(netbus)])
    @variable(m, vmin0[i]<= varvm[i=keys(netbus)] <= vmax0[i])
    @variable(m, varfp_fr[i=keys(netbranch)])
    @variable(m, varfp_to[i=keys(netbranch)])
    @variable(m, varfq_fr[i=keys(netbranch)])
    @variable(m, varfq_to[i=keys(netbranch)])
    @variable(m, pertPDown[i] <= varloadP[i=keys(netload)] <= pertPUp[i])
    @variable(m, pertQDown[i] <= varloadQ[i=keys(netload)] <= pertQUp[i])
    for i=NonSlackGen
        JuMP.set_start_value.(varpg[i],PG0[i])
    end
    for i=PVbus
        JuMP.set_start_value.(varvm[i],V0[i])
    end
    @constraint(m, FixV0[i=PVbus], varvm[i]==V0[i])    
    @constraint(
        m, DefPOmega[g=NonSlackGen],
        varpg[g]==PG0[g]
        + alpha[g]*sum(varloadP[i] - pd[i] for i=keys(netload))
    )
    
    @constraint(
        m, LimAngDown[b=keys(netbranch)],
        varva[f_bus[b]] - varva[t_bus[b]] >= angmin[b]
    )
    @constraint(
        m, LimAngUp[b=keys(netbranch)],
        varva[f_bus[b]] - varva[t_bus[b]] <= angmax[b]
    )
    @constraint(m, RefAngle[i=ref_bus], varva[i]==network["bus"][i]["va"])

    @NLconstraint(
        m, Deffp_fr[i=keys(netbranch)],
        varfp_fr[i] == (g[i]+g_fr[i])*(varvm[f_bus[i]]/tap[i])^2
        - g[i]*varvm[f_bus[i]]/tap[i]*varvm[t_bus[i]]*(
            cos(varva[f_bus[i]]-varva[t_bus[i]]-shift[i]))
        - b[i]*varvm[f_bus[i]]/tap[i]*varvm[t_bus[i]]*(
            sin(varva[f_bus[i]]-varva[t_bus[i]]-shift[i]))
    )
    @NLconstraint(
        m, Deffq_fr[i=keys(netbranch)],
        varfq_fr[i] == -(b[i]+b_fr[i])*(varvm[f_bus[i]]/tap[i])^2
        + b[i]*varvm[f_bus[i]]/tap[i]*varvm[t_bus[i]]*(
            cos(varva[f_bus[i]]-varva[t_bus[i]]-shift[i]))
        - g[i]*varvm[f_bus[i]]/tap[i]*varvm[t_bus[i]]*(
            sin(varva[f_bus[i]]-varva[t_bus[i]]-shift[i]))
    )
    @NLconstraint(
        m, Deffp_to[i=keys(netbranch)],
        varfp_to[i] == (g[i]+g_to[i])*(varvm[t_bus[i]])^2
        - g[i]*varvm[t_bus[i]]*varvm[f_bus[i]]/tap[i]*(
            cos(varva[t_bus[i]]-varva[f_bus[i]]+shift[i]))
        - b[i]*varvm[t_bus[i]]*varvm[f_bus[i]]/tap[i]*(
            sin(varva[t_bus[i]]-varva[f_bus[i]]+shift[i]))
    )
    @NLconstraint(
        m, Deffq_to[i=keys(netbranch)],
        varfq_to[i] == -(b[i]+b_to[i])*(varvm[t_bus[i]])^2
        + b[i]*varvm[t_bus[i]]*varvm[f_bus[i]]/tap[i]*(
            cos(varva[t_bus[i]]-varva[f_bus[i]]+shift[i]))
        - g[i]*varvm[t_bus[i]]*varvm[f_bus[i]]/tap[i]*(
            sin(varva[t_bus[i]]-varva[f_bus[i]]+shift[i]))
    )
    
    @constraint(
        m, RealPowerBalance[i=keys(netbus)],
        sum(varfp_fr[b] for b=branch_f[i]) + sum(varfp_to[b] for b=branch_t[i])
        == sum(varpg[g] for g=gen_i[i]) - sum(varloadP[d] for d=load_i[i])
        - sum(netshunt[s]["gs"] for s=shunt_i[i])*varvm[i]^2
    )            

    @constraint(
        m, ReactivePowerBalance[i=keys(netbus)],
        sum(varfq_fr[b] for b=branch_f[i]) + sum(varfq_to[b] for b=branch_t[i])
        == sum(varqg[g] for g=gen_i[i]) - sum(varloadQ[d] for d=load_i[i]) 
        + sum(netshunt[s]["bs"] for s=shunt_i[i])*varvm[i]^2
    )

    return m
end

function solveOPF_line(line::String,
                       PG0::Dict{String,Float64},V0::Dict{String,Float64},
                       pertLoads::Set{String},
                       threshold::Float64,pert::Float64)::Scenario
    m=OPF_constraints(PG0,V0,pertLoads,pert)
    @variable(m, t >= 0)
    var=JuMP.all_variables(m)
    varfp_fr=VariableRef(m)
    varfq_fr=VariableRef(m)
    varfp_to=VariableRef(m)
    varfq_to=VariableRef(m)
    varloadP=Dict{String,VariableRef}()
    varloadQ=Dict{String,VariableRef}()
    for v=var
        if occursin(string("varfp_fr[",line,"]"),name(v)) varfp_fr=v end
        if occursin(string("varfq_fr[",line,"]"),name(v)) varfq_fr=v end
        if occursin(string("varfp_to[",line,"]"),name(v)) varfp_to=v end
        if occursin(string("varfq_to[",line,"]"),name(v)) varfq_to=v end
        if occursin("varloadP",name(v))
            get!(varloadP,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],
                 v)
        end
        if occursin("varloadQ",name(v))
            get!(varloadQ,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],
                 v)
        end
        
    end
    
    @constraint(m, tLimFr, t <= varfp_fr^2 + varfq_fr^2)
    @constraint(m, tLimTo, t <= varfp_to^2 + varfq_to^2)
    
    @objective(m, Max, t) 
    optimize!(m, with_optimizer(Ipopt.Optimizer,print_level=0,max_iter=250,
                                tol=1e-6,dual_inf_tol=1e-6,constr_viol_tol=1e-6,
                                warm_start_init_point="yes"
                                ))
    if termination_status(m) != LOCALLY_SOLVED &&
        termination_status(m) != ALMOST_LOCALLY_SOLVED
        println(">>> Issue when solving opf on line ", line)
        return  Scenario(
            Dict{String,Float64}(
                i=> value(varloadP[i]) - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )
    end
    if sqrt(value(t)) > (1-threshold)*rate_a[line]
        return Scenario(
            Dict{String,Float64}(
                i=> value(varloadP[i]) - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )
    end
    return Scenario(Dict{String,Float64}(),Dict{String,Float64}())
end

function solveOPF_VUp(bus::String,
                      PG0::Dict{String,Float64},V0::Dict{String,Float64},
                      pertLoads::Set{String},
                      threshold::Float64,pert::Float64)::Scenario
    m=OPF_constraints(PG0,V0,pertLoads,pert)
    var=JuMP.all_variables(m)
    varvm=VariableRef(m)
    varloadP=Dict{String,VariableRef}()
    varloadQ=Dict{String,VariableRef}()
    for v=var
        if occursin(string("varvm[",bus,"]"),name(v))
            varvm=v
        end
         if occursin("varloadP",name(v))
            get!(varloadP,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
        if occursin("varloadQ",name(v))
            get!(varloadQ,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
    end
    @objective(m, Max, varvm)
    
    optimize!(m, with_optimizer(Ipopt.Optimizer,print_level=0,max_iter=250,
                                tol=1e-6,dual_inf_tol=1e-6,constr_viol_tol=1e-6,
                                warm_start_init_point="yes"
                                ))
    if termination_status(m) != LOCALLY_SOLVED &&
        termination_status(m) != ALMOST_LOCALLY_SOLVED
        println(">>> Issue when solving opf on bus ", bus)
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i])  - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )    end
    if value(varvm) > (1-threshold)*vmax[bus]
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i])  - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )
    end
    return Scenario(Dict{String,Float64}(),Dict{String,Float64}())
end


function solveOPF_VDown(bus::String,
                        PG0::Dict{String,Float64},V0::Dict{String,Float64},
                        pertLoads::Set{String},
                        threshold::Float64,pert::Float64)::Scenario
    m=OPF_constraints(PG0,V0,pertLoads,pert)
    var=JuMP.all_variables(m)
    varvm=VariableRef(m)
    varloadP=Dict{String,VariableRef}()
    varloadQ=Dict{String,VariableRef}()
    for v=var
        if occursin(string("varvm[",bus,"]"),name(v))
            varvm=v
        end
         if occursin("varloadP",name(v))
            get!(varloadP,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
        if occursin("varloadQ",name(v))
            get!(varloadQ,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
    end
    @objective(m, Min, varvm)

    optimize!(m, with_optimizer(Ipopt.Optimizer,print_level=0,max_iter=250,
                                tol=1e-6,dual_inf_tol=1e-6,constr_viol_tol=1e-6,
                                warm_start_init_point="yes"
                                ))
    if termination_status(m) != LOCALLY_SOLVED &&
        termination_status(m) != ALMOST_LOCALLY_SOLVED
        println(">>> Issue when solving opf on bus ", bus)
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i])  - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )
    end
    if value(varvm) < (1+threshold)*vmin[bus]
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i]) - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )
    end
    return Scenario(Dict{String,Float64}(),Dict{String,Float64}())
end

function solveOPF_QUp(gen::String,
                      PG0::Dict{String,Float64}, V0::Dict{String,Float64},
                      pertLoads::Set{String},
                      threshold::Float64,pert::Float64)::Scenario
    m=OPF_constraints(PG0,V0,pertLoads,pert)
    var=JuMP.all_variables(m)
    varqg=VariableRef(m)
    varloadP=Dict{String,VariableRef}()
    varloadQ=Dict{String,VariableRef}()
    for v=var
        if occursin(string("varqg[",gen,"]"),name(v))
            varqg=v
        end
         if occursin("varloadP",name(v))
            get!(varloadP,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
        if occursin("varloadQ",name(v))
            get!(varloadQ,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
    end
    @objective(m, Max, varqg)
    
    optimize!(m, with_optimizer(Ipopt.Optimizer,print_level=0,max_iter=250,
                                tol=1e-6,dual_inf_tol=1e-6,constr_viol_tol=1e-6,
                                warm_start_init_point="yes"
                                ))
    if termination_status(m) != LOCALLY_SOLVED &&
        termination_status(m) != ALMOST_LOCALLY_SOLVED
        println(">>> Issue when solving opf on gen ", gen)
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i])  - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )
    end
    if value(varqg) > (1-sign(qmax[gen])*threshold)*qmax[gen]
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i])  - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) -  qd[i] for i=keys(netload))
        )
    end
    return Scenario(Dict{String,Float64}(),Dict{String,Float64}())
end


function solveOPF_QDown(gen::String,
                        PG0::Dict{String,Float64},V0::Dict{String,Float64},
                        pertLoads::Set{String},
                        threshold::Float64,pert::Float64)::Scenario
    m=OPF_constraints(PG0,V0,pertLoads,pert)
    var=JuMP.all_variables(m)
    varqg=VariableRef(m)
    varloadP=Dict{String,VariableRef}()
    varloadQ=Dict{String,VariableRef}()
    for v=var
        if occursin(string("varqg[",gen,"]"),name(v))
            varqg=v
        end
         if occursin("varloadP",name(v))
            get!(varloadP,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
        if occursin("varloadQ",name(v))
            get!(varloadQ,
                 name(v)[findfirst(isequal('['),
                                   name(v))+1:findfirst(isequal(']'),
                                                        name(v))-1],v)
        end
    end
    @objective(m, Min, varqg)

    optimize!(m, with_optimizer(Ipopt.Optimizer,print_level=0,max_iter=250,
                                tol=1e-6,dual_inf_tol=1e-6,constr_viol_tol=1e-6,
                                warm_start_init_point="yes"
                                ))
    if termination_status(m) != LOCALLY_SOLVED &&
        termination_status(m) != ALMOST_LOCALLY_SOLVED
        println(">>> Issue when solving opf on gen", gen)
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i])  - pd[i] for i=keys(netload)),
            Dict{String,Float64}(
                i => value(varloadQ[i]) - qd[i] for i=keys(netload))
        )    end
    if value(varqg) < (1+sign(qmin[gen])*threshold)*qmin[gen]
        return Scenario(
            Dict{String,Float64}(
                i => value(varloadP[i]) - pd[i] for i=keys(netload)),
            Dict{String,Float64}(i => value(varloadQ[i]) - qd[i]
                                 for i=keys(netload))
        )
    end
    return Scenario(Dict{String,Float64}(),Dict{String,Float64}())
end


function FindCriticalScenarios(PG0::Dict{String,Float64},
                               V0::Dict{String,Float64},
                               pertLoads::Set{String},
                               threshold::Float64,pert::Float64)
    println("Finding critical scenarios...")
    criticalScen=Array{Scenario}(undef,0)
    weightScen=Array{Float64}(undef,0)
    print("From lines...")
    for i=keys(netbranch)
        scen=solveOPF_line(i,PG0,V0,pertLoads,threshold,pert)
        if length(scen)==0 continue end
        if !(scen in criticalScen)
            push!(criticalScen, scen)
            push!(weightScen,1.0)
        else
            index_s=1
            for s=criticalScen
                if s==scen continue end
                index_s+=1
            end
            weightScen[index_s]+=1
        end
    end
    println("OK.")
    print("From voltages...")
    for i=PQbus
        scen=solveOPF_VUp(i,PG0,V0,pertLoads,threshold,pert)
        if length(scen)==0 continue end
        if !(scen in criticalScen)
            push!(criticalScen, scen)
            push!(weightScen,1.0)
        else
            index_s=1
            for s=criticalScen
                if s==scen continue end
                index_s+=1
            end
            weightScen[index_s]+=1
        end
    end
    for i=PQbus
        scen=solveOPF_VDown(i,PG0,V0,pertLoads,threshold,pert)
        if length(scen)==0 continue end
        if !(scen in criticalScen)
            push!(criticalScen, scen)
            push!(weightScen,1.0)
        else
            index_s=1
            for s=criticalScen
                if s==scen continue end
                index_s+=1
            end
            weightScen[index_s]+=1
        end
    end
    println("OK.")
    print("From reactive power generation...")
    for g=NonSlackGen
        scen=solveOPF_QUp(g,PG0,V0,pertLoads,threshold,pert)
        if length(scen)==0 continue end
        if !(scen in criticalScen)
            push!(criticalScen, scen)
            push!(weightScen,1.0)
        else
            index_s=1
            for s=criticalScen
                if s==scen continue end
                index_s+=1
            end
            weightScen[index_s]+=1
        end
    end
    for g=NonSlackGen
        scen=solveOPF_QDown(g,PG0,V0,pertLoads,threshold,pert)
        if length(scen)==0 continue end
        if !(scen in criticalScen)
            push!(criticalScen, scen)
            push!(weightScen,1.0)
        else
            index_s=1
            for s=criticalScen
                if s==scen continue end
                index_s+=1
            end
            weightScen[index_s]+=1
        end
    end
    println("OK.")
    println("#################################")
    println("### Found ",length(criticalScen)
            , " / ", sum(weightScen), " / ",
            length(netbranch)+2*(length(PQbus)+length(NonSlackGen)),
            " critical scenarios ###")
    println("#################################")
    return criticalScen,weightScen
end


function checkCriticalScen(CriticalScen::Array{Scenario})
    print("Checking critical scenarios...")
    copyCriticalScen=Array{Scenario}(CriticalScen)
    count=1
    for scen=copyCriticalScen
        for i=keys(network["load"])
            network["load"][i]["pd"]=pd[i]+scen.loadP[i]
            network["load"][i]["qd"]=qd[i]+scen.loadQ[i]
        end
        res=run_ac_opf(network,with_optimizer(Ipopt.Optimizer,
                                              print_level=0,
                                              max_iter=250))
        if res["termination_status"] != LOCALLY_SOLVED
            deleteat!(CriticalScen,count)
            count-=1
        end
        count+=1
    end
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    println("OK. ", length(CriticalScen), " critical scenarios remaining.") 
end


function ReduceCriticalScen(CriticalScen::Array{Scenario},
                            nbScen::Int64
                            )::(Array{Scenario},Float64)
    arrayScen=Array{Float64}(undef, 2*length(netload),length(CriticalScen))
    for i=1:length(CriticalScen)
        arrayScen[:,i]=vcat(
            [values(CriticalScen[i].loadP[j]) for j=keys(netload)],
            [values(CriticalScen[i].loadQ[j]) for j=keys(netload)]
        )
    end
    R=kmeans(arrayScen,nbScen;maxiter=200)
    M=R.centers
    MP=sum(M[i,1:nbScen] for i=1:length(netload))
    MQ=sum(M[i,1:nbScen] for i=length(netload)+1:2*length(netload))
    
    arrayScenP=sum(arrayScen[i,1:length(CriticalScen)]
                   for i=1:length(netload))
    arrayScenQ=sum(arrayScen[i,1:length(CriticalScen)]
                   for i=length(netload)+1:2*length(netload))  
    NewCriticalScen=Array{Scenario}(undef,nbScen)
    for i=1:nbScen
        NewCriticalScen[i]=Scenario(
            Dict{String,Float64}(j => 0.0 for j=keys(netload)),
            Dict{String,Float64}(i => 0.0 for i=keys(netload))
        )
        k=1
        for j=keys(netload)
            NewCriticalScen[i].loadP[j]=M[k,i]
            NewCriticalScen[i].loadQ[j]=M[k+length(netload),i]
            k+=1
        end    
    end
    return NewCriticalScen, R.totalcost
end

function ReduceCriticalScen(CriticalScen::Array{Scenario},
                            WeightScen::Array{Float64},
                            nbScen::Int64)::(Array{Scenario},Array{Scenario},
                                             Float64)
    arrayScen=Array{Float64}(undef, 2*length(netload),length(CriticalScen))
    for i=1:length(CriticalScen)
        arrayScen[:,i]=vcat(
            [values(CriticalScen[i].loadP[j]) for j=keys(netload)],
            [values(CriticalScen[i].loadQ[j]) for j=keys(netload)]
        )
    end
    R=kmeans(arrayScen,nbScen;maxiter=200)
    Mbis=Array{Float64}(undef, 2*length(netload),nbScen)
    for s=1:nbScen
        a=assignments(R)
        index_s=collect(i for i=1:length(CriticalScen) if a[i]==s)
        indexmax=0
        maxi=0.0
        for i=index_s
            if WeightScen[i] > maxi
                maxi=WeightScen[i]
                indexmax=i
            end
        end
        Mbis[:,s]=arrayScen[:,indexmax]
    end
    M=R.centers
    MbisP=sum(Mbis[i,1:nbScen] for i=1:length(netload))
    MbisQ=sum(Mbis[i,1:nbScen] for i=length(netload)+1:2*length(netload))
    MP=sum(M[i,1:nbScen] for i=1:length(netload))
    MQ=sum(M[i,1:nbScen] for i=length(netload)+1:2*length(netload))
    arrayScenP=sum(arrayScen[i,1:length(CriticalScen)]
                   for i=1:length(netload))
    arrayScenQ=sum(arrayScen[i,1:length(CriticalScen)]
                   for i=length(netload)+1:2*length(netload))
    NewCriticalScen=Array{Scenario}(undef,nbScen)
    for i=1:nbScen
        NewCriticalScen[i]=Scenario(
            Dict{String,Float64}(j => 0.0 for j=keys(netload)),
            Dict{String,Float64}(j => 0.0 for j=keys(netload))
        )
        k=1
        for j=keys(netload)
            NewCriticalScen[i].loadP[j]=Mbis[k,i]
            NewCriticalScen[i].loadQ[j]=Mbis[k+length(netload),i]
            k+=1
        end    
    end
    CentersScen=Array{Scenario}(undef,nbScen)
    for i=1:nbScen
        CentersScen[i]=Scenario(
            Dict{String,Float64}(j => 0.0 for j=keys(netload)),
            Dict{String,Float64}(j => 0.0 for j=keys(netload))
        )
        k=1
        for j=keys(netload)
            CentersScen[i].loadP[j]=M[k,i]
            CentersScen[i].loadQ[j]=M[k+length(netload),i]
            k+=1
        end    
    end
    return NewCriticalScen,CentersScen, R.totalcost
end


function ComputeWeightScen(Scen::Array{Scenario},
                           PG0::Dict{String,Float64},V0::Dict{String,Float64},
                           method::String)::Array{Float64}
    weight=Array{Float64}(undef,length(Scen))
    for i=keys(netgen)
        network["gen"][i]["vg"]=V0[string(netgen[i]["gen_bus"])]
    end
    for i=keys(netbus)
        if netbus[i]["bus_type"]!=1
            network["bus"][i]["vm"]=V0[i]
        end
    end
    for s=1:length(Scen)
        for i=pertLoads
            network["load"][i]["pd"]=pd[i] + Scen[s].loadP[i]
            network["load"][i]["qd"]=qd[i] + Scen[s].loadQ[i]
        end
        w=sum(Scen[s].loadP[i] for i=keys(netload))
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
        ConstrViol,ViolConstr=checkSolution(results)
        if method=="NbConstr"
            weight[s]=Float64(length(ConstrViol))
        else
            weight[s]=maximum(ViolConstr)
        end
    end
    
    for i=keys(network["load"])
        network["load"][i]["pd"]=pd[i]
        network["load"][i]["qd"]=qd[i]
    end
    return weight
end
