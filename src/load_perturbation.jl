
function DetectEndLoads()::Set{String}
    EndLoads=Set{String}()
    for i=keys(netload)
        curbus=string(netload[i]["load_bus"])
        if length(union(branch_f[curbus],branch_t[curbus]))==1
            push!(EndLoads,i)
        end
    end
    ### Particular case of case1354pegase
    if "513" in EndLoads delete!(EndLoads,"513") end
    return EndLoads
end

function BoundsLoads(pertLoads::Set{String},pert::Float64)
    pertPUp=Dict{String,Float64}(pd)
    pertPDown=Dict{String,Float64}(pd)
    pertQUp=Dict{String,Float64}(qd)
    pertQDown=Dict{String,Float64}(qd)
    for i=pertLoads
        pertPUp[i]+=pert*abs(pd[i])
        pertPDown[i]-=pert*abs(pd[i])
        pertQUp[i]+=pert*abs(qd[i])
        pertQDown[i]-=pert*abs(qd[i])
    end
    return (pertPUp,pertPDown,pertQUp,pertQDown)
end

function DefAlpha(pertLoads::Set{String},pert::Float64)::Dict{String,Float64}
    (pertPUp,pertPDown,pertQUp,pertQDown) = BoundsLoads(pertLoads,pert)
    maxPertP=sum(values(pertPUp))-sum(values(pd))
    maxPertQ=sum(values(pertQUp))-sum(values(qd))
    alpha=Dict{String,Float64}(i=>0.0 for i=keys(netgen))
    initAlpha=length(netgen)/10
    normAlpha=0.0
    for i=keys(alpha)
        alpha[i]=1.0+1.0*sign(min(
            (pmax[i]-pmin[i])*initAlpha/4 - maxPertP,
            (qmax[i]-qmin[i])*initAlpha/4 - maxPertQ
        ))
    end
    while normAlpha!=sum(values(alpha))
        normAlpha=sum(values(alpha))
        for i=keys(alpha)
            alpha[i]=1.0+1.0*sign(min(
                (pmax[i]-pmin[i])*initAlpha/4 - maxPertP,
                (qmax[i]-qmin[i])*initAlpha/4 - maxPertQ
            ))
        end
    end
    for i=keys(alpha)
        alpha[i]/=normAlpha
        if pmax[i]-pmin[i] != 0.0
            @assert alpha[i]*maxPertP <= (pmax[i]-pmin[i])/2
        end
        if qmax[i]-qmin[i] != 0.0
            @assert alpha[i]*maxPertQ <= (qmax[i]-qmin[i])/2
        end
    end
    return alpha
end
