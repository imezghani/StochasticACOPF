struct Scenario
    loadP::Dict{String,Float64}
    loadQ::Dict{String,Float64}
end

Base.:(==)(x::Scenario, y::Scenario) = (
    if (length(keys(x.loadP))!=length(keys(y.loadP))
        || length(keys(x.loadQ))!=length(keys(y.loadQ)))
    return false
    else
    for i=keys(x.loadP)
    if abs(x.loadP[i]-y.loadP[i]) > tol return false end
    if abs(x.loadQ[i]-y.loadQ[i]) > tol return false end
    end
    return true
    end
)

Base.:in(x::Scenario, y::Set{Scenario}) = if isempty(y) return false
else for i=y if x==i return true; end; end; return false; end;

Base.length(x::Scenario) = return length(x.loadP)+length(x.loadQ)

