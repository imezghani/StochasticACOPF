netgen=network["gen"]
netbus=network["bus"]
netbranch=network["branch"]
netload=network["load"]
netshunt=network["shunt"]
PQbus=Set{String}([i for i=keys(netbus) if netbus[i]["bus_type"]==1])
PVbus=Set{String}([i for i=keys(netbus) if netbus[i]["bus_type"]==2])
f_bus=Dict{String, String}(
    b=> string(netbranch[b]["f_bus"]) for b=keys(netbranch)
)
t_bus=Dict{String, String}(
    b=> string(netbranch[b]["t_bus"]) for b=keys(netbranch)
)
branch_f=Dict{String, Set{String}}(
    i => Set{String}(b for b=keys(netbranch) if f_bus[b]==i) for i=keys(netbus)
)
branch_t=Dict{String, Set{String}}(
    i => Set{String}(b for b=keys(netbranch) if t_bus[b]==i) for i=keys(netbus)
)
gen_i=Dict{String, Set{String}}(
    i => Set{String}(g for g=keys(netgen) if string(netgen[g]["gen_bus"])==i)
    for i=keys(netbus)
)
load_i=Dict{String, Set{String}}(
    i => Set{String}(d for d=keys(netload) if string(netload[d]["load_bus"])==i)
    for i=keys(netbus)
)
shunt_i=Dict{String, Set{String}}(
    i => Set{String}(s for s=keys(netshunt)
                     if string(netshunt[s]["shunt_bus"])==i)
    for i=keys(netbus)
)
vmax=Dict{String,Float64}(i=> netbus[i]["vmax"] for i=keys(netbus))
vmin=Dict{String,Float64}(i=> netbus[i]["vmin"] for i=keys(netbus))
pmax=Dict{String,Float64}(i=> netgen[i]["pmax"] for i=keys(netgen))
pmin=Dict{String,Float64}(i=> netgen[i]["pmin"] for i=keys(netgen))
qmax=Dict{String,Float64}(i=> netgen[i]["qmax"] for i=keys(netgen))
qmin=Dict{String,Float64}(i=> netgen[i]["qmin"] for i=keys(netgen))
angmax=Dict{String,Float64}(i=> netbranch[i]["angmax"] for i=keys(netbranch))
angmin=Dict{String,Float64}(i=> netbranch[i]["angmin"] for i=keys(netbranch))
rate_a=Dict{String,Float64}(i=> netbranch[i]["rate_a"] for i=keys(netbranch))
ref_bus=Set{String}(i for i=keys(netbus) if netbus[i]["bus_type"]==3)
g=Dict{String,Float64}(
    i=> netbranch[i]["br_r"]/(netbranch[i]["br_r"]^2+netbranch[i]["br_x"]^2)
    for i=keys(netbranch)
)
b=Dict{String,Float64}(
    i=> -netbranch[i]["br_x"]/(netbranch[i]["br_r"]^2+netbranch[i]["br_x"]^2)
    for i=keys(netbranch)
)
g_to=Dict{String,Float64}(i=> netbranch[i]["g_to"] for i=keys(netbranch))
g_fr=Dict{String,Float64}(i=> netbranch[i]["g_fr"] for i=keys(netbranch))
b_to=Dict{String,Float64}(i=> netbranch[i]["b_to"] for i=keys(netbranch))
b_fr=Dict{String,Float64}(i=> netbranch[i]["b_fr"] for i=keys(netbranch))
tap=Dict{String,Float64}(i=> netbranch[i]["tap"] for i=keys(netbranch))
shift=Dict{String,Float64}(i=> netbranch[i]["shift"] for i=keys(netbranch))
NonSlackGen=[g for g=keys(netgen)
             if netbus[string(netgen[g]["gen_bus"])]["bus_type"]!=3]
RefGen=setdiff(keys(netgen), NonSlackGen)
pmin0=Dict{String,Float64}(pmin)
pmax0=Dict{String,Float64}(pmax)
for i=RefGen
    pmin0[i]=-Inf
    pmax0[i]=+Inf
end
vmin0=Dict{String,Float64}(vmin)
vmax0=Dict{String,Float64}(vmax)
for i=keys(netbus)
    if netbus[i]["bus_type"]==1
        vmin0[i]-=.1*vmax[i]
        vmax0[i]+=.1*vmax[i]
    end
end
qmin0=Dict{String,Float64}(qmin)
qmax0=Dict{String,Float64}(qmax)
for i=NonSlackGen
    qmin0[i]-=.5*(qmax[i]-qmin[i])
    qmax0[i]+=.5*(qmax[i]-qmin[i])
end
for i=RefGen
    qmin0[i]=-Inf
    qmax0[i]=+Inf
end
pd=Dict{String,Float64}(i=>network["load"][i]["pd"] for i=keys(netload))
qd=Dict{String,Float64}(i=>network["load"][i]["qd"] for i=keys(netload))
