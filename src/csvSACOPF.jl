
function ScenariosToCSV(Scen::Array{Scenario},fName::String)::String
    df=DataFrame(
        vcat([[Scen[i].loadP[j] for j=keys(Scen[i].loadP)] for i=keys(Scen)],
             [[Scen[i].loadQ[j] for j=keys(Scen[i].loadQ)] for i=keys(Scen)]))
    names!(df,vcat([Symbol("P",i) for i=1:length(Scen)],
                   [Symbol("Q",i) for i=1:length(Scen)]))
    insertcols!(df,1,LOAD=collect(keys(Scen[1].loadP)))
    CSV.write(string(fName,".csv"),df)
    return fName
end

function WeightToCSV(Weight::Array{Float64},fName::String)::String
    df=DataFrame(ID=1:length(Weight), Weight=Weight)
    CSV.write(string(fName,".csv"),df)
    return fName
end

function CSVToScenarios(fName::String)::Array{Scenario}
    df=CSV.read(fName)
    Scen=Array{Scenario}(undef,Int((size(df,2)-1)/2))
    for i=1:length(Scen)
        Scen[i]=Scenario(
            Dict{String,Float64}(string(df.LOAD[j]) => df[!,Symbol("P",i)][j]
                                 for j=1:size(df,1)),
            Dict{String,Float64}(string(df.LOAD[j]) => df[!,Symbol("Q",i)][j]
                                 for j=1:size(df,1))
        )
    end
    return Scen
end

function CSVToWeight(fName::String)::Array{Float64}
    df=CSV.read(fName)
    return Array{Float64}(df.Weight)
end
