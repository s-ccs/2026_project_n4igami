using Unfold
using WGLMakie, UnfoldMakie, MakieThemes, CairoMakie # plotting
using Revise
using TidierData
using Missings
using Statistics
using Plots
using DelimitedFiles
using DataFrames
using JLD2

includet("../src/io_eeglab.jl")

dc_model_path = raw"/Users/saramilligan/N4igami/QL_model_output/deconv"
epoch_model_path = raw"/Users/saramilligan/N4igami/QL_model_output/epoch"
effects_filename = raw"/Users/saramilligan/N4igami/QL_model_output/effects/QL_all_subjects_effects.csv"
save_path_effects = raw"/Users/saramilligan/N4igami/QL_model_output/effects"
mkpath(save_path_dc_effects)
dc_ModelFiles = filter(f -> endswith(f, ".jld2"), readdir(abspath(dc_model_path), join=true))
num_dcModelFiles = length(dc_ModelFiles)
epoch_ModelFiles = filter(f -> endswith(f, ".jld2"), readdir(abspath(epoch_model_path), join=true))
num_epochModelFiles = length(epoch_ModelFiles)
channel_labels_ROI = Dict(1 => "C3", 2 => "Cz", 3 => "C4", 4 => "CP1", 5 => "CP2", 6 => "P3", 7 => "Pz", 8 => "P4", 9 => "O1", 10 => "Oz", 11 => "O2")
eff_list = DataFrame[]
for i in 1:num_dcModelFiles
    println("Read in $num_dcModelFiles datasets in $dc_model_path")
    modelName = SubString(basename(dc_ModelFiles[i]), 1, 6)  # subject ID is the 1st to 6th characters of the filename
    dc_model_effects = joinpath(save_path_dc_effects, modelName * "_dc_effects.jld2")
    effects_filename = joinpath(save_path_dc_effects, modelName * "_effects.csv")
   
    m_dc = Unfold.load(dc_ModelFiles[i], UnfoldModel)
    c_dc = Unfold.coef(m_dc)
    co_dc = Unfold.coeftable(m_dc)
    
    println("Converting deconvolution model output to effect coded estimates for $i of $numDatasets: $modelName")
    _eff_dc = effects(Dict(:Condition =>["EX", "ON", "NN"]),m_dc)
    eff_dc_4plot = subset(_eff_dc, :eventname => ByRow(x -> x ∉ ["S250", "S251", "Beginning"])) #subset to only include events of interest for plotting
    eff_dc_4plot.Channel_Name = get.(Ref(channel_labels_ROI), eff_dc_4plot.channel, "Unknown")    #add channel names to the effect coded estimates for plotting
    eff_dc_4plot.group .= "deconv"
    eff_dc_4plot.subject .= modelName

    println("Read in $num_epochModelFiles datasets in $epoch_model_path")
    #modelName = SubString(basename(epoch_model_path[i]), 1, 9)  # subject ID is the 1st to 6th characters of the filename
    #epoch_model_effects = joinpath(save_path_epoch_effects, modelName * "_epoch_effects.jld2")
    m_epoch = Unfold.load(epoch_ModelFiles[i], UnfoldModel)
    c_epoch = Unfold.coef(m_epoch)
    co_epoch = Unfold.coeftable(m_epoch)

    println("Converting epoched model output to effect coded estimates for $i of $numDatasets: $modelName")
    _eff_epoch = effects(Dict(:Condition =>["EX", "ON", "NN"]),m_epoch)
    eff_epoch_4plot = subset(_eff_epoch,:eventname =>  ByRow(x -> x ∉ ["S250", "S251", "Beginning"])) #subset to only include events of interest for plotting
    eff_epoch_4plot.Channel_Name = get.(Ref(channel_labels_ROI), eff_epoch_4plot.channel, "Unknown") #add channel names to the effect coded estimates for plotting
    eff_epoch_4plot.group .= "epoch"
    eff_epoch_4plot.subject .= modelName

    eff_all = vcat(eff_dc_4plot,eff_epoch_4plot)

    println("Adding effect coded estimates for $i of $numDatasets: $modelName to effects all subjects dataframe")
    push!(eff_list, eff_all)
end

println("Saving combined effects coded estimates for all subjects to CSV file: $effects_filename")
all_effects = vcat(eff_list...)
using CSV
CSV.write(effects_filename, all_effects)

all_effects_GA = @chain all_effects begin
    @group_by(Channel_Name, eventname, group, Condition,time)
    @summarize(
        estimate = mean(yhat)
    )
end

FRP_plot_Oz = plot_erp(subset(all_effects_GA, :Channel_Name => ByRow(x -> x == "Oz")), 
    mapping=(;color=:Condition,row=:eventname,col=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])),
    axis = (; yreversed = true),
    #legend = (; position = :bottom, orientation = :horizontal)
    )
FRP_plot_Cz = plot_erp(subset(all_effects_GA, :Channel_Name => ByRow(x -> x == "Cz")),  
    mapping=(;color=:Condition,row=:eventname,col=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])),
    axis = (; yreversed = true),
    #legend = (; position = :bottom, orientation = :horizontal)
    )