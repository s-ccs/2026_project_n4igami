using Unfold
using WGLMakie, UnfoldMakie, MakieThemes, CairoMakie # plotting
using Revise
using TidierData
using Missings
using Statistics
using Plots
using DelimitedFiles
using DataFrames

includet("../src/io_eeglab.jl")
datasetFolder = "/Volumes/Z Slim/Eyesorted"
datasetFiles = filter(f -> endswith(f, ".set") && !startswith(basename(f), "._"), readdir(abspath(datasetFolder), join=true))
numDatasets = length(datasetFiles)
#setName_list = [SubString(basename(f), 1, 6) for f in datasetFiles]  # subject ID is the 1st to 6th characters of the filename
save_path_dc = raw"/Users/saramilligan/N4igami/QL_model_output/deconv"
mkpath(save_path_dc)
save_path_epoch = raw"/Users/saramilligan/N4igami/QL_model_output/epoch"
mkpath(save_path_epoch)
artrej_summary_folder = raw"/Users/saramilligan/N4igami/QL_model_output/artifact_rejection_summary"
mkpath(artrej_summary_folder)
save_path_artrej_frp = raw"/Users/saramilligan/N4igami/QL_model_output/artifact_rejection_summary/artrej_FRP_plots"
mkpath(save_path_artrej_frp)
save_path_frp = raw"/Users/saramilligan/N4igami/QL_model_output/FRP_plots"
mkpath(save_path_frp)

for i in 1:numDatasets
    setName = SubString(basename(datasetFiles[i]), 1, 6)  # subject ID is the 1st to 6th characters of the filename
    artrej_filename = joinpath(artrej_summary_folder, setName * "_artrej_summary.csv")
    dc_model_filename = joinpath(save_path_dc, setName * "_dc_model.jld2")
    epoch_model_filename = joinpath(save_path_epoch, setName * "_epoch_model.jld2")
    artrej_plot_filename = joinpath(save_path_artrej_frp, setName * "_artrej_plot.png")
    frp_plot_filename = joinpath(save_path_frp, setName * "_FRP_plot.png")
    all_channel_names = ["Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "FT9", "FC5", "FC1", "FC2", "FC6", "FT10", "T7", "C3", "Cz", "C4", "T8", "CP5", "CP1", "CP2", "CP6", "P7", "P3", "Pz", "P4", "P8", "O1", "Oz", "O2", "biPolHEOG", "R-GAZE-X"]
    #below are all channels in analysis ROIs (central/parietal/occipital channels)
    ch_ix_ROI = [15, 16, 17, 20, 21, 24, 25, 26, 28, 29, 30]
    #channel names in ROI [C3, Cz, C4, CP1, CP2, P3, Pz, P4, O1, Oz, O2]
    #all channel names [Fp1, Fp2, F7, F3, Fz, F4, F8, FT9, FC5, FC1, FC2, FC6, FT10, T7, C3, Cz, C4, T8, CP5, CP1, CP2, CP6, P7, P3, Pz, P4, P8, O1, Oz, O2]
    #channel name list: 1 = Fp1, 2 = Fp2, 3 = F7, 4 = F3, 5 = Fz, 6 = F4, 7 = F8, 8 = FT9, 9 = FC5, 10 = FC1, 11 = FC2, 12 = FC613 = FT10, 14 = T7, 15 = C3, 16 = Cz, 17 = C4, 18 = T8, 19 = CP5, 20 = CP1, 21 = CP2, 22 = CP6, 23 = P7,24 = P3, 25 = Pz, 26 = P4, 27 = P8, 28 = O1, 29 = Oz, 30 = O2, 31 = biPolHEOG, 32 = R-GAZE-X
    #frontal: 4, 5, 6, 10, 11
    channel_labels = Dict(1 => "Fp1", 2 => "Fp2", 3 => "F7", 4 => "F3", 5 => "Fz", 6 => "F4", 7 => "F8", 8 => "FT9", 9 => "FC5", 10 => "FC1", 11 => "FC2", 12 => "FC6", 13 => "FT10", 14 => "T7", 15 => "C3", 16 => "Cz", 17 => "C4", 18 => "T8", 19 => "CP5", 20 => "CP1", 21 => "CP2", 22 => "CP6", 23 => "P7", 24 => "P3", 25 => "Pz", 26 => "P4", 27 => "P8", 28 => "O1", 29 => "Oz", 30 => "O2", 31 => "biPolHEOG", 32 => "R-GAZE-X")
    channel_labels_ROI = Dict(1 => "C3", 2 => "Cz", 3 => "C4", 4 => "CP1", 5 => "CP2", 6 => "P3", 7 => "Pz", 8 => "P4", 9 => "O1", 10 => "Oz", 11 => "O2")

    println("Processing dataset $i of $numDatasets: $setName")  
    data, evts, eeg_dict = import_eeglab(datasetFiles[i])
    evts_clean = DataFrame(@chain evts begin
        @mutate(Sac_amp_in = lag(sac_amplitude,1),
                Sac_amp_out = lead(sac_amplitude,1))
        @filter(current_region in ("Beginning", "PreTarget", "Target_word", "Ending") || type in ("S250", "S251"))
        @mutate(event_type = ifelse.(type in ("S250", "S251"), type, current_region))
        @filter(region1_start != 0)
        @filter(region_pass_number == 1 || type in ("S250", "S251"))
        @mutate(event_type = ifelse.(event_type == "PreTarget", bdf_label_description, event_type))
        @mutate(
            Condition = case_when(
                condition_number in [211, 221, 231, 241] => "EX",
                condition_number in [212, 222, 232, 242] => "ON",
                condition_number in [213, 223, 233, 243] => "NN",
                true  => ""
            ))
    end)

    #artifact detection and exclusion
    mask = detectbad_peak_to_peak(data; threshold = 100, sfreq = 500, window = 0.5)
    #check proportion of bad data per channel
    col_mask = mean(mask, dims=2)

    #replace bad segments with missing values
    data_clean = allowmissing(data)
    data_clean[mask] .= missing

    #make sure bad segments were replced with missing values - proportion values should be the same as proportions in mask
    local col_missing = mean(ismissing.(data_clean), dims=2)
    local missing_by_chan = hcat(all_channel_names, col_missing)
    
    println("Saving artifact rejection summary for $setName to $artrej_filename")
    writedlm(artrej_filename, missing_by_chan, ',')
    
    #CURRENTLY USING ONLY ONE CHANNEL (Cz) FOR MODELING, BUT CAN CHANGE TO MULTIPLE CHANNELS IF DESIRED
    #ch_ix = 16

    #run intercept only models to plot/check for differences in the data with and without bad data
    
    m =
        fit(UnfoldModel, @formula(0 ~ 1), evts_clean, data[ch_ix_ROI,:], firbasis((-0.2, 1.0), 500))
    c = coeftable(m)
    m_missing =
        fit(UnfoldModel, @formula(0 ~ 1), evts_clean, data_clean[ch_ix_ROI,:], firbasis((-0.2, 1.0), 500))
    c.group .= "with bad data"
    c_missing = coeftable(m_missing)
    c_missing.group .= "bad data removed"
    #plot modeled data with and without bad data
    eff_all_missing = vcat(c,c_missing)
    eff_all_missing.Channel_Name = get.(Ref(channel_labels_ROI), eff_all_missing.channel, "Unknown") #add channel names to the effect coded estimates for plotting

    println("Saving waveform plots showing effect of artifact rejection at Cz, Pz, Oz for $i of $numDatasets: $setName")
    artrej_FRP_plot = plot_erp(subset(eff_all_missing, :Channel_Name => ByRow(x -> x .∈ Ref(["Cz", "Pz", "Oz"]))), mapping = (; color = :group, row = :Channel_Name))
    #artrej_FRP_plot = plot_erp(eff_all_missing, mapping = (; color = :group, row = :channel))

    save(artrej_plot_filename, artrej_FRP_plot)

    #evts_clean[evts_clean.fix_type .=="Beginning","type"] .= "other_fixation"
    #deconvolution model design
    design = [
        "S250" => (@formula(0 ~ 1),firbasis([-0.2,1.5],500)),
        "S251" => (@formula(0 ~ 1),firbasis([-0.5,0.5],500)),
        "Beginning" => (@formula(0 ~ 1),firbasis([-0.2,1],500)),
        "Pretarget Skip" => (@formula(0 ~ 1 + Condition),firbasis([-0.2,1],500)),
        "Pretarget NoSkip" => (@formula(0 ~ 1 + Condition),firbasis([-0.2,1],500)),
        "Target_word" => (@formula(0 ~ 1 + Condition),firbasis([-0.2,1],500)),
        "Ending" => (@formula(0 ~ 1 + Condition),firbasis([-0.2,1],500))
    ]
    #model design for epoched data without deconvolution to compare with deconvolution results
    data_epoch,times=Unfold.epoch(;data=data_clean,evts=evts_clean,τ=(-0.2,1),sfreq=500)
    design_epoch = [
        "S250" => (@formula(0 ~ 1),times),
        "S251" => (@formula(0 ~ 1),times),
        "Beginning" => (@formula(0 ~ 1),times),
        "Pretarget Skip" => (@formula(0 ~ 1 + Condition),times),
        "Pretarget NoSkip" => (@formula(0 ~ 1 + Condition),times),
        "Target_word" => (@formula(0 ~ 1 + Condition),times),
        "Ending" => (@formula(0 ~ 1 + Condition),times)
    ]

    #actually running deconvolution model with bad segments removed
    println("Unfolding dataset $i of $numDatasets: $setName")
    m_dc = fit(UnfoldModel, design,evts_clean,data_clean[ch_ix_ROI,:];eventcolumn=:event_type)
    c_dc = coef(m_dc)
    co_dc = coeftable(m_dc)

    #running model on epoched data without deconvolution to compare with deconvolution results
    println("Modeling epoched dataset $i of $numDatasets: $setName")
    m_epoch = fit(UnfoldModel,design_epoch,evts_clean,data_epoch[ch_ix_ROI,:,:];eventcolumn=:event_type)
    c_epoch = coef(m_epoch)
    co_epoch = coeftable(m_epoch)

    println("Converting deconvolution model output to effect coded estimates for $i of $numDatasets: $setName")
    _eff_dc = effects(Dict(:Condition =>["EX", "ON", "NN"]),m_dc)
    eff_dc_4plot = subset(_eff_dc, :eventname => ByRow(x -> x ∉ ["S250", "S251", "Beginning"])) #subset to only include events of interest for plotting
    eff_dc_4plot.Channel_Name = get.(Ref(channel_labels_ROI), eff_dc_4plot.channel, "Unknown")    #add channel names to the effect coded estimates for plotting

    println("Converting epoched model output to effect coded estimates for $i of $numDatasets: $setName")
    _eff_epoch = effects(Dict(:Condition =>["EX", "ON", "NN"]),m_epoch)
    eff_epoch_4plot = subset(_eff_epoch,:eventname =>  ByRow(x -> x ∉ ["S250", "S251", "Beginning"])) #subset to only include events of interest for plotting
    eff_epoch_4plot.Channel_Name = get.(Ref(channel_labels_ROI), eff_epoch_4plot.channel, "Unknown") #add channel names to the effect coded estimates for plotting

    #combining deconvolution and epoched model output for plotting
    eff_dc_4plot.group .= "deconv"
    eff_epoch_4plot.group .= "epoch"
    eff_all = vcat(eff_dc_4plot,eff_epoch_4plot)
    eff_all.subject .= setName

    println("Saving waveform plots for $i of $numDatasets: $setName")
    
    FRP_plot_Oz = plot_erp(subset(eff_all, :Channel_Name => ByRow(x -> x == "Oz")), mapping=(;color=:Condition,col=:eventname,row=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])))
    FRP_plot_Cz = plot_erp(subset(eff_all, :Channel_Name => ByRow(x -> x == "Cz")), mapping=(;color=:Condition,col=:eventname,row=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])))
    
    FRP_plot = Figure(size = (1000, 800))
    with_theme(theme_ggthemr(:fresh)) do
        # Left Plot: Cz Channel (Legend Disabled safely via layout NamedTuple)
        plot_erp!(
            FRP_plot[1, 1], 
            subset(eff_all, :Channel_Name => ByRow(x -> x == "Cz")); 
            mapping = (; 
                color = :Condition, 
                row = :eventname, 
                col = :group => UnfoldMakie.AlgebraOfGraphics.sorter(["epoch", "deconv"]) 
            ), 
            axis = (; title = "Cz"),
            layout = (; use_legend = true, titleposition = :left),
            legend = (; position = :bottom, orientation = :horizontal)
        )

        # Right Plot: Oz Channel (Legend configuration)
        plot_erp!(
            FRP_plot[1, 2], 
            subset(eff_all, :Channel_Name => ByRow(x -> x == "Oz")); 
            mapping = (; 
                color = :Condition, 
                row = :eventname, 
                col = :group => UnfoldMakie.AlgebraOfGraphics.sorter(["epoch", "deconv"])
            ), 
            axis = (; title = "Oz"),
            layout = (; use_legend = true, titleposition = :left), # Keeps the legend enabled
            legend = (; position = :bottom, orientation = :horizontal) # Pass Makie.Legend attributes as a NamedTuple
        )
    end
    save(frp_plot_filename, FRP_plot)

    println("Saving $setName deconvolution model results")
    save(dc_model_filename, m_dc; compress = true);

    println("Saving $setName epoch model results")
    save(epoch_model_filename, m_epoch; compress = true);
    
end
#plot_erp(co)   
#plot_erp(co2)
#plot_erp(co2,mapping=(;col=:eventname),legend=(;position=:bottom,orientation=:horizontal))




# ERPExplorer.jl - discover(m)
# TidierPlots.jl (ggplot)
# TidierData.jl (or Tidier.jl for everything)

