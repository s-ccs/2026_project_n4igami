using Unfold
using WGLMakie, UnfoldMakie # plotting
using Revise
using JLD2
using FileIO

includet("../src/io_eeglab.jl")
#datasetFolder = "/Users/saramilligan/WTS/DC_modeling/output_deconvolution/high_constraint_WTS/continuous_processed_datasets"
datasetFolder = "/Volumes/cas/PSY Schotter Lab/Experiments/WhatTheSpell/WTS_deconv/output_deconvolution/high_constraint_WTS/UF_modified_datasets_WTS/"
datasetFiles = filter(contains(r".*\.set"),readdir(abspath(datasetFolder), join=true))
numDatasets = length(datasetFiles)
setName_list = [SubString(basename(f), 4, 9) for f in datasetFiles]  # subject ID is the 4th to 9th characters of the filename
save_path_dc = raw"/Users/saramilligan/N4igami/WTS_model_output/deconv"
save_path_epoch = raw"/Users/saramilligan/N4igami/WTS_model_output/epoch"

for i in 1:numDatasets
    setName = SubString(basename(datasetFiles[i]), 4, 9)  # subject ID is the 4th to 9th characters of the filename
    println("Processing dataset $i of $numDatasets: $setName")
    data, evts, eeg_dict = import_eeglab(datasetFiles[i])
    # Further processing can be done here for each dataset
    eff_save_name = outputFolder * setName * raw"_effects.jld2"  # name of output effects file

    remove_onset = []
    onset_active = false
    for (idx,e) in reverse(collect(enumerate(evts.fix_type)))
        if onset_active && e == "onset"
            push!(remove_onset,idx)
        elseif e == "onset" 
                onset_active = true
        else
            onset_active = false
        end
    end
    evts_clean = evts[Not(remove_onset),:]

    #evts_clean[evts_clean.fix_type .=="Beginning","type"] .= "other_fixation"
    design = [
        "onset" => (@formula(0 ~ 1),firbasis([-0.2,1.5],500)),
        "Beginning" => (@formula(0 ~ 1),firbasis([-0.2,1],500)),
        "Pretarget" => (@formula(0 ~ 1 + PT_plausibility*T_plausibility),firbasis([-0.2,1],500)),
        "Target" => (@formula(0 ~ 1 + PT_plausibility*T_plausibility),firbasis([-0.2,1],500)),
        "Ending" => (@formula(0 ~ 1 + PT_plausibility*T_plausibility),firbasis([-0.2,1],500))
    ]

    data_epoch,times=Unfold.epoch(;data=data,evts=evts_clean,τ=(-0.2,1),sfreq=500)
    design_epoch = [
        "onset" => (@formula(0 ~ 1),times),
        "Beginning" => (@formula(0 ~ 1),times),
        "Pretarget" => (@formula(0 ~ 1 + PT_plausibility*T_plausibility),times),
        "Target" => (@formula(0 ~ 1 + PT_plausibility*T_plausibility),times),
        "Ending" => (@formula(0 ~ 1 + PT_plausibility*T_plausibility),times)
    ]
    
    ch_ix = [15, 16, 17, 20, 21, 24, 25, 26, 28, 29, 30]
    #frontal: 4, 5, 6, 10, 11

    println("Unfolding dataset $i of $numDatasets: $setName")

    m_dc = fit(UnfoldModel, design,evts_clean,data[ch_ix,:];eventcolumn=:fix_type)
    c_dc = coef(m_dc)
    co_dc = coeftable(m_dc)
    #c2 = coef(m2)
    #co2 = coeftable(m2)

    println("Modeling epoched dataset $i of $numDatasets: $setName")

    m_epoch = fit(UnfoldModel,design_epoch,evts_clean,data_epoch[ch_ix,:,:];eventcolumn=:fix_type)
    c_epoch = coef(m_epoch)
    co_epoch = coeftable(m_epoch)

    #plot_erp(co)
    #plot_erp(co2)
    #plot_erp(co2,mapping=(;col=:eventname),legend=(;position=:bottom,orientation=:horizontal))
    
    println("Converting model output to effect coded estimates for $i of $numDatasets: $setName")

    _eff_dc = effects(Dict(:PT_plausibility=>["A","P"],:T_plausibility=>["A","P"]),m_dc)
    eff_dc = subset(_eff_dc,:eventname =>  x -> x .!=["onset"])
    #plot_erp(eff,mapping=(;color=:PT_plausibility,col=:fix_type))
    #plot_erp(eff_epoch,mapping=(;color=:PT_plausibility,col=:fix_type))

    _eff_epoch = effects(Dict(:PT_plausibility=>["A","P"],:T_plausibility=>["A","P"]),m_epoch)
    eff_epoch = subset(_eff_epoch,:eventname =>  x -> x .!=["onset"])

    eff_dc.group .= "deconv"
    eff_epoch.group .= "epoch"
    eff_all = vcat(eff_dc,eff_epoch)
    eff_all.subject .= setName

    println("Plotting model results for $i of $numDatasets: $setName")
    plot_erp(eff_all,mapping=(;color=:PT_plausibility,row=:T_plausibility,col=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])))
    
    println("Saving $setName deconvolution model results")
    save(joinpath(save_path_dc, setName * "_dc_model.jld2"), m_dc; compress = true);

    println("Saving $setName epoch model results")
    save(joinpath(save_path_epoch, setName * "_epoch_model.jld2"), m_epoch; compress = true);
    
end

#eff_all_onset = eff_all
# ERPExplorer.jl - discover(m)
# TidierPlots.jl (ggplot)
# TidierData.jl (or Tidier.jl for everything)

