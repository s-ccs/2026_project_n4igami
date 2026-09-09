using Unfold
using WGLMakie, UnfoldMakie # plotting
using Revise

includet("../src/io_eeglab.jl")
filepath = "/Users/saramilligan/WTS/DC_modeling/output_deconvolution/high_constraint_WTS/continuous_processed_datasets_beg/UF_WTS032_ICA_eyesort_ia_processed.set"
filepath = "/Volumes/cas/PSY Schotter Lab/Experiments/WhatTheSpell/WTS_deconv/output_deconvolution/high_constraint_WTS/UF_modified_datasets_WTS/UF_WTS011_ICA_eyesort_ia_processed.set"



data, evts, eeg_dict = import_eeglab(filepath)



remove_onset = []
onset_active = false
for (idx,e) in reverse(collect(enumerate(evts.type)))
    if onset_active && e == "onset"
        push!(remove_onset,idx)
    elseif e == "onset" 
               onset_active = true
    else
        onset_active = false
    end
end
evts_clean = evts[Not(remove_onset),:]

evts_clean[evts_clean.fix_type .=="Beginning","type"] .= "other_fixation"
design = [
    "other_fixation" => (@formula(0 ~ 1),firbasis([-0.2,0.8],500)),
    "R_fixation" => (@formula(0 ~ 1 + PT_plausibility*fix_type),firbasis([-0.2,0.8],500)),
    "onset" => (@formula(0 ~ 1),firbasis([-0.2,1.5],500))
]


data_epoch,times=Unfold.epoch(;data=data,evts=evts_clean,τ=(-0.2,0.8),sfreq=500)
design_epoch = [
    "other_fixation" => (@formula(0 ~ 1),times),
    "R_fixation" => (@formula(0 ~ 1 + PT_plausibility*fix_type),times),
    "onset" => (@formula(0 ~ 1),times)
]
ch_ix = 16
m_epoch = fit(UnfoldModel,design_epoch,evts_clean,data_epoch[ch_ix,:,:];eventcolumn=:type)

#m = fit(UnfoldModel, @formula(0~1),evts,data[16,:],firbasis([-0.2,0.8],500))
#c = coef(m)
#co = coeftable(m)

m2 = fit(UnfoldModel, design,evts_clean,data[ch_ix,:];eventcolumn=:type)
#c2 = coef(m2)
#co2 = coeftable(m2)


#plot_erp(co)
#plot_erp(co2)
#plot_erp(co2,mapping=(;col=:eventname),legend=(;position=:bottom,orientation=:horizontal))


_eff = effects(Dict(:PT_plausibility=>["A","P"],:fix_type=>["Pretarget","Target"]),m_epoch)
eff_epoch = subset(_eff,:eventname =>  x -> x .==["R_fixation"])

_eff = effects(Dict(:PT_plausibility=>["A","P"],:fix_type=>["Pretarget","Target"]),m2)
eff = subset(_eff,:eventname =>  x -> x .==["R_fixation"])
#plot_erp(eff,mapping=(;color=:PT_plausibility,col=:fix_type))
#plot_erp(eff_epoch,mapping=(;color=:PT_plausibility,col=:fix_type))


eff.group .= "deconv"
eff_epoch.group .= "epoch"
eff_all = vcat(eff,eff_epoch)
plot_erp(eff_all,mapping=(;color=:PT_plausibility,row=:fix_type,col=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])))



eff_all_onset = eff_all
# ERPExplorer.jl - discover(m)
# TidierPlots.jl (ggplot)
# TidierData.jl (or Tidier.jl for everything)

