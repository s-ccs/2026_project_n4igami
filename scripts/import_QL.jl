using Unfold
using WGLMakie, UnfoldMakie # plotting
using Revise

includet("../src/io_eeglab.jl")
filepath = "/Users/saramilligan/QL_processing/EEG_files/EyeSorted/QLF102_ICA_clean_processed.set"


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


design_epoch = [
    "other_fixation" => (@formula(0 ~ 1),times),
    "R_fixation" => (@formula(0 ~ 1 + PT_plausibility*fix_type),times),
    "onset" => (@formula(0 ~ 1),times)
]
#m = fit(UnfoldModel, @formula(0~1),evts,data[16,:],firbasis([-0.2,0.8],500))
#c = coef(m)
#co = coeftable(m)

m2 = fit(UnfoldModel, design,evts_clean,data[16,:];eventcolumn=:type)
c2 = coef(m2)
co2 = coeftable(m2)


plot_erp(co)
plot_erp(co2)
plot_erp(co2,mapping=(;col=:eventname),legend=(;position=:bottom,orientation=:horizontal))


_eff = effects(Dict(:PT_plausibility=>["A","P"],:fix_type=>["Pretarget","Target"]),m2)
eff = subset(_eff,:eventname =>  x -> x .==["R_fixation"])
plot_erp(eff,mapping=(;color=:PT_plausibility,col=:fix_type))

# ERPExplorer.jl - discover(m)
# TidierPlots.jl (ggplot)
# TidierData.jl (or Tidier.jl for everything)

