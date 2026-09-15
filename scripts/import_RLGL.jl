using Unfold
using WGLMakie, UnfoldMakie # plotting
using Revise
using TidierData
includet("../src/io_eeglab.jl")

datasetFolder = "/Users/saramilligan/N4igami/RLGL_cont_data/";
datasetFiles = dir(fullfile(datasetFolder, '*.set')); % Get all .set files
numDatasets = length(datasetFiles);
filepath = "/Users/saramilligan/N4igami/RLGL_cont_data/RLGL033_BPF_AR_processed.set"

data, evts, eeg_dict = import_eeglab(filepath)

_evts = DataFrame(@chain evts begin
    @mutate(Sac_amp_in = lag(sac_amplitude,1),
            Sac_amp_out = lead(sac_amplitude,1))
    @filter(bdf_condition_description in ("ID", "ON", "TH") || type in ("S221", "S222", "S223"))
    @mutate(type = ifelse.(type in ("S221", "S222", "S223"), "onset", current_region))
end)

#=remove_onset = []
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
evts_clean = evts[Not(remove_onset),:]=#

#evts_clean[evts_clean.fix_type .=="Beginning","type"] .= "other_fixation"
design = [
    "onset" => (@formula(0 ~ 1),firbasis([-0.2,1.5],500))
    "Beginning" => (@formula(0 ~ 1),firbasis([-0.2,1],500)),
    "PreTarget" => (@formula(0 ~ 1 + bdf_condition_description * bdf_label_description + Sac_amp_in + Sac_amp_out),firbasis([-0.2,1],500)),
    "Target_word" => (@formula(0 ~ 1 + bdf_condition_description + Sac_amp_in + Sac_amp_out),firbasis([-0.2,1],500)),
    "Ending" => (@formula(0 ~ 1 + bdf_condition_description + Sac_amp_in + Sac_amp_out),firbasis([-0.2,1],500))
]

#"PreTarget" => (@formula(0 ~ 1 + bdf_condition_description * bdf_label_description + Sac_amp_in + Sac_amp_out),firbasis([-0.2,0.8],500)),

data_epoch,times=Unfold.epoch(;data=data,evts=_evts,τ=(-0.2,0.8),sfreq=500)
design_epoch = [
    "Beginning" => (@formula(0 ~ 1),times),
    "PreTarget" => (@formula(0 ~ 1 + bdf_condition_description * bdf_label_description + Sac_amp_in + Sac_amp_out),times),
    "Target_word" => (@formula(0 ~ 1 + bdf_condition_description + Sac_amp_in + Sac_amp_out),times),
    "Ending" => (@formula(0 ~ 1 + bdf_condition_description + Sac_amp_in + Sac_amp_out),times),
    "onset" => (@formula(0 ~ 1),times)
]


#"PreTarget" => (@formula(0 ~ 1 + bdf_condition_description * bdf_label_description + Sac_amp_in + Sac_amp_out),times),

ch_ix = 16
m_epoch = fit(UnfoldModel,design_epoch,_evts,data_epoch[ch_ix,:,:];eventcolumn=:type)


#m = fit(UnfoldModel, @formula(0~1),evts,data[16,:],firbasis([-0.2,0.8],500))
#c = coef(m)
#co = coeftable(m)

m_dc = fit(UnfoldModel, design,_evts, data[ch_ix,:];eventcolumn=:type)

save("test.jld2", m_dc; compress=true)
c_dc = coef(m_dc)
co_dc = coeftable(m_dc)

_eff = effects(Dict(:bdf_condition_description=>["ID","ON", "TH"],:bdf_label_description=>["PT_Fix", "PT_Skip"]),m_dc)
eff = subset(_eff,:eventname =>  x -> x .==["PreTarget"])
#plot_erp(eff,mapping=(;color=:PT_plausibility,col=:fix_type))
#plot_erp(eff_epoch,mapping=(;color=:PT_plausibility,col=:fix_type))

_eff_epoch = effects(Dict(:bdf_condition_description=>["ID","ON", "TH"],:bdf_label_description=>["PT_Fix", "PT_Skip"]),m_epoch)
eff_epoch = subset(_eff_epoch,:eventname =>  x -> x .==["PreTarget"])

eff.group .= "deconv"
eff_epoch.group .= "epoch"
eff_all = vcat(eff,eff_epoch)
plot_erp(eff_all,mapping=(;color=:bdf_condition_description,row=:bdf_label_description,col=:group=>UnfoldMakie.AlgebraOfGraphics.sorter(["epoch","deconv"])))


#eff_all_onset = eff_all

# ERPExplorer.jl - discover(m)
# TidierPlots.jl (ggplot)
# TidierData.jl (or Tidier.jl for everything)

