using MAT, DataFrames
function import_eeglab(filepath)

    eeg_dict = MAT.matread(filepath)


    function converter(v)
        v2 = try 
                    Float64.(v)
             catch
                    string.(v)
             end
        return dropdims(v2,dims=1) 
    end
    evts = DataFrame([k=>converter(v) for (k,v) in eeg_dict["event"]])


    fdt_filename = eeg_dict["data"]

    # Construct full path to .fdt file (same directory as .set file)
    set_dir = dirname(filepath)
    fdt_path = joinpath(set_dir, fdt_filename)

    if !isfile(fdt_path)
        error("External data file not found: $fdt_path")
    end

    @info "Loading external data from: $fdt_filename"

    # Read binary .fdt file
    n_channels = Int(eeg_dict["nbchan"])
    n_timepoints = Int(eeg_dict["pnts"])
    n_trials = Int(eeg_dict["trials"])

    # Total expected elements
    n_elements = n_channels * n_timepoints * n_trials

    # Read binary file as Float32 array
    data_vec = open(fdt_path, "r") do io
        read!(io, Vector{Float32}(undef, n_elements))
    end

    # Reshape to [channels × timepoints × trials] or [channels × timepoints]
    if n_trials > 1
        return reshape(data_vec, n_channels, n_timepoints, n_trials), evts, eeg_dict
    else
        return reshape(data_vec, n_channels, n_timepoints), evts, eeg_dict
    end

end