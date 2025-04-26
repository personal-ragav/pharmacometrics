function ct_latex(C₀::Number, k::Number, Time::Number, Cu::AbstractString)

    post = (x -> string(round(x, digits = 3), "\\ ", Cu))
    fmt = (x -> format(round(x, digits=3)))
  
    eqn = @latexdefine Cₜ = $C₀ * ℯ ^ (-$k * $Time) env=:equation post=post fmt=fmt
    side_eqn = @latexdefine Time post = (x -> string(x, "\\ ", "hr"))
    lstr = LaTeXString(side_eqn * eqn)
  
    return Cₜ, lstr
  end
  
  function c0_latex(Cₜ::Number, k::Number, Time::Number, Cu::AbstractString)
    
    post = (x -> string(round(x, digits = 3), "\\ ", Cu))
    fmt = (x -> format(round(x, digits=3)))
  
    eqn = @latexdefine C₀ = $Cₜ * ℯ ^ ($k * $Time) env=:equation post=post fmt=fmt 
    side_eqn = @latexdefine Time post = (x -> string(x, "\\ ", "hr"))
    lstr = LaTeXString(side_eqn * eqn) # latexify("Time=$t") *
  
    return C₀, lstr
  end

  function k_latex(tₕ::Number)

    post = (x -> string(round(x, digits = 3), "\\ hr^{-1}"))
    fmt = (x -> format(round(x, digits=3)))
  
    eqn = @latexdefine k = 0.693 / $tₕ env=:eq post=post fmt=fmt 
    lstr = LaTeXString(latexify("k = 0.693 / t_h") * eqn)
  
    return k, lstr
  end
  
  function th_latex(Cₜ₁::Number, t_1::Number, Cₜ₂::Number, t_2::Number, Cu::String)
    post = (x -> string(round(x, digits = 3), "\\ ", Cu))
    fmt = (x -> format(round(x, digits = 3)))

    post1 = (x -> string(x, "\\ ", "hr", ", \\ ")) 

    eqn = @latexdefine tₕ = 0.693 * ($t_2 - $t_1) / (log($Cₜ₁) - log($Cₜ₂)) env=:equation post=post fmt=fmt
    side_eqn1 = @latexdefine t_1 post=post1
    side_eqn2 = @latexdefine t_2 post=post1
    
    lstr = LaTeXString(side_eqn1 * side_eqn2 * eqn)
  
    return tₕ, lstr
  end
  
  function vd_latex(Amount::Number,C₀::Number, Au::String, Vu::String)
    post = (x -> string(round(x, digits = 3), "\\ ", Vu))
    fmt = (x -> format(round(x, digits=3)))

    post1 = (x -> string(x, "\\ ", Au))

    eqn = @latexdefine V_d = $Amount / $C₀ env=:equation post=post fmt=fmt
    side_eqn = @latexdefine Amount post=post1
    lstr = LaTeXString(side_eqn * eqn)
  
    return V_d, lstr
  end
  
  function vd_cl_latex(Cl::Number,k::Number, Vu::String)
    post = (x -> string(round(x, digits=3), "\\ ", Vu))
    fmt = (x -> round(x, digits=3))

    eqn = @latexdefine V_d = $Cl / $k env=:equation post=post fmt=fmt
    lstr = LaTeXString(latexify("V_d = Cl / k") * eqn)
  
    return V_d, lstr
  end
  
  function cl_latex(vd::Number, k::Number, Vu::String)
    post = (x -> string(round(x, digits=3), "\\ ", "$Vu/hr"))
    fmt = (x -> format(round(x, digits=3)))
  
    eqn = @latexdefine Cl = $vd * $k env=:equation post=post fmt=fmt
    lstr = LaTeXString(latexify("Cl = V_d * k") * eqn)
  
    return Cl, lstr
  end
  
  function cl_auc_latex(Dose::Number, AUC::Number, Clu)
    post = (x -> string(round(x, digits=3), "\\ ", Clu))
    fmt = (x -> round(x, digits=3))
  
    eqn = @latexdefine Cl = $Dose / $AUC env=:equation post=post fmt=fmt
    lstr = LaTeXString(latexify("Cl = Dose / AUC") * eqn)
  
    return Cl, lstr
  end
  
  function ct_miss_latex(df::DataFrame, C₀::Number, k::Number, Cu::String)
  
    df_dict = Vector{Dict{Symbol, Float64}}()
    latexstrs = Vector{LaTeXString}()
  
    dfrow(row) = Dict(pairs(row))
  
    for row in eachrow(df)
        if ismissing(row.conc)
            ct, ct_lat = ct_latex(C₀, k, row.time, Cu)
            row.conc = ct
            push!(df_dict, dfrow(row))
            push!(latexstrs, ct_lat)
        else
            push!(df_dict, dfrow(row))
        end
    end
  
    df_new = @chain begin
        DataFrame(df_dict)
        select(propertynames(df))
        @orderby :time
    end
  
    return df_new, latexstrs
  
  end
  
  # Function to compute AUC using trapezoidal rule
  function trapz_auc(time::Vector{Float64}, conc::Vector{Float64})

    auci = Vector{Float64}()
    post = (x -> string(round(x, digits=3), "\\ ", "mg*h/L"))
    fmt = (x -> format(round(x, digits=3)))

    if length(time) != length(conc)
        error("Time and concentration vectors must be the same length")
    end
    
    auc_total = 0.0
    for i in 2:length(time)
        t1 = time[i]
        t2 = time[i-1]
        c1 = conc[i]
        c2 = conc[i-1]
        
        lstr = @latexdefine AUC = ($c1 + $c2) * ($t1-$t2) / 2 env=:equation

        push!(auci, AUC)

        auc_total += AUC
    end

    expr = :($(join(string.(round.(auci, digits=2)), " + "))) 
    equation = Meta.parse(string(expr))
    AUC_total = eval(equation)
    
    exp_eqn = latexify(string("AUC_total => ", equation))
    fin_eqn = @latexdefine AUC_total post=post fmt=fmt
    lstr = LaTeXString(exp_eqn * fin_eqn)

    return lstr, auc_total
end

function table(mdf::DataFrame, cell_border_bottom = true)
    C(value; kwargs...) =
        Cell(value; halign = :left, border_bottom = cell_border_bottom, kwargs...)

    headercells = C.(reshape(names(mdf), 1, ncol(mdf)), bold = true, border_bottom = true)

    # C(Multiline(ele)) 
    df_cells = @chain begin
        map(ele -> isa(ele, Vector) ? C(Multiline(ele...)) : C(ele), Matrix(mdf))
        reshape(nrow(mdf), ncol(mdf))
    end

    colgap_list = [i => 10.0 for i = 1:ncol(mdf)-1]
    rowgap_list = [i => 8.0 for i = 1:nrow(mdf)]

    cellmatdf = vcat(headercells, df_cells)
    table = Table(cellmatdf; colgaps = colgap_list, rowgaps = rowgap_list)

    return table
end