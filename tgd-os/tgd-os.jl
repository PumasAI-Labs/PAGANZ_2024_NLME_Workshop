using Pumas, CSV, CairoMakie, PumasUtilities

# change the working directory
if isdir("PAGANZ2024")
    cd("PAGANZ2024")
end
# Current working directory
pwd()

df = CSV.read("tgd-os/tgd-os2.csv", DataFrame; missingstring=[".", ""])
tgd_os_pop = read_pumas(
    df,
    observations = [:SLD, :Death],
    covariates = [:WT, :AGE, :SEX],
    event_data = false,
)

tgd_os_model = @model begin
    @param begin
        # Typical values of the coefficients of the predictors of SLD
        β ∈ VectorDomain(2, init = [1.0, -1.0])
        # IIV of the coefficients of the predictors of SLD
        Ω ∈ PSDDomain(init = [16.0 0.0; 0.0 1.0])
        # SLD residual error
        σ ∈ RealDomain(lower = 0.0, init = 1.0)

        # Log logistic hazard function parameters
        h0 ∈ RealDomain(lower = 0.0, init = 0.001)
        κ ∈ RealDomain(; lower=0.0, init = 1.1)
        # Effect of SLD on hazard
        α ∈ RealDomain(init = 0.001)
        # Coefficients of all the other predictors of hazard
        γ ∈ VectorDomain(2, init = [0.01, 0.01])
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @covariates begin
        AGE
        SEX
    end
    @pre begin
        # mean SLD
        m = (β[1] + η[1]) + (β[2] + η[2]) * t
        # log logistic instanenous hazard
        sexn = SEX == "Female" ? 0.0 : 1.0
        h = h0 * exp(γ[1] * AGE + γ[2] * sexn + α * m)
        λ = h * κ * (h * t + 1e-10)^(κ - 1) / (1 + (h * t + 1e-10)^κ)
    end
    @dynamics begin
        # cumulative hazard
        Λ' = λ
    end
    @derived begin
        SLD ~ @. Normal(m, σ)
        Death ~ @. TimeToEvent(λ, Λ)
    end
end

nsubj = round(Int, 0.8 * length(tgd_os_pop))
loglikelihood(tgd_os_model, tgd_os_pop[1:nsubj], init_params(tgd_os_model), LaplaceI())

tgd_os_fpm = fit(
  tgd_os_model,
  tgd_os_pop[1:nsubj],
  init_params(tgd_os_model),
  LaplaceI(),
)

validation_ll = loglikelihood(
  tgd_os_model,
  tgd_os_pop[nsubj+1:end],
  coef(tgd_os_fpm),
  LaplaceI(),
)

# Simulation/VPC workaround until next Pumas version

tgd_model = @model begin
    @param begin
        # Typical values of the coefficients of the predictors of SLD
        β ∈ VectorDomain(2, init = [1.0, -1.0])
        # IIV of the coefficients of the predictors of SLD
        Ω ∈ PSDDomain(init = [16.0 0.0; 0.0 1.0])
        # SLD residual error
        σ ∈ RealDomain(lower = 0.0, init = 1.0)

        # Log logistic hazard function parameters
        h0 ∈ RealDomain(lower = 0.0, init = 0.001)
        κ ∈ RealDomain(; lower=0.0, init = 1.1)
        # Effect of SLD on hazard
        α ∈ RealDomain(init = 0.001)
        # Coefficients of all the other predictors of hazard
        γ ∈ VectorDomain(2, init = [0.01, 0.01])
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @covariates begin
        AGE
        SEX
    end
    @pre begin
        # mean SLD
        m = (β[1] + η[1]) + (β[2] + η[2]) * t
        # log logistic instanenous hazard
        sexn = SEX == "Female" ? 0.0 : 1.0
        h = h0 * exp(γ[1] * AGE + γ[2] * sexn + α * m)
        λ = h * κ * (h * t + 1e-10)^(κ - 1) / (1 + (h * t + 1e-10)^κ)
    end
    @derived begin
        SLD ~ @. Normal(m, σ)
    end
end
tgd_pop = read_pumas(
    df,
    observations = [:SLD],
    covariates = [:WT, :AGE, :SEX],
    event_data = false,
)
tgd_fpm = fit(
    tgd_model,
    tgd_pop,
    coef(tgd_os_fpm),
    LaplaceI();
    checkidentification = false,
    optim_options = (iterations = 0,),
)
tgd_insp = inspect(tgd_fpm)
sf_sld = subject_fits(
    tgd_insp,
    separate = true,
    ids = string.(1:16),
    facet = (combinelabels = true,),
    axis=(xlabel = "Time (days)", xticklabelrotation=pi / 4,)
)
# figurelegend(sf_sld, orientation=:horizontal)
sf_sld

tgd_vpc_res = vpc(tgd_fpm, 
                  observations = [:SLD], 
                  bandwidth = 20,
                ensemblealg = EnsembleThreads())
# not a very good TGD model!
vpc_plot(
    tgd_vpc_res;
    simquantile_medians = true,
    observations = false,
    axis = (xlabel = "Time (days)", ylabel  = "Sum of Longest Diameter"),

)

os_model = @model begin
    @param begin
        # Typical values of the coefficients of the predictors of SLD
        β ∈ VectorDomain(2, init = [1.0, -1.0])
        # IIV of the coefficients of the predictors of SLD
        Ω ∈ PSDDomain(init = [16.0 0.0; 0.0 1.0])
        # SLD residual error
        σ ∈ RealDomain(lower = 0.0, init = 1.0)

        # Log logistic hazard function parameters
        h0 ∈ RealDomain(lower = 0.0, init = 0.001)
        κ ∈ RealDomain(; lower=0.0, init = 1.1)
        # Effect of SLD on hazard
        α ∈ RealDomain(init = 0.001)
        # Coefficients of all the other predictors of hazard
        γ ∈ VectorDomain(2, init = [0.01, 0.01])
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @covariates begin
        AGE
        SEX
    end
    @pre begin
        # mean SLD
        m = (β[1] + η[1]) + (β[2] + η[2]) * t
        # log logistic instanenous hazard
        sexn = SEX == "Female" ? 0.0 : 1.0
        h = h0 * exp(γ[1] * AGE + γ[2] * sexn + α * m)
        λ = h * κ * (h * t + 1e-10)^(κ - 1) / (1 + (h * t + 1e-10)^κ)
    end
    @dynamics begin
        # cumulative hazard
        Λ' = λ
    end
    @derived begin
        Death ~ @. TimeToEvent(λ, Λ)
    end
end

os_pop = read_pumas(
    df,
    observations = [:Death],
    covariates = [:WT, :AGE, :SEX],
    event_data = false,
)
os_fpm = fit(
    os_model,
    os_pop,
    coef(tgd_os_fpm),
    LaplaceI();
    checkidentification = false,
    optim_options = (iterations = 0,),
)
vpc_res_os = vpc(os_fpm, 
                observations = [:Death],
                ensemblealg = EnsembleThreads())
vpc_plot(vpc_res_os,
        axis=(xlabel="Time (days)",),
        )
