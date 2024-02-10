using Pumas, CSV, CairoMakie

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

fpm = fit(
  tgd_os_model,
  tgd_os_pop[1:nsubj],
  init_params(tgd_os_model),
  LaplaceI(),
)

validation_ll = loglikelihood(
  tgd_os_model,
  tgd_os_pop[nsubj+1:end],
  coef(fpm),
  LaplaceI(),
)
