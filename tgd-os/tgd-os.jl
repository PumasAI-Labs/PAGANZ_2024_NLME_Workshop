using Pumas, CSV

# Ref: https://www.sciencedirect.com/science/article/pii/S2451865421001289

df = CSV.read("tgd-os/tgd-os-data.csv", DataFrame; missingstring=".")

tgd_os_model = @model begin
    @params begin
        # Typical values of the coefficients of the predictors of SLD
        β ∈ VectorDomain(3)
        # IIV of the coefficients of the predictors of SLD
        Ω ∈ PSDDomain(2)
        # SLD residual error
        σ ∈ RealDomain(lower = 0.0)

        # Log logistic hazard function parameters
        h0 ∈ RealDomain(lower = 0.0)
        κ ∈ RealDomain(; lower=0.0)
        # Effect of SLD on hazard
        α ∈ RealDomain()
        # Coefficients of all the other predictors of hazard
        γ ∈ VectorDomain(2)
    end
    @random begin
        b ~ MvNormal(Ω)
    end
    @covariates begin
        NILSN # number of lesions
        AGE
        SEX
    end
    @pre begin
        # mean SLD
        m = dot(β, [1, t, NILSN]) + dot(b, [1, t])
        # 
        h = h0 * exp(dot(γ, [AGE, SEX]) + α * m)
        # log logistic hazard instanenous hazard
        λ := h * κ * (h * t + 1e-10)^(κ - 1) / (1 + (h * t + 1e-10)^κ)
    end
    @dynamics begin
        # cumulative hazard
        Λ' = λ
    end
    @derived begin
        SLD ~ @. Normal(m, σ)
        Death ~ TimeToEvent(λ, Λ)
    end
end

