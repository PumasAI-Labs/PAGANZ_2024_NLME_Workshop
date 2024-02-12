# Pumas Docs: https://docs.pumas.ai/stable/
using Pumas, CSV, DataFrames, DataFramesMeta, PumasUtilities, AlgebraOfGraphics, CategoricalArrays, Random, CairoMakie
const AOG = AlgebraOfGraphics

## Data wrangling ##
# Docs: https://dataframes.juliadata.org/stable/

# change the working directory
if isdir("PAGANZ2024")
    cd("PAGANZ2024")
end
# Current working directory
pwd()

df = CSV.read("warfarin/warfarin_data.csv", DataFrame; missingstring=".")
# avoid duplicate time points for the same DVID = 1
@. df[[133, 135, 137, 139],:TIME] += 1e-6 

@rtransform! df begin
    :FSZV = :WEIGHT / 70
    :FSZCL = (:WEIGHT / 70)^0.75
    :DVNAME = "DV$(:DVID)"
    :CMT = ismissing(:AMOUNT) ? missing : 1
    :EVID = ismissing(:AMOUNT) ? 0 : 1
end
df2 = unstack(
    df,
    Not([:DVID,:DVNAME,:DV]),
    :DVNAME,
    :DV,
)
df2
print(names(df2))
rename!(df2, :DV1 => :conc, :DV2 => :pca)

## Reading data ##

pop = read_pumas(
    df2;
    id = :ID,
    time = :TIME,
    amt = :AMOUNT,
    cmt = :CMT,
    evid = :EVID,
    covariates = [:SEX, :WEIGHT, :FSZV, :FSZCL],
    observations = [:conc, :pca],
)

## Model definition ##

model = @model begin
    @param begin
        # PK parameters
        "Clearance (L/h/70kg)"
        pop_CL   ∈ RealDomain(lower = 0.0, init = 0.134) 
        "Central Volume L/70kg"
        pop_V    ∈ RealDomain(lower = 0.0, init = 8.11)  
        "Absorption time (h)"
        pop_tabs ∈ RealDomain(lower = 0.0, init = 0.523) 
        "Lag time (h)"
        pop_lag  ∈ RealDomain(lower = 0.0, init = 0.1)   
        # PD parameters
        "Baseline"
        pop_e0   ∈ RealDomain(lower = 0.0, init = 100.0)
        "Emax"
        pop_emax ∈ RealDomain(init = -1.0)
        "EC50"
        pop_c50  ∈ RealDomain(lower = 0.0, init = 1.0)
        "Turnover"
        pop_tover  ∈ RealDomain(lower = 0.0, init = 14.0)
        # Inter-individual variability
        """
        - ΩCL
        - ΩVc
        - ΩTabs
        """
        pk_Ω     ∈ PDiagDomain([0.01, 0.01, 0.01])
        "Ωlag"
        lag_ω    ∈ RealDomain(lower = 0.0, init = 0.1) # unitless
        """
        - Ωe0
        - Ωemax
        - Ωec50
        - Ωturn
        """
        pd_Ω     ∈ PDiagDomain([0.01, 0.01, 0.01, 0.01]) 
        # Residual variability
        "Proportional Residual Error"
        σ_prop   ∈ RealDomain(lower = 0.0, init = 0.00752)
        "Additive Residual Error (mg/L)"
        σ_add    ∈ RealDomain(lower = 0.0, init = 0.0661)
        "Additive Error for FX"
        σ_fx     ∈ RealDomain(lower = 0.0, init = 0.01) 
    end

    @random begin
        # mean = 0, covariance = pk_Ω
        pk_η ~ MvNormal(pk_Ω)
        # mean = 0, standard deviation = lag_ω
        lag_η ~ Normal(0.0, lag_ω)
        # mean = 0, covariance = pd_Ω
        pd_η ~ MvNormal(pd_Ω)
    end

    @covariates FSZV FSZCL

    @pre begin
        # PK
        CL = FSZCL * pop_CL * exp(pk_η[1])
        Vc = FSZV * pop_V * exp(pk_η[2])
        tabs = pop_tabs * exp(pk_η[3])
        Ka = log(2) / tabs
        # PD
        e0 = pop_e0 * exp(pd_η[1])
        emax = pop_emax * exp(pd_η[2])
        c50 = pop_c50 * exp(pd_η[3])
        tover = pop_tover * exp(pd_η[4])
        kout = log(2) / tover
        rin = e0 * kout
        time = t
    end

    @dosecontrol begin
        lags = (Depot = pop_lag * exp(lag_η),)
    end

    @init begin
        Turnover = e0
    end

    # aliases for use in @dynamics and @derived
    @vars begin
        cp := Central / Vc
        ratein := Ka * Depot
        pd := 1 + emax * cp / (c50 + cp)
    end

    @dynamics begin
        Depot'    = -ratein
        Central'  =  ratein - CL * cp
        Turnover' =  rin * pd - kout * Turnover
    end

    @derived begin
        "Warfaring Concentration (mg/L)"
        conc ~ @. Normal(cp, sqrt((σ_prop * cp)^2 + σ_add^2))
        "PCA"
        pca  ~ @. Normal(Turnover, σ_fx)
    end
end

## Fit ##

# For most functions, the order of arguments is one of the following:
# model, data, params, [alg]
# model, data, params, randeffs

# Algorithms: FO(), FOCE(), LaplaceI()
fpm = fit(model, pop, init_params(model), FOCE())
# resets the inverse Hessian approximation of BFGS
fpm = fit(
    model,
    pop,
    coef(fpm),
    FOCE();
    init_randeffs = empirical_bayes(fpm),
)
coef(fpm)
coeftable(fpm)
coefficients_table(fpm)
ics = icoef(fpm)
ebes = empirical_bayes(fpm)

nlls = findinfluential(fpm)
@rsubset df2 :ID == "15"

## NONMEM OFV ##

loglikelihood(fpm)
n = nobs(pop)
OFV = (-2 * loglikelihood(fpm) - n * log(2π))

# Remove the lag random effect and use LaplaceI()
# Why remove the lag random effect?
# See the refs in https://www.page-meeting.org/?abstract=10625
fpm_laplace = fit(
    model,
    pop,
    coef(fpm),
    LaplaceI();
    init_randeffs = empirical_bayes(fpm),
    constantcoef = (lag_ω = 0.0,)
)

OFV_laplace = (-2 * loglikelihood(fpm_laplace) - n * log(2π))

## Visual predictive check ##

vpc_res_conc = vpc(fpm; 
                    observations = [:conc],
                    ensemblealg = EnsembleThreads())
warf_vpc = vpc_plot(
    vpc_res_conc;
    simquantile_medians = true,
    observations = false,
    axis = (xlabel = "Time (h)", ylabel = "Warfarin Concentration (mg/L)",
            xticks = 0:12:120)
)
figurelegend(warf_vpc, position = :b, orientation = :horizontal, nbanks = 3, tellwidth = true)
warf_vpc

vpc_res_pca = vpc(fpm; 
                  observations = [:pca],
                  ensemblealg = EnsembleThreads())

pca_vpc = vpc_plot(
    vpc_res_pca;
    simquantile_medians = true,
    observations = false,
    axis=(xlabel="Time (h)", ylabel="PCA",
        xticks=0:12:150)
)
figurelegend(pca_vpc, position=:b, orientation=:horizontal, nbanks=3, tellwidth=true)
pca_vpc

## Model diagnostics ##

insp = inspect(fpm)
insp_df = DataFrame(insp)
@transform! insp_df begin
    :SEXC = recode(:SEX, 0 => "female", 1 => "male")
end

## AlgebraOfGraphics ##
# Tutorials: https://tutorials.pumas.ai/html/PlottingInJulia/
# Docs: https://aog.makie.org/stable/

plt_df_conc = dropmissing(insp_df, [:conc])
data(plt_df_conc) * mapping(
    :conc_ipred => "Individual prediction",
    :conc => "Observation";
    layout = :SEXC,
) * (
    AOG.linear(interval = nothing) + 
    visual(AOG.Scatter)
) |> draw

plt_df_pca = dropmissing(insp_df, [:pca])
data(plt_df_pca) * mapping(
    :pca_ipred => "Individual prediction",
    :pca => "Observation";
    color = :SEXC,
    layout = :SEXC,
) * (
    AOG.linear(interval = nothing) + 
    visual(AOG.Scatter)
) |> draw

# Docs: https://docs.pumas.ai/dev/docstrings/pumasplots_docstrings/
# Observation plots
observations_vs_time(insp)
observations_vs_time(pop[1])
observations_vs_time(pop)

# Prediction plots
subject_fits(
    insp,
    separate = true,
    ids = unique(df2.ID)[1:12],
    observations = [:conc],
    facet = (combinelabels = true, )
)
subject_fits(
    insp,
    separate = true,
    ids = unique(df2.ID)[1:12],
    observations = [:pca],
    facet=(combinelabels=true,)
)
observations_vs_predictions(insp)
observations_vs_ipredictions(insp)

# Residual plots
wresiduals_vs_time(insp)
wresiduals_vs_predictions(insp)
wresiduals_vs_covariates(insp)
wresiduals_dist(insp)

## Non-gaussian random effects ##
# Docs: https://docs.pumas.ai/stable/basics/models/#@random:-Random-effects

model2 = @model begin
    @param begin
        # PK parameters
        pop_CL   ∈ RealDomain(lower = 0.0, init = 0.134) # L/h/70kg
        pop_V    ∈ RealDomain(lower = 0.0, init = 8.11)  # L/70kg
        pop_tabs ∈ RealDomain(lower = 0.0, init = 0.523) # h
        pop_lag  ∈ RealDomain(lower = 0.0, init = 0.1)   # h
        # PD parameters
        pop_e0   ∈ RealDomain(lower = 0.0, init = 100.0)
        pop_emax ∈ RealDomain(init = -1.0)
        pop_c50  ∈ RealDomain(lower = 0.0, init = 1.0)
        pop_tover  ∈ RealDomain(lower = 0.0, init = 14.0)
        # Inter-individual variability
        pk_Ω     ∈ PDiagDomain([0.01, 0.01, 0.01]) # unitless
        lag_ω    ∈ RealDomain(lower = 0.0, init = 0.1) # unitless
        pd_Ω     ∈ PDiagDomain([0.01, 0.01, 0.01, 0.01]) # unitless
        # Residual variability
        σ_prop   ∈ RealDomain(lower = 0.0, init = 0.00752) # unitless
        σ_add    ∈ RealDomain(lower = 0.0, init = 0.0661) # mg/L
        σ_fx     ∈ RealDomain(lower = 0.0, init = 0.01) # unitless
    end

    @random begin
        pk_η ~ MvLogNormal(pk_Ω)
        lag_η ~ LogNormal(0.0, lag_ω)
        pd_η ~ MvLogNormal(pd_Ω)
    end

    @covariates FSZV FSZCL

    @pre begin
        # PK
        CL = FSZCL * pop_CL * pk_η[1]
        Vc = FSZV * pop_V * pk_η[2]
        tabs = pop_tabs * pk_η[3]
        Ka = log(2) / tabs
        # PD
        e0 = pop_e0 * pd_η[1]
        emax = pop_emax * pd_η[2]
        c50 = pop_c50 * pd_η[3]
        tover = pop_tover * pd_η[4]
        kout = log(2) / tover
        rin = e0 * kout
    end

    @dosecontrol begin
        lags = (Depot = pop_lag * lag_η,)
    end

    @init begin
        Turnover = e0
    end

    # aliases for use in @dynamics and @derived
    @vars begin
        cp := Central / Vc
        ratein := Ka * Depot
        pd := 1 + emax * cp / (c50 + cp)
    end

    @dynamics begin
        Depot'    = -ratein
        Central'  =  ratein - CL * cp
        Turnover' =  rin * pd - kout * Turnover
    end

    @derived begin
        conc ~ @. Normal(cp, sqrt((σ_prop * cp)^2 + σ_add^2))
        pca  ~ @. Normal(Turnover, σ_fx)
    end
end

model3 = @model begin
    @param begin
        # PK parameters
        pop_CL   ∈ RealDomain(lower = 0.0, init = 0.134) # L/h/70kg
        pop_V    ∈ RealDomain(lower = 0.0, init = 8.11)  # L/70kg
        pop_tabs ∈ RealDomain(lower = 0.0, init = 0.523) # h
        pop_lag  ∈ RealDomain(lower = 0.0, init = 0.1)   # h
        # PD parameters
        pop_e0   ∈ RealDomain(lower = 0.0, init = 100.0)
        pop_emax ∈ RealDomain(init = -1.0)
        pop_c50  ∈ RealDomain(lower = 0.0, init = 1.0)
        pop_tover  ∈ RealDomain(lower = 0.0, init = 14.0)
        # Inter-individual variability
        pk_Ω     ∈ PDiagDomain([0.01, 0.01, 0.01]) # unitless
        lag_ω    ∈ RealDomain(lower = 0.0, init = 0.1) # unitless
        pd_Ω     ∈ PDiagDomain([0.01, 0.01, 0.01, 0.01]) # unitless
        # Residual variability
        σ_prop   ∈ RealDomain(lower = 0.0, init = 0.00752) # unitless
        σ_add    ∈ RealDomain(lower = 0.0, init = 0.0661) # mg/L
        σ_fx     ∈ RealDomain(lower = 0.0, init = 0.01) # unitless
    end

    @random begin
        pk_η ~ MvLogNormal(
            log.([pop_CL, pop_V, pop_tabs]),
            pk_Ω,
        )
        lag_η ~ LogNormal(
            log(pop_lag),
            lag_ω,
        )
        pd_η ~ MvLogNormal(
            log.([pop_e0, 1.0, pop_c50, pop_tover]),
            pd_Ω,
        )
    end

    @covariates FSZV FSZCL

    @pre begin
        # PK
        CL = FSZCL * pk_η[1]
        Vc = FSZV * pk_η[2]
        tabs = pk_η[3]
        Ka = log(2) / tabs
        # PD
        e0 = pd_η[1]
        emax = pop_emax * pd_η[2]
        c50 = pd_η[3]
        tover = pd_η[4]
        kout = log(2) / tover
        rin = e0 * kout
    end

    @dosecontrol begin
        lags = (Depot = lag_η,)
    end

    @init begin
        Turnover = e0
    end

    # aliases for use in @dynamics and @derived
    @vars begin
        cp := Central / Vc
        ratein := Ka * Depot
        pd := 1 + emax * cp / (c50 + cp)
    end

    @dynamics begin
        Depot'    = -ratein
        Central'  =  ratein - CL * cp
        Turnover' =  rin * pd - kout * Turnover
    end

    @derived begin
        conc ~ @. Normal(cp, sqrt((σ_prop * cp)^2 + σ_add^2))
        pca  ~ @. Normal(Turnover, σ_fx)
    end
end

## Uncertainty quantification ##

# Sampling distribution covariance matrix estimation
inf = infer(fpm) # fails
# Bootstrap
bts_inf = infer(fpm, Bootstrap(samples = 20))

## Simulation ##

# Sample random effects from their priors
sims1 = simobs(model, pop, coef(fpm))

# Simulate 100 populations
sims2 = [simobs(model, pop, coef(fpm)) for _ in 1:100]

# Small time step simulation
fine_sims = simobs(model, pop, coef(fpm), obstimes = 0.0:0.5:150.0)
fine_sims[1]

# Plot the simulation curves
sim_plot(fine_sims, observations = [:conc])
sim_plot(fine_sims, observations = [:pca])

# Use the empirical Bayes estimates
sims3 = simobs(model, pop, coef(fpm), empirical_bayes(fpm))
# or
sims3 = simobs(fpm)

# Fix the random number generator
rng = Pumas.default_rng()
Random.seed!(rng, 123)
simobs(fpm; rng = rng)
# or
simobs(fpm; rng)

# Sample population parameters from the bootstrap samples first
# Then sample random effects from their priors
sims4 = simobs(bts_inf, pop, samples = 10)
# or
@time sims4 = [simobs(bts_inf, pop, samples = 1)[1] for _ in 1:1000]

# VPC from simulations
vpc_res_conc = vpc(sims4; observations = [:conc])
vpc_plot(
    vpc_res_conc;
    simquantile_medians = true,
    observations = false,
)

vpc_res_pca = vpc(sims4; observations = [:pca])
vpc_plot(
    vpc_res_pca;
    simquantile_medians = true,
    observations = false,
)

## Prediction ##

ipreds = predict(fpm)
ipreds_df = DataFrame(ipreds)
print(names(ipreds_df))
ipreds_df[:,[:id, :time, :conc_ipred, :pca_ipred]]

## Bayesian workflow ##
# Reference: https://arxiv.org/abs/2304.04752

bayes_model = @model begin
    @param begin
        # PK parameters
        pop_CL   ~ LogNormal(log(0.134), 1.0) # L/h/70kg
        pop_V    ~ LogNormal(log(8.11), 1.0)  # L/70kg
        pop_tabs ~ LogNormal(log(0.523), 1.0) # h
        # PD parameters
        pop_e0   ~ LogNormal(log(100.0), 1.0)
        pop_emax ~ Normal(-1.0, 1.0)
        pop_c50  ~ LogNormal(log(1.0), 1.0)
        pop_tover ~ LogNormal(log(14.0), 1.0)
        # Inter-individual variability
        pk_Ω     ~ Constrained(MvNormal([0.01, 0.01, 0.01], 1.0), lower = 0.0) # unitless
        pd_Ω     ~ Constrained(MvNormal([0.01, 0.01, 0.01, 0.01], 1.0), lower = 0.0) # unitless
        # Residual variability
        σ_prop   ~ Constrained(Normal(0.00752, 0.1), lower = 0.0) # unitless
        σ_add    ~ Constrained(Normal(0.0661, 1.0), lower = 0.0) # mg/L
        σ_fx     ~ Constrained(Normal(0.01, 1.0), lower = 0.0) # unitless
    end

    @random begin
        # mean = 0, covariance = pk_Ω
        pk_η ~ MvNormal(Diagonal(pk_Ω))
        # mean = 0, covariance = pd_Ω
        pd_η ~ MvNormal(Diagonal(pd_Ω))
    end

    @covariates FSZV FSZCL

    @pre begin
        # PK
        CL = FSZCL * pop_CL * exp(pk_η[1])
        Vc = FSZV * pop_V * exp(pk_η[2])
        tabs = pop_tabs * exp(pk_η[3])
        Ka = log(2) / tabs
        # PD
        e0 = pop_e0 * exp(pd_η[1])
        emax = pop_emax * exp(pd_η[2])
        c50 = pop_c50 * exp(pd_η[3])
        tover = pop_tover * exp(pd_η[4])
        kout = log(2) / tover
        rin = e0 * kout
        time = t
    end

    # @dosecontrol begin
    #     lags = (Depot = pop_lag * exp(lag_η),)
    # end

    @init begin
        Turnover = e0
    end

    # aliases for use in @dynamics and @derived
    @vars begin
        cp := Central / Vc
        ratein := Ka * Depot
        pd := 1 + emax * cp / (c50 + cp)
    end

    @dynamics begin
        Depot'    = -ratein
        Central'  =  ratein - CL * cp
        Turnover' =  rin * pd - kout * Turnover
    end

    @derived begin
        conc ~ @. Normal(cp, sqrt((σ_prop * cp)^2 + σ_add^2))
        pca  ~ @. Normal(Turnover, σ_fx)
    end
end

bayes_fpm = fit(
    bayes_model,
    pop,
    init_params(bayes_model),
    BayesMCMC(
        nsamples = 500,
        nadapts = 250,
        nchains = 4,
        # remove the following option in a real run!
        alg = GeneralizedNUTS(max_depth = 3),
        diffeq_options = (; alg = Pumas.Rodas5())
    ),
);
bayes_fpm_samples = Pumas.discard(bayes_fpm, burnin = 250)
mean(bayes_fpm_samples)
coef(fpm)

density_plot(bayes_fpm_samples, parameters = [:pop_tover])

## Post-processing, AUC and Cmax ##

# Get AUC and NCA given simulations
nca_params = postprocess(reduce(vcat, sims2)) do gen, obs
    pk_auc = NCA.auc(gen.conc, gen.time)
    pk_cmax = NCA.cmax(gen.conc, gen.time)
    pd_auc = NCA.auc(gen.pca, gen.time)
    pd_cmax = NCA.cmax(gen.pca, gen.time)
    (; pk_auc, pk_cmax, pd_auc, pd_cmax)
end
# Estimate the probability of being in a therapeutic range
# Can choose the dose to maximize this probability
prob = mean(nca_params) do p
    !ismissing(p.pk_auc) && 500.0 < p.pk_auc < 1000.0 && p.pk_cmax <= 15
end

getproperty.(nca_params, :pk_auc) |> skipmissing |> collect |> hist
getproperty.(nca_params, :pk_cmax) |> skipmissing |> collect |> hist
getproperty.(nca_params, :pd_auc) |> skipmissing |> collect |> hist
getproperty.(nca_params, :pd_cmax) |> skipmissing |> collect |> hist

# Counter-factual simulation
# New dose simulation for an existing subject
function new_dose_sim(bayes_fpm_samples, i, nsim, dose_multiplier, obstimes = 0.0:1.0:150.0)
    ids = unique(df2.ID)
    subj_df = @rsubset df2 :ID == ids[i]
    subj_df[1, :AMOUNT] = subj_df[1, :AMOUNT] * dose_multiplier
    subj_df
    df2
    new_subj = read_pumas(
        subj_df;
        id = :ID,
        time = :TIME,
        amt = :AMOUNT,
        cmt = :CMT,
        evid = :EVID,
        covariates = [:SEX, :WEIGHT, :FSZV, :FSZCL],
        observations = [:conc, :pca],
    )[1]

    # Samples the populaton parameters from the posterior
    # Samples the random effects from the posterior of subject i
    # Uses the dose and covariates in `new_subj`
    sims = simobs(
        bayes_fpm_samples,
        new_subj;
        subject = i,
        samples = nsim,
        obstimes = obstimes,
    )
    bayes_nca_params = postprocess(reduce(vcat, sims)) do gen, obs
        pk_auc = NCA.auc(gen.conc, gen.time)
        pk_cmax = NCA.cmax(gen.conc, gen.time)
        pd_auc = NCA.auc(gen.pca, gen.time)
        pd_cmax = NCA.cmax(gen.pca, gen.time)
        (; pk_auc, pk_cmax, pd_auc, pd_cmax)
    end
    prob = mean(bayes_nca_params) do p
        # define desired therapeutic ranges or effects
        # the following numbers are made up!
        !ismissing(p.pk_auc) && 1000.0 < p.pk_auc < 2000.0 && p.pk_cmax <= 50
    end
    return sims, bayes_nca_params, prob
end

sims5, bayes_nca_params, prob = new_dose_sim(
    bayes_fpm_samples, 1, 100, 0.6,
)
prob

# Simulate a new subject altogether
# Samples the population parameters from the posterior
# Samples the random effects from the prior
sims6 = simobs(
    bayes_fpm_samples,
    # assume this is a new subject
    pop[1];
    samples = 100,
)

# Calculate probability of successful treatment for different doses
probs = map(0.1:0.1:3.0) do dose_multiplier
    new_dose_sim(bayes_fpm_samples, 1, 500, dose_multiplier)[3]
end
data((dose_multiplier = 0.1:0.1:3.0, probs = probs)) * 
    mapping(:dose_multiplier => "Dose multiple", :probs => "Probability of successful treatment") * 
    AOG.visual(Lines) |> draw
