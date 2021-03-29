library(tidyverse)
library(SmoothHte)
library(SimulateHte)
library(SimulationEvaluationHte)
library(dplyr)

set.seed(19910930)

analysisIds <- readr::read_csv(
  "~/Documents/Projects/arekkas_HteSimulation_XXXX_2021/data/processed/analysisIds.csv",
  col_types = "iffiiddddddddddddd"
)

analysisIdsInteractions <- readr::read_csv(
  "~/Documents/Projects/arekkas_HteSimulation_XXXX_2021/data/processed/analysisIdsInteractions.csv",
  col_types = "icfiidddddddddddddd"
)

selectedScenario <- 27

if (selectedScenario <= 63) {
  idSettings <- analysisIds %>%
    filter(scenario == selectedScenario)
  createF1 <- function(c) function(x) x - c
  createF2 <- function(c) function(x) (x - c)^2
  treatmentEffectSettings <- SimulateHte::createTreatmentEffectSettings(
    type = "lp",
    modelSettings = SimulateHte::createModelSettings(
      constant = idSettings$g0,
      modelMatrix = matrix(c(1, 1)),
      transformationSettings = list(
        createF1(idSettings$c),
        createF2(idSettings$c)
      ),
      coefficients = c(
        idSettings$g1,
        idSettings$g2
      )
    )
  )
} else {
  idSettings <- analysisIdsInteractions %>%
    filter(scenario == selectedScenario)
  treatmentEffectSettings <- SimulateHte::createTreatmentEffectSettings(
    type = "covariates",
    modelSettings = SimulateHte::createModelSettings(
      constant = idSettings$g0,
      modelMatrix = rbind(
        c(1, rep(0, 7)),
        c(0, 1, rep(0, 6)),
        c(rep(0, 4), 1, 0, 0, 0),
        c(rep(0, 5), 1, 0, 0)
      ),
      transformationSettings = list(
        identity,
        identity,
        identity,
        identity
      ),
      coefficients = idSettings %>% select(matches("g[1-9]")) %>% unlist()
    )
  )
}
# idSettings$sampleSize <- 5e4

databaseSettings <- SimulateHte::createDatabaseSettings(
  numberOfObservations = as.numeric(as.character(idSettings$sampleSize)),
  numberOfCovariates = 8,
  covariateDistributionSettings = list(
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createBinomialDistributionSettings(prob = .2),
    SimulateHte::createBinomialDistributionSettings(prob = .2),
    SimulateHte::createBinomialDistributionSettings(prob = .2),
    SimulateHte::createBinomialDistributionSettings(prob = .2)
  )
)

baselineRiskSettings <- SimulateHte::createBaselineRiskSettings(
  type = "binary",
  modelSettings = SimulateHte::createModelSettings(
    constant = idSettings %>% pull(b0),
    modelMatrix = diag(8),
    transformationSettings = list(
      identity,
      identity,
      identity,
      identity,
      identity,
      identity,
      identity,
      identity
    ),
    coefficients = idSettings %>% select(paste0("b", 1:8)) %>% unlist()
  )
)


propensitySettings <- SimulateHte::createPropensitySettings(
  type = "binary",
  modelSettings = SimulateHte::createModelSettings(
    constant = 0,
    modelMatrix = diag(0),
    transformationSettings = NULL
  )
)

simulationSettings <- list(
  databaseSettings        = databaseSettings,
  propensitySettings      = propensitySettings,
  baselineRiskSettings    = baselineRiskSettings,
  treatmentEffectSettings = treatmentEffectSettings
)

analysisSettings <- SimulationEvaluationHte::createAnalysisSettings(
  threads        = 4,
  seed           = 19910930,
  replications   = 4,
  validationSize = 1e3,
  analysisId     = paste(
    "scenario",
    idSettings$scenario,
    sep = "_"
  ),
  description = "description",
  saveDirectory = "data/raw"
)

smoothSettings <- list(
  constant = list(
    type = "modelBased",
    label = "constant_treatment_effect",
    settings = createModelBasedSettings(
      type = "treatment",
      model = "logistic",
      adjustmentCovariates = paste0("x", 1:8)
    )
  ),
  stratified = createSmoothSettings(
    type = "stratified",
    settings = SmoothHte::createStratifiedSettings(),
    label = "stratified"
  ),
  constantLp = list(
    type = "modelBased",
    label = "linear_predictor",
    settings = createModelBasedSettings()
  ),
  loess = createSmoothSettings(
    type = "loess",
    settings = SmoothHte::createLoessSettings(),
    label = "loess"
  ),
  rcs = createSmoothSettings(
    type = "rcs",
    settings = SmoothHte::createRcsSettings(),
    label = "rcs"
  ),
  locfit = createSmoothSettings(
    type = "locfit",
    settings = SmoothHte::createLocfitSettings(),
    label = "locfit"
  )
)

predictionSettings <- createPredictionSettings(
  args = list(
    formula = "outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + treatment",
    family = "binomial"
  ),
  fun = "glm"
)

dat <- SimulateHte::runDataGeneration(
  databaseSettings,
  propensitySettings,
  baselineRiskSettings,
  treatmentEffectSettings
) %>%
  dplyr::tibble()
predictionSettings$args$data <- dat

predictionModel <- do.call(
  eval(parse(text = predictionSettings$fun)),
  predictionSettings$args
)

dat <- dat %>%
  dplyr::mutate(
    riskLinearPredictor = predict(
      predictionModel,
      newdata = dat %>%
        dplyr::mutate(treatment = 0)
    )
  )

simulatedDataset0 <- dat %>%
  dplyr::filter(treatment == 0)
simulatedDataset1 <- dat %>%
  dplyr::filter(treatment == 1)

# ==============================================================================
# Locfit
# ==============================================================================
s0 <- SmoothHte::fitLocfitHte(
  data = simulatedDataset0,
  settings = smoothSettings$locfit$settings
)
s1 <- SmoothHte::fitLocfitHte(
  data = simulatedDataset1,
  settings = smoothSettings$locfit$settings)

predictLocfit <- SmoothHte::predictSmoothBenefit(
  p               = plogis(dat$riskLinearPredictor),
  smoothControl   = s0,
  smoothTreatment = s1
)

# ==============================================================================
# Restricted cubic splines
# ==============================================================================
s0 <- SmoothHte::fitRcsHte(
  data = simulatedDataset0,
  settings = smoothSettings$rcs$settings
)
s1 <- SmoothHte::fitRcsHte(
  data = simulatedDataset1,
  settings = smoothSettings$rcs$settings)

predictRcs <- SmoothHte::predictSmoothBenefit(
  p               = plogis(dat$riskLinearPredictor),
  smoothControl   = s0,
  smoothTreatment = s1
)
# ==============================================================================
# Loess
# ==============================================================================
s0 <- SmoothHte::fitLoessHte(
  data = simulatedDataset0,
  settings = smoothSettings$loess$settings
)
s1 <- SmoothHte::fitLoessHte(
  data = simulatedDataset1,
  settings = smoothSettings$loess$settings)

predictLoess <- SmoothHte::predictSmoothBenefit(
  p               = plogis(dat$riskLinearPredictor),
  smoothControl   = s0,
  smoothTreatment = s1
)


# ==============================================================================
# Risk stratification
# ==============================================================================
stratifiedHte <- SmoothHte::fitStratifiedHte(
  data = dat,
  settings = smoothSettings$stratified$settings
)

stratifiedHte$data <- stratifiedHte$data %>%
  mutate(
    estimate = 100 * estimate,
    lower    = 100 * lower,
    upper    = 100 * upper,
    meanRisk = 100 * meanRisk
  )

predictStratified <- SmoothHte::predictStratifiedBenefit(
  p             = plogis(dat$riskLinearPredictor),
  stratifiedHte = stratifiedHte
)

# ==============================================================================
# Plots
# ==============================================================================
res <- tibble(
  risk = 100 * plogis(dat$riskLinearPredictor),
  predictLoess      = 100 * predictLoess,
  predictLocfit     = 100 * predictLocfit,
  predictRcs        = 100 * predictRcs,
  predictStratified = 100 * predictStratified
)

dataRect <- function(xmin, xmax) {
  data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf
  )
}

ggplot() +
  geom_rect(
    data = dataRect(-Inf, 3.205238),
    inherit.aes = FALSE,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "#f0f9e8",
    alpha = .2
  ) +
  geom_rect(
    data = dataRect(3.205238, 10.302985),
    inherit.aes = F,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "#bae4bc",
    alpha = .2
  ) +
  geom_rect(
    data = dataRect(10.302985, 29.368621),
    inherit.aes = F,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "#7bccc4",
    alpha = .2
  ) +
  geom_rect(
    data = dataRect(29.368621, Inf),
    inherit.aes = F,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "#2b8cbe",
    alpha = .2
  ) +
  geom_pointrange(
    data = stratifiedHte$data,
    aes(
      x = meanRisk,
      y = estimate,
      ymin = lower,
      ymax = upper,
      color = "Risk stratification"
    )
  ) +
  geom_line(
    data = res,
    aes(x = risk, y = predictLoess, color = "Loess"),
    size = .8
  ) +
  geom_line(
    data = res,
    aes(x = risk, y = predictRcs, color = "Restricted cubic spline"),
    size = .8
  ) +
  geom_line(
    data = res,
    aes(x = risk, y = predictLocfit, color = "Locfit"),
    size = .8
  ) +
  scale_color_manual(
    values = c("black", "#66c2a5", "#fc8d62", "#8da0cb"),
    breaks = c("Risk stratification", "Loess", "Restricted cubic spline", "Locfit")
  ) +
  xlab("Predicted risk (%)") +
  ylab("Predicted benefit (%)") +
  xlim(c(0, 60)) +
  ylim(c(-.5, 20)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  ggsave("smooth.png")
