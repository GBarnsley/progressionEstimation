# Load dependencies
require(tidyverse)
require(magrittr)
require(rstan)
require(bridgesampling)
require(loo)
require(cowplot)
require(ggrepel)
require(xlsx)

# Load stan package
devtools::load_all()
pkgbuild::compile_dll()

set.seed(1234)

studies_by_size <-
  S_pneumoniae_infant_serotype %>%
    group_by(study) %>%
    dplyr::mutate(total_count = sum(carriage)+sum(disease)) %>%
    dplyr::select(study,total_count) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(total_count)) %>%
    dplyr::select(study) %>%
    dplyr::pull()

S_pneumoniae_infant_serotype %<>%
  dplyr::mutate(study = factor(study, levels = studies_by_size))


serotype_by_count <-
  S_pneumoniae_infant_serotype %>%
    group_by(type) %>%
    dplyr::mutate(total_count = sum(carriage)+sum(disease)) %>%
    dplyr::select(type,total_count) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(total_count)) %>%
    dplyr::select(type) %>%
    dplyr::pull()

S_pneumoniae_infant_serotype %<>%
  dplyr::mutate(type = factor(type, levels = serotype_by_count))

threshold_count <- 5

filtered_S_pneumoniae_infant_serotype <-
  S_pneumoniae_infant_serotype %>%
    dplyr::filter(carriage >= threshold_count | disease >= threshold_count)
# %>%
#     dplyr::group_by(type)%>%
#     dplyr::mutate(study_count = n()) %>%
#     dplyr::ungroup() %>%
#     dplyr::filter(study_count > threshold_count) %>%
#     dplyr::group_by(study) %>%
#     dplyr::mutate(serotype_count = n()) %>%
#     dplyr::ungroup() %>%
#     dplyr::filter(serotype_count > threshold_count) %>%
#     dplyr::mutate(type = factor(type))

filtered_S_pneumoniae_infant_serotype %<>%
    dplyr::mutate(study = factor(study, levels = studies_by_size[studies_by_size %in% unique(filtered_S_pneumoniae_infant_serotype$study)]))

filtered_S_pneumoniae_infant_serotype %<>%
    dplyr::mutate(type = factor(type,
                                          levels = serotype_by_count[serotype_by_count %in% unique(filtered_S_pneumoniae_infant_serotype$type)]))

infant_serotype_invasiveness <- 
  process_input_data(filtered_S_pneumoniae_infant_serotype)

n_chains <- 2
n_iter <- 2.5e4
n_core <- 2

full_infant_serotype_invasiveness <-
  process_input_data(
    S_pneumoniae_infant_serotype %<>%
      dplyr::mutate(study = factor(study, levels = studies_by_size)
      )
    )

full_adjusted_type_specific_negbin_fit <- 
  fit_progression_rate_model(full_infant_serotype_invasiveness,
                                                    type_specific = TRUE,
                                                    location_adjustment = TRUE,
                                                    stat_model = "negbin",
                                                    num_chains = n_chains,
                                                    num_iter = n_iter,
                                                    num_cores = n_core,
                                                    seed = 1200)

rhat_plot <- plot(full_adjusted_type_specific_negbin_fit, plotfun = "rhat", binwidth = 0.00005)

trace_plot <- rstan::traceplot(full_adjusted_type_specific_negbin_fit, pars = "lp__") +
                theme(axis.text.x = element_text(angle = 90))

best_fitting_model_validation_plot <- cowplot::plot_grid(plotlist = list(rhat_plot,trace_plot),
                   nrow = 1,
                   ncol = 2,
                   labels = "AUTO")

ggsave(best_fitting_model_validation_plot,
       file = "infant_model_validation_plot.pdf",
       width = 8,
       height = 6)

#extract posterior samples

# Extract factor levels
j_levels <- levels(S_pneumoniae_infant_serotype$type)

# Calculate invasiveness values
nu_name <- "nu"
if ("nu_j" %in% full_adjusted_type_specific_negbin_fit@model_pars) {
  nu_name <- "nu_j"
}

invasiveness <- rstan::extract(full_adjusted_type_specific_negbin_fit, pars=c(nu_name))[[nu_name]]
colnames(invasiveness) <- j_levels
write_csv(as_tibble(invasiveness), "child_invasiveness_posterior.csv")

#for comparison to the original output
pivot_longer(as_tibble(invasiveness), cols = everything()) |>
summarise(
  mean = median(value),
  lower = quantile(value, 0.025),
  upper = quantile(value, 0.975),
  .by = name
) |>
print(n=Inf)
