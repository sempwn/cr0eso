
# cr0eso 0.0.9

* New `prior_list` list of parameterizations of priors that can be added to `seir_model_fit()`. 

* In `seir_model_fit()`, also takes optional `priors` variable (#1).

* New `hom_plot_incidence_by_location()` plots posterior predictive distribution of incidence with incidence data by location.

* New `hom_plot_zeta_by_location()` plots model marginal posterior of intervention strength zeta by location with predictive marginal distribution.

* New `hom_plot_r0_by_location()` plots model marginal posterior of r0 by location with predictive marginal distribution.

* New `hom_extract_posterior_draws()` simplifies extraction of posterior from rstan fit.

* New `seir_model_fit()` provides mechanism for fitting multiple outbreaks.
