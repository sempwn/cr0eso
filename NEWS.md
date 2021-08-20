# cr0eso 0.0.0.9002

* In `seir_model_fit()`, also takes optional `likelihood_type` variable.



# cr0eso 0.0.0.9001

* In `seir_model_fit()`, also takes optional `chains` and `fit_type` variable (#28).

* New `BC_OSABC_facilitydata` meta-data on facilities in BC (#24).

* New `BC_LTHC_outbreaks_100Imputs` data on facility outbreaks in BC (#24).

# cr0eso 0.0.0.9

* New `create_pub_tables()` create tibble of $R_0$, $\zeta$ and critical time for variable number of models (#21).

* New `hom_plot_counterfactual_by_location()` plot counterfactual and predictive distribution (#17).

* New `prior_list` list of parameterizations of priors that can be added to `seir_model_fit()`. 

* In `seir_model_fit()`, also takes optional `priors` variable (#1).

* New `hom_plot_incidence_by_location()` plots posterior predictive distribution of incidence with incidence data by location.

* New `hom_plot_zeta_by_location()` plots model marginal posterior of intervention strength zeta by location with predictive marginal distribution.

* New `hom_plot_r0_by_location()` plots model marginal posterior of r0 by location with predictive marginal distribution.

* New `hom_extract_posterior_draws()` simplifies extraction of posterior from rstan fit.

* New `seir_model_fit()` provides mechanism for fitting multiple outbreaks.
