data
{
    int<lower=0> n; // Max number of observations (MSM or HSX)
    int<lower=0> r; // Number of migrant groups
    int<lower=1> idx_to_obs_array[n]; // group index of each observation
    vector[n] y;  // time to diagnosis for all observations
    real log_wb_quantiles_priormean[2]; // priors for MSM/HSX models
    int N_diag[r];
}

transformed data
{
    real log_weibull_constants[2] = { log(log(2)), log(log(5)) };
}

parameters
{
    real<lower=0> wb_log_q50_sd;
    real<lower=0> wb_log_q80_q50_sd;
    
    real<lower=-5,upper=5> wb_log_q50_overall;
    real<lower=0> wb_log_q80_q50_overall;
    vector<lower=-5,upper=5>[r] wb_log_q50_grp;
    vector<lower=0>[r] wb_log_q80_q50_grp;
    
}

transformed parameters
{
    real<lower=0> wb_scale_overall;
    real<lower=0> wb_shape_overall;
    vector<lower=0>[r] wb_scale_grp;
    vector<lower=0>[r] wb_shape_grp;

    // define wb_scale from wb_log_q50_grp and wb_log_q80_grp
    // define wb_shape from wb_log_q50_grp and wb_log_q80_grp
    wb_shape_grp = 
        ( log_weibull_constants[2] - log_weibull_constants[1] ) ./ wb_log_q80_q50_grp;
    wb_scale_grp =
        exp( wb_log_q50_grp - log_weibull_constants[1] ./ wb_shape_grp );
    wb_shape_overall = 
        ( log_weibull_constants[2] - log_weibull_constants[1]) ./ wb_log_q80_q50_overall;
    wb_scale_overall =
        exp( wb_log_q50_overall - log_weibull_constants[1] ./ wb_shape_overall );

}

model
{
    // hierarchical prior on Weibull distribution parameter, for each migrant group
    target += normal_lpdf( wb_log_q50_overall | log_wb_quantiles_priormean[1], 0.5); 
    target += normal_lpdf( wb_log_q80_q50_overall | log_wb_quantiles_priormean[2], 0.5); 
    target += exponential_lpdf( wb_log_q50_sd | 2);
    target += lognormal_lpdf( wb_log_q80_q50_sd | 0, 1);
    target += normal_lpdf( wb_log_q50_grp |  wb_log_q50_overall, wb_log_q50_sd );
    target += normal_lpdf( wb_log_q80_q50_grp | wb_log_q80_q50_overall, wb_log_q80_q50_sd );
    
    // likelihood of time to diagnoses for each observation
    target += weibull_lpdf( y | wb_shape_grp[ idx_to_obs_array ], wb_scale_grp[ idx_to_obs_array ] );
}


 generated quantities
{
    
		real p_undiag_year_trunc[r,5]; 
		real p_undiag_av[r];
	  vector[r] undiagnosed;
	  
    for(i in 1:r){
	    for(t in 1:5){
				    p_undiag_year_trunc[i,t] = 1 - weibull_cdf( t , wb_shape_grp[i], wb_scale_grp[i]);
	    	}
			p_undiag_av[i] = mean(p_undiag_year_trunc[i,]);
			undiagnosed[i] = round(N_diag[i]/(1-p_undiag_av[i]));
    }
}

