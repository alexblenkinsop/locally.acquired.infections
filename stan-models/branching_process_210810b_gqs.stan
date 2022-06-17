functions{
  
    row_vector lpmf_actual_cs_for_given_index(
        // data
        int N_cs_actual,
        real index_cases_local,
        // parameters
        real r0,
        real vmr_minus_one
    )
    {
        row_vector[N_cs_actual] cs_actual_lpmf;
        
        real kappa = r0/vmr_minus_one;
        real log_r0_div_kappa = log( vmr_minus_one );
        real log_vmr = log( 1+vmr_minus_one );
        
        vector[6] ones = rep_vector(1., 6);
        matrix[N_cs_actual,6] tmp;
        
        for(i in 1:N_cs_actual)
        {
            // i-1 because index 1 is 0 cases generated
            // first term is normalising factor for prob(extinction)
            tmp[i,1] = (1.*(index_cases_local)) / (1.*(index_cases_local+(i-1)));
            tmp[i,2]= kappa*index_cases_local + kappa*(i-1) + (i-1) ;
            tmp[i,3]= kappa*index_cases_local + kappa*(i-1);
            tmp[i,4]= ((i-1))+1;
            tmp[i,5]= (i-1);
            tmp[i,6]= kappa*index_cases_local + kappa*(i-1) + (i-1);
        }
        tmp[,1] = log( tmp[,1] );
        tmp[,2] = lgamma( tmp[,2] );
        tmp[,3] = -lgamma( tmp[,3] );
        tmp[,4] = -lgamma( tmp[,4] );
        tmp[,5] *= log_r0_div_kappa;
        tmp[,6] *= -log_vmr;
        
        cs_actual_lpmf = (tmp * ones)';
            
        return(cs_actual_lpmf);
    }
    
   row_vector lpmf_obs_cs(
        // data
        int N_cs_actual,
        int N_jcases_sbt,
        // lpmf of actual chain
        row_vector cs_actual_lpmf,
        matrix bin_lpmf,
        // index
        int c_local,
        int renormalise
    )
    {
        // lpmf of observed chains, adjust actual chain sizes for sampling probability
        // col indices of cs_actual_lpmf represent generated cases
        // indices of bin_lpmf represent chain size so only consider chains > c_local (index cases)
        row_vector[N_jcases_sbt] cs_obs_lpmf;
		    
        // pr(observe 0 new cases)
        cs_obs_lpmf[1] = log_sum_exp( cs_actual_lpmf[1:(N_cs_actual)] + (bin_lpmf[1, 1:(N_cs_actual)]) );
        for(i in 1:(N_jcases_sbt-1))
        {
            // pr(observe i new cases)
            cs_obs_lpmf[i+1] = log_sum_exp( cs_actual_lpmf[(i+1):N_cs_actual] + (bin_lpmf[(i+1), (i+1):(N_cs_actual)]) );
        }
		
        if(renormalise==1){
            // renormalise conditional on 1 - prob nothing sampled, where cols of cs_obs_lpmf are generated cases
            cs_obs_lpmf[ 2:(N_jcases_sbt) ] -= log( 1-exp(cs_obs_lpmf[1]) );
        }
        return(cs_obs_lpmf);
    }
       row_vector lpmf_obs_cs_m1(
        // data
        int N_cs_actual,
        int N_jcases_sbt,
        // lpmf of actual chain
        row_vector cs_actual_lpmf,
        matrix bin_lpmf,
        // index
        int c_local
    )
    {
        // lpmf of observed chains, adjust actual chain sizes for sampling probability
        // col indices of cs_actual_lpmf represent generated cases
        // indices of bin_lpmf represent chain size so only consider chains > c_local (index cases)
        row_vector[N_jcases_sbt] cs_obs_lpmf;
		    
        // pr(observe 0 new cases)
        cs_obs_lpmf[1] = log_sum_exp( cs_actual_lpmf[1:N_cs_actual] + (bin_lpmf[1, 2:(N_cs_actual+1)]) );
        for(i in 1:N_jcases_sbt-1)
        {
            cs_obs_lpmf[i+1] = log_sum_exp( cs_actual_lpmf[i:N_cs_actual] + (bin_lpmf[(i+1), (i+1):(N_cs_actual+1)]) );
        }

        //    renormalise conditional on 1 - prob nothing sampled
        cs_obs_lpmf[ 2:N_jcases_sbt ] -= log( 1-exp(cs_obs_lpmf[1]) );

        return(cs_obs_lpmf);
    }
}


data
{
    int<lower=1> N_sbts;                    // total number of HIV subtypes
    int<lower=1> N_cs_obs;				    // max obs chain size
    int<lower=1> N_cs_actual[N_sbts];	    // max size of actual chain
    int<lower=1> N_icases;					// max number of unique index cases in data
    int<lower=1> N_jcases[N_sbts];		    // max number of generated cases in data
    int<lower=1> sampling_n;                            // sampling total
    int<lower=1,upper=sampling_n> sampling_k;           // sampling success
    int<lower=0> cs_obs[N_icases+1,max(N_jcases),N_sbts]; // observed chain sizes, rows=index cases, cols=generated cases
    int<lower=0> N_origins;					// number of geographic origins
    int<lower=0> N_bplace;					// number of birthplaces
    int<lower=0> N_sgs_m;					// number of obs subgraphs for index case m
    matrix<lower=0>[N_sbts,N_origins] pr_origins;		  // proportion of subgraphs from each origin
    matrix<lower=0>[N_sbts,N_bplace] pr_bplace;		  // proportion of subgraphs from each origin
    int<lower=1, upper=N_icases> INDEX_CASE_IDX;
    int<lower=0> index_flag;			// whether index cases should be assumed part of observed subgraph
    int<lower=1> z[max(N_cs_actual)];
}

transformed data
{
    int<lower=1> max_N_cs_actual = max(N_cs_actual);      
    int<lower=0> renormalise;
    //	get Binomial sampling successes 0,...,N_cs_obs.
    vector[N_cs_obs+1] bin_successes;
    
    // if index_cases not counted as part of subgraph set flag to renormalise pmf
    if(index_flag==0) {
        renormalise = 1;
    } else{
        renormalise = 0;
    }
    
    bin_successes[1] = 0;
    for(i in 1:N_cs_obs)
    {
        bin_successes[i+1] = bin_successes[i]+1;
    }
}

parameters
{
    real<lower=-23, upper=0> log_r0_overall;
    real<lower=0> log_r0_sd;
    vector<lower=-23, upper=0>[N_sbts] log_r0_sbts;
    real<lower=0> inv_vmr_minus_one_overall;
    vector<lower=0>[N_sbts] vmr_minus_one_sbts;
    real<lower=0, upper=1> rho;
    matrix[max_N_cs_actual+1, max_N_cs_actual+1] bin_lpmf;
}

generated quantities
{
    real<lower=0> N_sgs_pre_st[N_sbts];
    real<lower=0> N_sgs_post_st[N_sbts];
    int<lower=0> N_sgs_unobs[N_sbts];
    real<lower=0> N_sgs_e[N_sbts];
    int<lower=0> N_inf_x[N_sbts];
    int<lower=0> N_inf_e[N_sbts];
    int actual_cs[N_sbts,N_sgs_m];
    int actual_cs_unobs[N_sbts,1000];
    matrix[N_sbts,N_sgs_m] obs_cs_pre;
    matrix[N_sbts,1000] obs_cs_post;
    matrix[N_sbts,N_sgs_m] origins_subgraphs_x;
    matrix[N_sbts,1000] origins_subgraphs_e;
    vector[max_N_cs_actual] unsampled_ch;
    real<lower=0> pr_unsampled[N_sbts];
    int<lower=0> origins_ind_x[N_sbts,N_bplace];
    int<lower=0> origins_ind_e[N_sbts,N_bplace];

    // generate expected cases by age + expected deaths by age
    {
        int c = INDEX_CASE_IDX;
        matrix[N_sbts,max_N_cs_actual] cs_actual_lpmf;
        matrix[N_sbts,max_N_cs_actual] cs_actual_pmf;
        vector[max_N_cs_actual+1] ones = rep_vector(1., max_N_cs_actual+1);
        vector[N_sbts] r0 = exp( log_r0_sbts );
        matrix[N_sbts,max_N_cs_actual] cs_obs_pre_lpmf;
        matrix[N_sbts,max_N_cs_actual] cs_obs_post_lpmf;
        matrix[N_sbts,max_N_cs_actual] cs_obs_pre_pmf;
        matrix[N_sbts,max_N_cs_actual-1] cs_obs_post_pmf;
        //vector[max_N_cs_actual] unsampled_ch;
        //real pr_unsampled[N_sbts];

        for(s in 1:N_sbts)
        {
            // number of subgraphs to predict for pre-exisiting/emergent subgraphs
            N_sgs_pre_st[s] = sum(cs_obs[c+1,,s]);
            N_sgs_post_st[s] = sum(cs_obs[1,,s]);
          
            // calculate lpmf of actual chain sizes given R0, kappa and initial cases (icases)
            cs_actual_lpmf[s,] =
                lpmf_actual_cs_for_given_index(
                    max_N_cs_actual,
                    c,
                    r0[s],
                    vmr_minus_one_sbts[s]
                    );
          
            // renormalise lpmf
            cs_actual_lpmf[s,] -=  log( exp(cs_actual_lpmf[s,]) * ones[1:max_N_cs_actual]);

            cs_actual_pmf[s,] = exp(cs_actual_lpmf[s,]);

            // when m=1, predict number of unobserved emergent chains
            if(c==1){
            	  // prob a chain is unsampled
            	  for(i in 1:max_N_cs_actual){
            	  	unsampled_ch[i] = (1-rho)^i / i*1.0;
            	  }
                pr_unsampled[s] = unsampled_ch' * cs_actual_pmf[s,]';
                // number of unsampled chains
                if(N_sgs_post_st[s]>0){
                    N_sgs_unobs[s] = neg_binomial_rng(N_sgs_post_st[s],(1-pr_unsampled[s])/pr_unsampled[s]);
                }else{
                	  N_sgs_unobs[s] = 0;
                }
            }else{
                unsampled_ch = rep_vector(0.,max_N_cs_actual);
                pr_unsampled[s] = 0;
                N_sgs_unobs[s] = 0;
            }
            
            N_inf_x[s] = 0;
            for(i in 1:N_sgs_m)
            {
                if(i<=N_sgs_pre_st[s]){
                    actual_cs[s,i] =  categorical_rng(to_vector(cs_actual_pmf[s,]));
                    N_inf_x[s] = N_inf_x[s] + (actual_cs[s,i] - 1);
                }else{
                    actual_cs[s,i] =  0;
                }
            }

            N_sgs_e[s] = N_sgs_post_st[s]+N_sgs_unobs[s];

            // predict emergent subgraph sizes only when m=1
            N_inf_e[s] = 0;
            if(c==1){
                for(i in 1:1000)
                {
                    if(i<=N_sgs_e[s]){
                        actual_cs_unobs[s,i] =  categorical_rng(to_vector(cs_actual_pmf[s,]));
                        N_inf_e[s] = N_inf_e[s] + (actual_cs_unobs[s,i]);
                    }else{
                        actual_cs_unobs[s,i] =  0;
                    }
                }
            }
            else{
                actual_cs_unobs[s,] =  rep_array(0,1000);
            }

             //	calculate lpmf of observed chain sizes given R0 and kappa and initial cases (icases)
            cs_obs_pre_lpmf[s,1:(max_N_cs_actual)] =
                lpmf_obs_cs(
                    max_N_cs_actual,
                    max_N_cs_actual,
                    cs_actual_lpmf[s,],
                    bin_lpmf,
                    c,
                    0);
            // renormalise lpmf
            cs_obs_pre_lpmf[s,1:(max_N_cs_actual)] -=  log( exp(cs_obs_pre_lpmf[s,1:(max_N_cs_actual)]) * ones[1:(max_N_cs_actual)]);

						cs_obs_pre_pmf[s,] = exp(cs_obs_pre_lpmf[s,]);

            //	calculate lpmf of observed chain sizes given R0 and kappa and initial cases (icases)
            cs_obs_post_lpmf[s,1:(max_N_cs_actual)] =
                lpmf_obs_cs_m1(
                    max_N_cs_actual,
                    max_N_cs_actual,
                    cs_actual_lpmf[s,],
                    bin_lpmf,
                    c);
           cs_obs_post_lpmf[s,2:(max_N_cs_actual)] -=  log( exp(cs_obs_post_lpmf[s,2:(max_N_cs_actual)]) * ones[2:(max_N_cs_actual)]);

           cs_obs_post_pmf[s,] = exp(cs_obs_post_lpmf[s,2:(max_N_cs_actual)]);

            for(i in 1:N_sgs_m)
            {
              if(i<=N_sgs_pre_st[s]){
                obs_cs_pre[s,i] =  categorical_rng(to_vector(cs_obs_pre_pmf[s,]));
              }else{
                obs_cs_pre[s,i] =  positive_infinity();
              }
            }
            
            // predict emergent subgraph sizes only when m=1
            if(c==1){
              for(i in 1:1000)
              {
                if(i<=N_sgs_post_st[s]){
                  obs_cs_post[s,i] =  categorical_rng(to_vector(cs_obs_post_pmf[s,]));
                }else{
                  obs_cs_post[s,i] =  positive_infinity();
                }
              }
            }else{
                obs_cs_post[s,] =  rep_row_vector(positive_infinity(),1000);
            }

            for(i in 1:N_sgs_m)
            {
            	if(i<=N_sgs_pre_st[s]){
                origins_subgraphs_x[s,i] = categorical_rng(to_vector(pr_origins[s,]));
            	}else{
            		origins_subgraphs_x[s,i] =  positive_infinity();
            	}
            }
        
            for(i in 1:1000)
            {
                if(i<=N_sgs_e[s]){
                    origins_subgraphs_e[s,i] = categorical_rng(to_vector(pr_origins[s,]));
                }else{
                	  origins_subgraphs_e[s,i] = positive_infinity();
                }
            }
            if(N_inf_x[s]>0){
                origins_ind_x[s,] = multinomial_rng(pr_bplace[s,]',N_inf_x[s]);
            }else{
            	origins_ind_x[s,] =	rep_array(0,N_bplace);
            }
            if(N_inf_e[s]>0){
                origins_ind_e[s,] = multinomial_rng(pr_bplace[s,]',N_inf_e[s]);
            }else{
            	  origins_ind_e[s,] = rep_array(0,N_bplace);
            }
        }
    }
}

