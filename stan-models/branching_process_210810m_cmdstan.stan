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

    matrix calculate_bin_lpmf_matrix(
        //data
        int max_N_cs_actual,
        int N_cs_obs,
        vector bin_successes,
        //pars
        real rho
    )
    {
        // declare
        matrix[N_cs_obs+1, max_N_cs_actual+1] bin_lpmf;
        real log_rho;
        real log_one_minus_rho;
        int tmp_int;
            
        // define
        log_rho = log( rho );
        log_one_minus_rho = log( 1 - rho );

        // calculate lpmf of sampling probabilities given rho
        // entry i,j is the probability of observed chain of size i having actual size j
        for(j in 1:(max_N_cs_actual+1))
        {
            tmp_int = j;
            if(j > N_cs_obs+1)
            {
                tmp_int = N_cs_obs+1;
            }
            bin_lpmf[1:tmp_int,j] = bin_successes[1:tmp_int] * log_rho;
            bin_lpmf[1:tmp_int,j] += ( (rep_vector(j-1, tmp_int) - bin_successes[1:tmp_int]) * log_one_minus_rho );
            for(i in 1:tmp_int)
            {
                bin_lpmf[i,j] += lchoose( 1.*(j-1), bin_successes[i] );
            }
            if(tmp_int<(N_cs_obs+1))
            {
                bin_lpmf[(tmp_int+1):(N_cs_obs+1), j ] = rep_vector( negative_infinity(), N_cs_obs+1-tmp_int);
            }
        }
        return(bin_lpmf);
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
		    
        vector[N_jcases_sbt+1] ones_n = rep_vector(1., N_jcases_sbt+1);
		
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

    real log_dens_for_given_index(
        int[,] index_cases,
       	int start,
        int end,
        // data
        int N_sbts,
        int max_N_cs_actual,
        int[] N_jcases,
        int[,,] cs_obs_pre,
        int[,,] cs_obs_post,
        int renormalise,
        // parameters
        vector log_r0,
        vector vmr_minus_one,
        //transformed parameters
        matrix bin_lpmf
    )
    {
        real lpmf = 0.0;
        int C_slice = end - start + 1;
        vector[N_sbts] r0 = exp( log_r0 );
    
        for(c_slice in 1:C_slice)
        {
            int c = c_slice + start - 1;
          
            for(s in 1:N_sbts)
            {
                vector[max_N_cs_actual] ones = rep_vector(1., max_N_cs_actual);
                row_vector[max_N_cs_actual] cs_actual_lpmf;
                row_vector[max_N_cs_actual] cs_obs_pre_lpmf;
                row_vector[max_N_cs_actual] cs_obs_post_lpmf;
         
                // next if no index cases corresponding to c_slice for this subtype
                if(sum(cs_obs_pre[c,,s]) == 0) continue;
            
                // calculate lpmf of actual chain sizes given R0, kappa and initial cases (icases)
                cs_actual_lpmf =
                    lpmf_actual_cs_for_given_index(
                        max_N_cs_actual,
                        c,
                        r0[s],
                        vmr_minus_one[s]);
            
                // renormalise lpmf to ensure probs sum to 1 to due infinite sum approximation
                cs_actual_lpmf[1:max_N_cs_actual] -=  log( exp(cs_actual_lpmf[1:max_N_cs_actual]) * ones[1:max_N_cs_actual]);
						
                //	calculate lpmf of observed chain sizes given R0 and kappa and initial cases (icases)
                cs_obs_pre_lpmf =
                    lpmf_obs_cs(
                        max_N_cs_actual,
                        max_N_cs_actual,
                        cs_actual_lpmf,
                        bin_lpmf,
                        c,
                        0);

                //	calculate lpmf of observed chain sizes given R0 and kappa and initial cases (icases)
                cs_obs_post_lpmf =
                    lpmf_obs_cs_m1(
                        max_N_cs_actual,
                        max_N_cs_actual,
                        cs_actual_lpmf,
                        bin_lpmf,
                        c);


                lpmf += to_row_vector(cs_obs_pre[c,1:N_jcases[s],s]) * cs_obs_pre_lpmf[1:N_jcases[s] ]';
                if(c==1){
                        lpmf += to_row_vector(cs_obs_post[1,2:N_jcases[s],s]) * cs_obs_post_lpmf[2:N_jcases[s] ]';
                }
            }
        }
        return(lpmf);
    }
}

data{
    int<lower=1> N_sbts;                // total number of HIV subtypes
    int<lower=1> N_cs_obs;			    // max obs chain size
    int<lower=1> N_cs_actual[N_sbts];	// max size of actual chain
    int<lower=1> N_icases;				// max number of unique index cases in data
    int<lower=1> N_jcases[N_sbts];		// max number of generated cases in data
    int<lower=-1> index_cases[N_icases,N_sbts];		 // individuals on ART by 2010
    int<lower=0> cs_obs[N_icases+1,max(N_jcases),N_sbts]; // observed chain sizes, rows=index cases, cols=generated cases
    int<lower=1> sampling_n;						    // sampling total
    int<lower=1,upper=sampling_n> sampling_k;		    // sampling success
    int<lower=0> index_flag;			// whether index cases should be assumed part of observed subgraph
}

transformed data{
    int<lower=1> max_N_cs_actual = max(N_cs_actual);
    // separate chain size data into pre/post cut-off date for index cases.
    // cs_obs_post are those which started after the cut-date with 0 obs index cases
    int<lower=0> cs_obs_pre[N_icases,max(N_jcases),N_sbts];
    int<lower=0> cs_obs_post[1,max(N_jcases),N_sbts];
    int<lower=0> renormalise;
    vector[max_N_cs_actual+1] bin_successes;	

    // if index_cases not counted as part of subgraph set flag to renormalise pmf
    if(index_flag==0) {
        renormalise = 1;
    } else{
        renormalise = 0;
    }
	
    // subgraphs started pre-2010 just keep m>=1
    cs_obs_pre = cs_obs[2:N_icases+1,,];
	
    // subgraphs started post-cutoff, just keep m=0
    cs_obs_post[1:1,,] = cs_obs[1:1,,];
	
    //	get Binomial sampling successes 0,...,N_cs_obs.
    bin_successes[1] = 0;
    for(i in 1:max_N_cs_actual)
    {
        bin_successes[i+1] = bin_successes[i]+1; 
    }
}

parameters{
    real<lower=-23, upper=0> log_r0_overall;
    real<lower=0> log_r0_sd;
    vector<lower=-23, upper=0>[N_sbts] log_r0_sbts;
    real<lower=0> inv_vmr_minus_one_overall;
    vector<lower=0>[N_sbts] vmr_minus_one_sbts;
    real<lower=0, upper=1> rho;
}

transformed parameters
{
    matrix[max_N_cs_actual+1, max_N_cs_actual+1] bin_lpmf;
	
    bin_lpmf =
        calculate_bin_lpmf_matrix(
            max_N_cs_actual,
            max_N_cs_actual,
            bin_successes,
            rho
            );
}

model
{
    target+= normal_lpdf( log_r0_overall | -1.47, 0.5);
    target+= exponential_lpdf( log_r0_sd | 2);
    target+= normal_lpdf( log_r0_sbts | log_r0_overall, log_r0_sd );
                                        
    target+= exponential_lpdf( inv_vmr_minus_one_overall | 1);
    target+= exponential_lpdf( vmr_minus_one_sbts | inv_vmr_minus_one_overall );
    
    target+= beta_lpdf(rho | sampling_k+0.5, sampling_n-sampling_k+0.5);
	
    // rstan version
    // target += log_dens_for_given_index(index_cases, 1, N_icases,
    // cmdstan version
    target += reduce_sum(log_dens_for_given_index, index_cases, 1,
        N_sbts,
        max_N_cs_actual,
        N_jcases,
        cs_obs_pre,
        cs_obs_post,
        renormalise,
        log_r0_sbts,
        vmr_minus_one_sbts,
        bin_lpmf);
}

