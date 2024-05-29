DATA_SECTION

!! ad_comm::change_datafile_name("log_r.dat");
  init_number inlog_r;

// !! ad_comm::change_datafile_name("log_K.dat");
//   init_number inlog_K;

// !! ad_comm::change_datafile_name("log_q.dat");
//   init_number inlog_q;
  
// !! ad_comm::change_datafile_name("log_sd_cpue.dat");
//   init_number inlog_sd_cpue;

!! ad_comm::change_datafile_name("surp_prod1.dat");
 init_int fyear;
 init_int lyear;
 init_vector cat(fyear,lyear);
 init_vector cpue(fyear,lyear);


PARAMETER_SECTION

//  init_number log_r;
 init_number log_q;
 init_number log_K;
 init_number log_sd_cpue;
 init_vector log_F(fyear,lyear);
 
 number log_r; // fixed parameter here
 number r;
 number q;
 number K;
 number sd_cat;
 number sd_cpue;

 vector bio(fyear,lyear+1);
 vector cat_hat(fyear,lyear);
 vector cpue_hat(fyear,lyear);
 vector expl_out(fyear,lyear);
 vector F(fyear,lyear);
 
 objective_function_value jnll;

!! log_r = inlog_r;
// !! log_K = inlog_K;
// !! log_q = inlog_q;
// !! log_sd_cpue = inlog_sd_cpue;


INITIALIZATION_SECTION

//  log_r -0.6
 log_q -9
 log_K 8.5 
 log_F 1

 
PROCEDURE_SECTION

 int t;
 dvariable expl;
 dvariable sum_sq;
 
 // Convert from log to normal space
 r = mfexp(log_r);
 q = mfexp(log_q);
 K = mfexp(log_K);
 F = mfexp(log_F);
 sd_cat = 0.05;
 sd_cpue = mfexp(log_sd_cpue);

 // Project the model forward
 bio(fyear) = K;
 for (t=fyear; t<=lyear; t++) {
   expl = 1.0/(1.0+F(t));
   bio(t+1) = bio(t) + r*bio(t)*(1.0-bio(t)/K) - expl*bio(t);
   cat_hat(t) = expl * bio(t);
   expl_out(t) = expl;
   cpue_hat(t) = q * bio(t);
  } 
 

 // Compute the likelihoods  
 jnll = 0;
 for (t=fyear; t<=lyear; t++) {
  jnll += 0.5 * square(log(cat(t)/cat_hat(t)) / sd_cat) + log(sd_cat);
  jnll += 0.5 * square(log(cpue(t)/cpue_hat(t)) / sd_cpue) + log(sd_cpue);
 }
 

GLOBALS_SECTION

  #include <admodel.h>
  #include <admb2r.cpp>


REPORT_SECTION

 open_r_file("out.rdat", 6, -999);
  wrt_r_complete_vector("obs_cat", cat);
  wrt_r_complete_vector("obs_cpue", cpue);
  wrt_r_complete_vector("est_bio", bio);
  wrt_r_complete_vector("est_cat", cat_hat);
  wrt_r_complete_vector("est_cpue", cpue_hat);
  wrt_r_complete_vector("est_expl", expl_out);
  wrt_r_item("jnll", jnll);
 close_r_file();


FINAL_SECTION

  ofstream myout("estpars.dat",ios::app);
    myout<< inlog_r << " " << log_r << " " << log_q << " " << log_K << " " << log_sd_cpue << " " << jnll << endl;
  myout.close();