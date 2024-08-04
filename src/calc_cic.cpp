// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------
// Utilities.
// ----------------------------------------------------------------------------


// Check if Value is in Vector
// 
// @param a Value to search for.
// @param b Vector to search.
// @return bool.
bool IsIn(const double &a, const arma::vec &b) {
  
  for(int i=0; i<b.size(); i++) {
    if(b(i) == a) {
      return true;
    }
  }
  return false;
};


// Union
// 
// @param a Vector.
// @param b Vector.
// @return Vector.
arma::colvec Union(const arma::colvec &a, arma::colvec b) {
  
  // Sets all elements of b that are in a to a value that is known
  // to be in a. Then concatenates a and b and filters to unique values.
  
  double a0 = a(0);
  for(int i=0; i<b.size(); i++) {
    if(IsIn(b(i), a)){
      b(i) = a0;
    }
  }

  arma::colvec out = arma::unique(arma::join_cols(a, b));
  return out;
};


// Truncate
//
// @param time Vector of time points.
// @param tau Truncation time.
// @return Truncated vector.
arma::colvec Truncate(const arma::colvec &time, const double tau) {
  arma::colvec unique_times = arma::unique(time);
  for (int i=0; i<unique_times.size(); i++) {
    if (unique_times(i) > tau) {
      unique_times(i) = tau;
    }
  }

  // Add truncation time, if not already present.
  arma::colvec out;
  if (!arma::any(unique_times == tau)) {
    int n = unique_times.n_elem;
    out = arma::zeros(n + 1);
    out.subvec(0, n - 1) = unique_times;
    out(n) = tau;
  } else {
    out = unique_times;
  }

  return arma::unique(out);
};


// Add Leading Value
//
// @param x Input vector.
// @param value Value to insert at leading position.
arma::colvec AddLeadVal(const arma::colvec &x, const double value) {
  const int n = x.n_elem;
  arma::colvec out;
  if (!arma::any(x == value)) {
    out = arma::zeros(n + 1);
    out(0) = value;
    out.subvec(1, n) = x;
  } else {
    out = x;
  }
  return out;
};


// ----------------------------------------------------------------------------
// Main functions.
// ----------------------------------------------------------------------------

//' Tabulate Events
//' 
//' Tabulate the number at risk and the number of events at each unique
//' observation time.
//' 
//' @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
//' @param time Event time.
//' @return Data.frame with the censorings, deaths, and events occurring
//'   at each distinct time point. 
//' @export
// [[Rcpp::export]]
SEXP TabulateEvents(
  const arma::colvec &status,
  const arma::colvec &time
){

  // Initial number at risk.
  int n = time.n_elem;

  // Unique observation times.
  arma::vec unique_times = arma::unique(time);
  unique_times = AddLeadVal(unique_times, 0);
  int n_time = unique_times.n_elem;

  // Initialize output vectors.
  arma::colvec n_cens  = arma::zeros(n_time);
  arma::colvec n_event = arma::zeros(n_time);
  arma::colvec n_death = arma::zeros(n_time);
  
  arma::colvec nar = arma::zeros(n_time);
  arma::colvec par = arma::zeros(n_time);

  // Set current number at risk.
  int current_nar = n;
  for(int i=0; i<n_time; i++) {
    
    double current_time = unique_times(i);
    arma::colvec current_status = status.elem(arma::find(time == current_time));

    // Counts.
    nar(i) = current_nar; // NAR at the beginning of the time step.
    par(i) = nar(i) / n;
    n_cens(i)  = arma::accu(current_status == 0.0);
    n_event(i) = arma::accu(current_status == 1.0);
    n_death(i) = arma::accu(current_status == 2.0);

    // Update number at risk for next interval.
    current_nar -= n_cens(i) + n_event(i) + n_death(i);

  };

  return Rcpp::DataFrame::create(
    Rcpp::Named("time")=unique_times,
    Rcpp::Named("censor")=n_cens,
    Rcpp::Named("event")=n_event,
    Rcpp::Named("death")=n_death,
    Rcpp::Named("nar")=nar
  );
};


// Structure to hold CIC.
struct CICTab {
  arma::vec time;
  arma::vec censor;
  arma::vec event;
  arma::vec death;
  arma::vec nar;
  arma::vec death_rate;
  arma::vec event_rate;
  arma::vec haz;
  arma::vec surv_init;
  arma::vec cic_event;
  arma::vec cic_death;
  arma::vec var_cic_event;
  arma::vec se_cic_event;
};


// Tabulate Events Cpp
// 
// Tabulate the number at risk and the number of events at each unique
// observation time.
// 
// @param eval_times Unique times at which to tabulate events.
// @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
// @param time Observation time.
// @return Data.frame with the censorings, deaths, and events occurring
//   at each distinct time point. 
CICTab TabulateEventsCpp(
  const arma::colvec &eval_times,
  const arma::colvec &status,
  const arma::colvec &time
){

  // Initial number at risk.
  int n = time.n_elem;

  // Unique observation times.
  arma::colvec unique_times = Union(eval_times, arma::unique(time));
  unique_times = AddLeadVal(unique_times, 0);
  int n_time = unique_times.n_elem;
   
  // Initialize output vectors.
  arma::colvec n_cens  = arma::zeros(n_time);
  arma::colvec n_event = arma::zeros(n_time);
  arma::colvec n_death = arma::zeros(n_time);
  arma::colvec nar = arma::zeros(n_time);
   
  // Set current number at risk.
  int current_nar = n;
  for(int i=0; i<n_time; i++) {
     
    double current_time = unique_times(i);
    arma::colvec current_status = status.elem(arma::find(time == current_time));
     
    // Counts.
    nar(i) = current_nar; // NAR at the beginning of the time step.
    n_cens(i)  = arma::accu(current_status == 0.0);
    n_event(i) = arma::accu(current_status == 1.0);
    n_death(i) = arma::accu(current_status == 2.0);
     
    // Update number at risk for next interval.
    current_nar -= n_cens(i) + n_event(i) + n_death(i);
     
  };

  // Restrict to evaluation times.
  int n_eval = eval_times.size();
  arma::colvec censor_out(n_eval);
  arma::colvec event_out(n_eval);
  arma::colvec death_out(n_eval);
  arma::colvec nar_out(n_eval);

  int j = 0;
  for(int i=0; i<n_time; i++) {
    
    double current_time = unique_times(i);
    if (IsIn(current_time, eval_times)) {
      censor_out(j) = n_cens(i);
      event_out(j) = n_event(i);
      death_out(j) = n_death(i);
      nar_out(j) = nar(i);
      j += 1;
    };

  };
   
  // Output structure.
  CICTab out;
  out.time = eval_times;
  out.censor = censor_out;
  out.event = event_out;
  out.death = death_out;
  out.nar = nar_out;
  return out;
 };


// ----------------------------------------------------------------------------

//' Calculate CIC
//' 
//' @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
//' @param time Observation time.
//' @return Tabulate cumulative incidence curve. 
//' @export 
// [[Rcpp::export]]
SEXP CalcCIC(const arma::vec &status, const arma::vec &time) {
  // Tabulate events and numbers at risk.
  // int n = time.n_elem;
  arma::colvec eval_times = AddLeadVal(arma::unique(time), 0);
  CICTab out = TabulateEventsCpp(eval_times, status, time);
  int n_time = out.time.n_elem;
  out.death_rate = out.death / out.nar;
  out.event_rate = out.event / out.nar;
  out.haz = (out.death + out.event) / out.nar;
  
  // Survival at the beginning of the interval.
  out.surv_init.set_size(n_time);
  out.surv_init(0) = 1.0;
  arma::vec cumprod_haz = arma::cumprod(1.0 - out.haz);
  out.surv_init.subvec(1, n_time - 1) = cumprod_haz.subvec(0, n_time - 2);
  
  // Cumulative incidence curves.
  out.cic_event = arma::cumsum(out.surv_init % out.event_rate);
  out.cic_death = arma::cumsum(out.surv_init % out.death_rate);
  
  // Variance calculation.
  // See equation (3) of <https://pubmed.ncbi.nlm.nih.gov/9160487/>.
  arma::vec var1 = arma::square(out.cic_event) % arma::cumsum(out.event / arma::square(out.nar)) +
    arma::cumsum(arma::square(1.0 - out.cic_death) % out.event / arma::square(out.nar)) -
    2.0 * out.cic_event % arma::cumsum((1.0 - out.cic_death) % out.event / arma::square(out.nar));
  
  arma::vec var2 = arma::square(out.cic_event) % arma::cumsum(out.death / arma::square(out.nar)) +
    arma::cumsum(arma::square(out.cic_event) % out.death / arma::square(out.nar)) -
    2.0 * out.cic_event % arma::cumsum(out.cic_event % out.death / arma::square(out.nar));
  
  // Output.
  out.var_cic_event = (var1 + var2);
  out.se_cic_event = arma::sqrt(var1 + var2);
  return Rcpp::DataFrame::create(
    Rcpp::Named("time")=out.time,
    Rcpp::Named("nar")=out.nar,
    Rcpp::Named("censor")=out.censor,
    Rcpp::Named("event")=out.event,
    Rcpp::Named("death")=out.death,
    Rcpp::Named("rate_event")=out.event_rate,
    Rcpp::Named("rate_death")=out.death_rate,
    Rcpp::Named("haz_total")=out.haz,
    Rcpp::Named("surv_init")=out.surv_init,
    Rcpp::Named("cic_event")=out.cic_event,
    Rcpp::Named("cic_death")=out.cic_death,
    Rcpp::Named("var_cic_event")=out.var_cic_event,
    Rcpp::Named("se_cic_event")=out.se_cic_event
  );
};


// ----------------------------------------------------------------------------
// Martingales.
// ----------------------------------------------------------------------------


// Calculate Martingales Cpp
//
// Construct a subject (row) by evaluation time (col) matrix where
// dMj[i, t] = dNj[i, t] - Y[i, t]dAj[i, t].
// 
// @param code Status value for the event of interest.
// @param cshaz Cause-specific hazard at each unique time.
// @param eval_times Unique times at which to evaluate the martingale.
// @param status Subject status.
// @param time Subject observation times.
// @return Matrix with subjects as rows and unique times as columns.
arma::mat CalcMartingale(
    const int code,
    const arma::colvec &cshaz,
    const arma::colvec &eval_times,    
    const arma::colvec &status,
    const arma::colvec &time
) {
  
  // Subjects.
  const int n = time.size();
  
  // Unique times.
  const int n_times = eval_times.size();
  
  // Create a subject by evaluation times matrix, where
  // dM[i, t] is the martingale increment for subject i at time t.
  arma::mat dm = arma::zeros(n, n_times);
  
  // Loop over subjects.
  for(int i=0; i<n; i++) {
    
    // Time and status for the focus subject.
    const double subj_time = time(i);
    const int subj_status = status(i);

    // Loop over times.
    for(int j=0; j<n_times; j++) {
      
      const double current_time = eval_times(j);
      const double current_haz = cshaz(j);
      
      // Add dN_{ji}(t).
      if(current_time == subj_time && subj_status == code) {
        dm(i, j) += 1;
      }
      
      // Add -Y_{i}(t)dAj(t).
      if(current_time <= subj_time) {
        dm(i, j) -= current_haz;
      } else {
        break;  
      }

    } // End loop over times.
  } // End loop over subjects.
  return dm;
}


// ----------------------------------------------------------------------------

//' Influence CIC
//' 
//' @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
//' @param time Observation time.
//' @param trunc_time Time at which to evaluate the influence function.
//' @return Tabulate cumulative incidence curve. 
//' @export 
// [[Rcpp::export]]
SEXP InfluenceCIC(
  const arma::colvec &status,
  const arma::colvec &time,
  const double trunc_time 
){

  // Evaluation times.
  arma::colvec eval_times = Truncate(arma::unique(time), trunc_time);

  // Tabulate events and numbers at risk.
  int n = time.n_elem;
  CICTab tab = TabulateEventsCpp(eval_times, status, time);

  int n_time = tab.time.n_elem;
  tab.death_rate = tab.death / tab.nar;
  tab.event_rate = tab.event / tab.nar;
  tab.haz = (tab.death + tab.event) / tab.nar;
  
  // Survival at the beginning of the interval.
  tab.surv_init.set_size(n_time);
  tab.surv_init(0) = 1.0;
  arma::vec cumprod_haz = arma::cumprod(1.0 - tab.haz);
  tab.surv_init.subvec(1, n_time - 1) = cumprod_haz.subvec(0, n_time - 2);
  
  // Cumulative incidence curves.
  tab.cic_event = arma::cumsum(tab.surv_init % tab.event_rate);
  tab.cic_death = arma::cumsum(tab.surv_init % tab.death_rate);

  // Calculate martingales.
  arma::mat dM1 = CalcMartingale(
    1, tab.event_rate, tab.time, status, time);
  arma::mat dM2 = CalcMartingale(
    2, tab.death_rate, tab.time, status, time);
  arma::mat dM = dM1 + dM2;

  // Cumulative incidence at truncation time.
  double ft = arma::as_scalar(tab.cic_event.elem(arma::find(eval_times == trunc_time)));

  // Calculate influence function.
  arma::colvec dMi;
  arma::colvec dM1i;
  arma::colvec influence = arma::zeros(n);

  for(int i=0; i<n; i++) {
    dMi = arma::trans(dM.row(i));
    dM1i = arma::trans(dM1.row(i));

    // Term 1.
    double t1 = -ft * n * arma::accu(dMi / tab.nar);

    // Term 2.
    double t2 = n * arma::accu(tab.cic_event % dMi / tab.nar);

    // Term 3.
    double t3 = n * arma::accu(tab.surv_init % dM1i / tab.nar);

    influence(i) = t1 + t2 + t3;
  }

  return Rcpp::wrap(influence);
};

