// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------
// Utilities.
// ----------------------------------------------------------------------------

// Add Leading Value
//
// @param x Input vector.
// @param value Value to insert at leading position.
arma::colvec AddLeadVal(const arma::colvec &x, const double value) {
  const int len = x.n_elem;
  arma::colvec out = arma::zeros(len + 1);
  out(0) = value;
  out.subvec(1, len) = x;
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
//' @param time Event time.
//' @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
//' @return Data.frame with the censorings, deaths, and events occurring
//'   at each distinct time point. 
//' @export
// [[Rcpp::export]]
SEXP TabulateEvents(
  const arma::colvec time,
  const arma::colvec status  
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
// @param time Event time.
// @param status Status, coded as 0 for censoring, 1 for an event, 2 for death.
// @return Data.frame with the censorings, deaths, and events occurring
//   at each distinct time point. 
CICTab TabulateEventsCpp(
  const arma::colvec time,
  const arma::colvec status  
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
   
   // Output structure.
   CICTab out;
   out.time = unique_times;
   out.censor = n_cens;
   out.event = n_event;
   out.death = n_death;
   out.nar = nar;
   return out;
 };


// ----------------------------------------------------------------------------

//' Calculate CIC
//' 
//' @param time Observation time.
//' @param status Event status. The event coded as 1 is assumed to be the event
//'   of interest.
//' @return Tabulate cumulative incidence curve. 
//' @export 
// [[Rcpp::export]]
SEXP CalcCIC(const arma::vec &time, const arma::vec &status) {
  // Tabulate events and numbers at risk.
  int n = time.n_elem;
  CICTab out = TabulateEventsCpp(time, status);
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
  out.var_cic_event = n * (var1 + var2);
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

