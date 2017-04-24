#include "hphase.h"

index_db create_name_index(const SEXP names);

void do_impute(const SEXP Obs, const int nrow, const int *female,
	       const int *subset, int nuse, index_db snp_names, SEXP Rule, 
	       GTYPE **gt2ht, double *value_a, double *value_d);
