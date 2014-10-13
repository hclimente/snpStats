/* 
   Moving window on a large covariance matrix 

   `size' is the size of a moving square "window" on a large covariance 
   matrix. The window starts at `start', so that the window has top
   left-hand corner at [start, start] and bottom right-hand corner at
   [start+size-1, start+size-1]. 

   `evaluate(i, j, ap)' is a function which evaluates the
   covariance between elements `i' and `j' in the large matrix. `ap' is
   a list of further arguments passed through the calling functions
   `get_row' and `get_diag'. Covariances are cached and only re-evaluated 
   if they are not in the cache.

   `new_window(size, start)' creates a new window object.
   `free_window(win)' destroys the window and returns storage to the pool
   `move_window(win, new_start)' moves the window to a new starting position
      (keeping what covariances can be preserved).
   `get_row(win, i, row, evaluate, ...)' gets that part of row `i' of the
      large covariance matrix which falls within the current window and
      stores it in `row'. If they are not already available, covariances 
      are evaluated (using the function `evaluate') and cached.
      Trailing arguments are passed to `evaluate' as a va_list. If row `i' 
      does not fall within the window, `NA_REAL' is returned for all row
      elements.
   `get_diag(win, diag, evaluate, ...)' returns the diagonal of the current
      window. Again, elements not already available are calculated and cached 
      using the `evaluate' function.

   Internally the lower triangle of the window is cached with rows recharged 
   on a circular basis in order to avoid moves of large number of covariances
   when a window is moved. The `start' row of the window will correspond with
   row `start_local' in the cache.

*/
 
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h> /* To define NA_REAL */
#include "covwin.h"

COV_WIN_PTR new_window(const int size, const int start) {
  COV_WIN_PTR res = (COV_WIN_PTR) Calloc(1, COV_WIN);
  res->size = size;
  res->start = start;
  res->start_local = 0;
  int ncov = (size*(size+1))/2;
  res->covariances = (double *) Calloc(ncov, double);
  for (int i=0; i<ncov; i++)
    res->covariances[i] = NA_REAL;
  return res;
}

void free_window(COV_WIN_PTR win) {
  Free(win->covariances);
  Free(win);
}

void move_window(COV_WIN_PTR win, const int new_start){
  int start_local = win->start_local;
  if (abs(new_start - win->start)>=win->size) {
    int ncov = (win->size*(win->size+1))/2;
    for (int i=0; i<ncov; i++)
      win->covariances[i] = NA_REAL;
    win->start = new_start;
    win->start_local = 0;
  }
  else if (new_start > win->start) {
    for(int i=win->start; i<new_start; i++) {
      for (int j=0, ij=start_local; j<win->size; j++) {
	win->covariances[ij] = NA_REAL;
	if (j<start_local)
	  ij += win->size-1-j;
	else
	  ij++;
      }
      start_local++;
      if (start_local==win->size)
	start_local = 0;
    }
    win->start = new_start;
    win->start_local = start_local;
  }
  else if (new_start < win->start) {
    for(int i=win->start; i>new_start; i--) { 
      start_local--;
      if (start_local<0)
	start_local = win->size-1;
      for (int j=0, ij=start_local; j<win->size; j++) {
	win->covariances[ij] = NA_REAL;
	if (j<start_local)
	  ij += win->size-1-j;
	else
	  ij++;
      }
    }
    win->start = new_start;
    win->start_local = start_local;
  }
  return;
}

void get_row(COV_WIN_PTR win, const int i, double *row, 
	     double(*evaluate)(int, int, va_list), ...) {
  va_list ap;
  va_start(ap, evaluate);
  if ((i<win->start) || (i>=win->start+win->size)) {
    for (int j=0; j<win->size; j++) 
      row[j] = NA_REAL;
    return;
  }
  int i_local = (i - win->start + win->start_local) % win->size;
  for (int j=0, ij=i_local, j_out=win->size-win->start_local; j<win->size; 
       j++, j_out++) {
    if (j_out==win->size)
      j_out = 0;
    double w = win->covariances[ij];
    if (ISNA(w)) {
      va_list aq;
      va_copy(aq, ap);
      row[j_out] = win->covariances[ij] = 
	(*evaluate)(i, win->start+j_out, aq);
      va_end(aq);
    }
    else {
      row[j_out] = w;
    }
    if (j<i_local)
      ij += win->size-j-1;
    else
      ij++;
  }
  va_end(ap);
}

void get_diag(COV_WIN_PTR win, double *diag, 
	      double(*evaluate)(int, int, va_list), ...) {
  va_list ap;
  va_start(ap, evaluate);
  for (int i=0, ij=0, i_out=win->size - win->start_local; i<win->size; 
       i++, i_out++){
    if (i_out==win->size)
      i_out = 0;
    double w = win->covariances[ij];
    if (ISNA(w)) {
      int ii =  win->start+i_out;
      va_list aq;
      va_copy(aq, ap);
      diag[i_out] = win->covariances[ij] = (*evaluate)(ii, ii, aq);
      va_end(aq);
    }
    else {
      diag[i_out] = w;
    }
    ij += (win->size - i); 
  }
  va_end(ap);
} 
  
