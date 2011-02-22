typedef struct cwin {
  int size, start, start_local;
  double *covariances;
} COV_WIN, *COV_WIN_PTR;

COV_WIN_PTR new_window(const int size, const int start);
void free_window(COV_WIN_PTR win);
void move_window(COV_WIN_PTR win, const int new_start);
void get_row(COV_WIN_PTR win, const int i, double *row, 
	     double(*evaluate)(int, int, va_list), ...);
void get_diag(COV_WIN_PTR win, double *diag, 
	      double(*evaluate)(int, int, va_list), ...);
