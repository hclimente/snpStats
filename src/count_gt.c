void count_gt(unsigned char *a, unsigned char *b, int *snp_count, int *sample_count, int *output, int *output_signed) {
  for(int j = 0; j < *snp_count ; j++) {
    for(int i = 0; i < *sample_count ; i++) {
      if ((*a) != (*b)) {
        (*output) ++ ;
        if ((*b) != 0x00) {
          (*output_signed)++;
        }
        if ((*a) != 0x00) {
          (*output_signed)--;
        }
      }
      a++;
      b++;
    }
    output++;
    output_signed++;
  }
}
