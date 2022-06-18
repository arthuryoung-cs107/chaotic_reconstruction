#include "MH_solvers.hh"

void MH_doctor::stage_diagnostics()
{
  for (int frame_i = 0; frame_i < Frames_test; frame_i++)
  {
    int foffset = 2*nbeads*frame_i;
    double *f = xs+(foffset);
    for (int bead_i = 0, j=0; bead_i < nbeads; bead_i++, j+=2)
    {TEST_refp[j+foffset]=f[j]; TEST_refp[j+1+foffset]=f[j+1];}
  }
  int MH_doc_header_len = 1;
  int header_ints[] = {MH_doc_header_len, Frames_test};
  sprintf(test_buffer+test_buf_end, "results_startspecs.redat");
  FILE * test_startspecs = fopen(test_buffer, "wb");
  write_MH_params(test_startspecs);
  fwrite(header_ints, sizeof(int), MH_doc_header_len+1, test_startspecs);
  fwrite(TEST_refp, sizeof(double), 2*nbeads*Frames_test, test_startspecs);
  fclose(test_startspecs);
}
void MH_doctor::close_diagnostics()
{
  int MH_doc_header_len = 0;
  int header_ints[] = {MH_doc_header_len};
  sprintf(test_buffer+test_buf_end, "results_endspecs.redat");
  FILE * test_endspecs = fopen(test_buffer, "wb");
  fwrite(header_ints, sizeof(int), MH_doc_header_len+1, test_endspecs);
  fwrite(evcount_bead_frame[0], sizeof(int), 2*nbeads*Frames_test, test_endspecs);
  fclose(test_endspecs);
}
