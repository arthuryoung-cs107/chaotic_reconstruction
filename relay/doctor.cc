#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "particle_relay.hh"

extern "C"
{
  #include "AYaux.h"
}


void doctor::init_test(int test_id_, int test_relay_id_)
{
  test_buffer = new char[strlen(rep->out_buf) + 100];
  sprintf(test_buffer, "%s%s.re%d_test%d.redat", rep->out_buf, rep->file_name, test_relay_id_, test_id_);
  FILE * test_file = fopen(test_buffer, "r");
  printf("reading matlab test file: %s\n", test_buffer);

  int header[2];
  fread_safe(header, sizeof(int), 2, test_file);
  double *raw_params = new double[header[0]*header[1]];
  fread_safe(param_chunk[nlead], sizeof(double), header[0]*header[1], test_file);
  fclose(test_file);

  sprintf(test_buffer, "%s%s.re%d_test%d_results/", rep->out_buf, rep->file_name, test_relay_id_, test_id_); mkdir(test_buffer, S_IRWXU);
  printf("Made test directory: %s\n", test_buffer);
}

void doctor::test_run(int Frame_end_)
{
  int int_param_len = 3;
  int int_params[] = {Frame_end_, n, 2};
  char * test_run_spec_name = new char[strlen(test_buffer)+50];
  sprintf(test_run_spec_name,"%sresults_specs.redat", test_buffer);
  FILE * test_specs = fopen(test_run_spec_name, "wb");
  fwrite(int_params, sizeof(int), int_param_len, test_specs);
  fclose(test_specs);

#pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    rt->init_test_run(Frame_end_);
    medic med(Frame_end_, rt, test_buffer);
    rt->clear_event_data();
#pragma omp for
    for (int i = 0; i < npool; i++)
    {
      med.write_test_run(rt->run_test(pool[i], 0, Frame_end_, i),i);
    }
    rt->conclude_test_run();
  }
}

void medic::write_test_run(double accres_, int i_)
{
  int header[] = {beads, Frame_end};
  double accres_local=accres_;

  sprintf(test_directory+buf_end, "par%d.redat", i_);
  FILE * test_file = fopen(test_directory, "wb");

  fwrite(header, sizeof(int), 2, test_file);
  fwrite(&accres_local, sizeof(double), 1, test_file);
  fwrite(rt->TEST_sim_pos, sizeof(double), 2*beads*Frame_end, test_file);
  fwrite(rt->TEST_pos, sizeof(double), 2*beads*Frame_end, test_file);
  fwrite(rt->TEST_ref_pos, sizeof(double), 2*beads*Frame_end, test_file);
  fwrite(rt->TEST_pos_res, sizeof(double), beads*Frame_end, test_file);
  fwrite(rt->TEST_alpha_INTpos_res, sizeof(double), beads*Frame_end, test_file);
  fwrite(rt->TEST_INTpos_res, sizeof(double), beads*Frame_end, test_file);

  fclose(test_file);
}

void runner::init_test_run(int Frame_end_)
{
  TEST_sim_pos=new double[2*n*Frame_end_]; TEST_pos=new double[2*n*Frame_end_];
  TEST_ref_pos=new double[2*n*Frame_end_]; TEST_pos_res=new double[n*Frame_end_];
  TEST_alpha_INTpos_res=new double[n*Frame_end_]; TEST_INTpos_res=new double[n*Frame_end_];
}

void runner::conclude_test_run()
{
  delete [] TEST_sim_pos; delete [] TEST_pos;
  delete [] TEST_ref_pos; delete [] TEST_pos_res;
  delete [] TEST_alpha_INTpos_res; delete [] TEST_INTpos_res;
}

double runner::run_test(record * rec_, int start_, int end_, int i_)
{
  double t_vec[3];
  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  int frame_local=start_; // t0
  double dur, res_acc_local = 0.0;
  double * f = xs+(2*n*frame_local);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref,rsq=xt*xt+yt*yt;
    TEST_sim_pos[j+(2*frame_local*n)]=x_sim; TEST_sim_pos[j+1+(2*frame_local*n)]=y_sim;
    TEST_pos[j+(2*frame_local*n)]=x_now; TEST_pos[j+1+(2*frame_local*n)]=y_now;
    TEST_ref_pos[j+(2*frame_local*n)]=x_ref; TEST_ref_pos[j+1+(2*frame_local*n)]=y_ref;

    TEST_pos_res[i+(frame_local*n)]=rsq;
    TEST_INTpos_res[i+(frame_local*n)]=0.0;
    TEST_alpha_INTpos_res[i+(frame_local*n)]=0.0; // doesn't matter
    res_acc_local+=rsq;
  }
  frame_local++; f = xs+(2*n*frame_local); //t1
  advance(dur=((ts[frame_local]-ts[frame_local-1])/t_phys), d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref,rsq=xt*xt+yt*yt, rsq_old=TEST_pos_res[i+((frame_local-1)*n)]; // technically unnecessary, but keeping for consistency

    TEST_sim_pos[j+(2*frame_local*n)]=x_sim; TEST_sim_pos[j+1+(2*frame_local*n)]=y_sim;
    TEST_pos[j+(2*frame_local*n)]=x_now; TEST_pos[j+1+(2*frame_local*n)]=y_now;
    TEST_ref_pos[j+(2*frame_local*n)]=x_ref; TEST_ref_pos[j+1+(2*frame_local*n)]=y_ref;

    TEST_pos_res[i+(frame_local*n)]=rsq;
    TEST_INTpos_res[i+(frame_local*n)]=TEST_INTpos_res[i+(frame_local*n)]+0.5*(rsq_old+rsq)*dur; // ditto above
    TEST_alpha_INTpos_res[i+(frame_local*n)]=0.0; // doesn't matter
    res_acc_local+=rsq;
  }
  frame_local++; f = xs+(2*n*frame_local); // t2
  advance(dur=((ts[frame_local]-ts[frame_local-1])/t_phys), d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref,rsq=xt*xt+yt*yt, rsq_old=TEST_pos_res[i+((frame_local-1)*n)];

    TEST_sim_pos[j+(2*frame_local*n)]=x_sim; TEST_sim_pos[j+1+(2*frame_local*n)]=y_sim;
    TEST_pos[j+(2*frame_local*n)]=x_now; TEST_pos[j+1+(2*frame_local*n)]=y_now;
    TEST_ref_pos[j+(2*frame_local*n)]=x_ref; TEST_ref_pos[j+1+(2*frame_local*n)]=y_ref;

    TEST_pos_res[i+(frame_local*n)]=rsq;
    TEST_INTpos_res[i+(frame_local*n)]=TEST_INTpos_res[i+((frame_local-1)*n)]+0.5*(rsq_old+rsq)*dur;
    TEST_alpha_INTpos_res[i+(frame_local*n)]=0.0; // rewritten next step
    res_acc_local+=rsq;
  }
  t_vec[1]=ts[frame_local]; t_vec[2]=ts[frame_local-1];
  int frame_body_start=frame_local+1;
  for (frame_local = frame_body_start; frame_local < end_; frame_local++)
  {
    advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
    t_vec[0]=ts[frame_local];
    for (int i=0, j=0; i < n; i++, j+=2)
    {
      double  x_sim = q[i].x, y_sim = q[i].y,
              x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
              x_ref=f[j], y_ref=f[j+1],
              xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt, rsq_old=TEST_pos_res[i+(frame_local-1)*n];

      TEST_sim_pos[j+(2*frame_local*n)]=x_sim; TEST_sim_pos[j+1+(2*frame_local*n)]=y_sim;
      TEST_pos[j+(2*frame_local*n)]=x_now; TEST_pos[j+1+(2*frame_local*n)]=y_now;
      TEST_ref_pos[j+(2*frame_local*n)]=x_ref; TEST_ref_pos[j+1+(2*frame_local*n)]=y_ref;

      TEST_pos_res[i+frame_local*n]=rsq;
      double Ii = TEST_INTpos_res[i+frame_local*n]=TEST_INTpos_res[i+((frame_local-1)*n)]+0.5*(rsq_old+rsq)*dur,
            Im2 = TEST_INTpos_res[i+(frame_local-2)*n];
      TEST_alpha_INTpos_res[i+(frame_local-1)*n]=log(Ii/Im2)/log(t_vec[0]/t_vec[2]);
      res_acc_local+=rsq;
    }
    t_vec[2]=t_vec[1]; t_vec[1]=t_vec[0];
  }
  printf("(thread %d) ran parameters %d\n", thread_id, i_);
  return res_acc_local;
}
