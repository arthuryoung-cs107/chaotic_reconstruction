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
  bool first2finish = true;

  TEST_ref_pos = new double[2*n*Frame_end_];
  for (int frame_i = 0; frame_i < Frame_end_; frame_i++)
  {
    int foffset = 2*n*frame_i;
    double *f = xs+(foffset);
    for (int bead_i = 0, j=0; bead_i < n; bead_i++, j+=2)
    {TEST_ref_pos[j+foffset]=f[j]; TEST_ref_pos[j+1+foffset]=f[j+1];}
  }

  int int_param_len = 4;
  int int_params[] = {Frame_end_, n, 2, Frames};
  char * test_run_spec_name = new char[strlen(test_buffer)+50];
  sprintf(test_run_spec_name,"%sresults_specs.redat", test_buffer);
  FILE * test_specs = fopen(test_run_spec_name, "wb");
  fwrite(int_params, sizeof(int), int_param_len, test_specs);
  fwrite(TEST_ref_pos, sizeof(double), 2*n*Frame_end_, test_specs);
  fclose(test_specs);
  delete [] TEST_ref_pos;

#pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    rt->init_test_run(Frame_end_);
    medic med(Frame_end_, rt, test_buffer);
    rt->clear_event_data();
#pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {
      med.write_test_run(rt->test_detection(pool[i], 0, Frame_end_, i),i);
    }
#pragma omp critical
    {
      if (first2finish)
      {
        for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] = rt->event_frame_count[0][i];
        first2finish=false;
      }
      else for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] += rt->event_frame_count[0][i];
    }
    rt->free_test_buffers();
  }

  sprintf(test_run_spec_name,"%sfinal_event_count.redat", test_buffer);
  FILE * test_event_count = fopen(test_run_spec_name, "wb");
  fwrite(global_event_frame_count, sizeof(int), n*Frames, test_event_count);

  fclose(test_event_count);
  delete [] test_run_spec_name;
}
int runner::start_test_detection(int start_, double * params_, double * t_history_)
{
  for (int i = 0; i < n*Frames; i++) pos_res[0][i]=alpha_INTpos_res[0][i]=0.0;

  reset_sim(params_, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  int frame_local = start_, poffset = 2*n*frame_local, foffset=n*frame_local;
  double res_acc_local = 0.0, *f = xs+poffset;

  // initialize t=tstart data
  t_history_[0]=t_history_[1]=ts[frame_local];
  for (int i=0,j=0; i < n; i++,j+=2)
  {
    kill_frames[i] = 0;
    pos_res[frame_local][i]=alpha_INTpos_res[frame_local][i]=INTpos_res[i][0]=INTpos_res[i][1]=INTpos_res[i][2]=0.0;

    // write test outputs
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im;

    TEST_sim_pos[j+poffset]=x_sim; TEST_sim_pos[j+1+poffset]=y_sim;
    TEST_pos[j+poffset]=x_now; TEST_pos[j+1+poffset]=y_now;

    TEST_pos_res[i+foffset]=pos_res[frame_local][i]; TEST_INTpos_res[i+foffset]=INTpos_res[i][0];
    TEST_alpha_INTpos_res[i+foffset]=alpha_INTpos_res[frame_local][i];
  }

  // step to frame 1, begin computing divergence from data
  frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
  advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[frame_local];
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt;

    res_acc_local+=pos_res[frame_local][i]=rsq;
    INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
    INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_INTpos_res[frame_local][i]=0.0;

    // write test outputs, being absolutely sure we are writing exactly what is being calculated
    TEST_sim_pos[j+poffset]=x_sim; TEST_sim_pos[j+1+poffset]=y_sim;
    TEST_pos[j+poffset]=x_now; TEST_pos[j+1+poffset]=y_now;

    TEST_pos_res[i+foffset]=pos_res[frame_local][i]; TEST_INTpos_res[i+foffset]=INTpos_res[i][0];
    TEST_alpha_INTpos_res[i+foffset]=alpha_INTpos_res[frame_local][i];
  }

  // step to frame 2, compute divergence from data as we will in the body of the detection
  frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
  advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[frame_local];
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt;

    res_acc_local+=pos_res[frame_local][i]=rsq;
    INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
    INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_INTpos_res[frame_local][i]=0.0;

    // write test outputs, being absolutely sure we are writing exactly what is being calculated
    TEST_sim_pos[j+poffset]=x_sim; TEST_sim_pos[j+1+poffset]=y_sim;
    TEST_pos[j+poffset]=x_now; TEST_pos[j+1+poffset]=y_now;

    TEST_pos_res[i+foffset]=pos_res[frame_local][i]; TEST_INTpos_res[i+foffset]=INTpos_res[i][0];
    TEST_alpha_INTpos_res[i+foffset]=alpha_INTpos_res[frame_local][i];
  }

  pos_res_acc=res_acc_local;
  return frame_local;
}

double runner::test_detection(record * rec_, int start_, int end_, int i_)
{
  double t_history[3];
  int frame_local = start_test_detection(start_, rec_->params, t_history);
  int poffset = 2*n*frame_local, foffset=n*frame_local;
  double res_acc_local = pos_res_acc, res_acc_local_full = pos_res_acc, *f;
  do
  {
    frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
    advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
    t_history[2]=t_history[1]; t_history[1]=t_history[0]; t_history[0]=ts[frame_local];

    for (int i=0, j=0; i < n; i++, j+=2)
    {
      double  x_sim = q[i].x, y_sim = q[i].y,
              x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im,
              x_ref=f[j], y_ref=f[j+1],
              xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt;

      TEST_sim_pos[j+poffset]=x_sim; TEST_sim_pos[j+1+poffset]=y_sim;
      TEST_pos[j+poffset]=x_now; TEST_pos[j+1+poffset]=y_now;
      res_acc_local_full+=rsq;
      if (!(kill_frames[i])) // if this bead does not yet have a kill frame, continue search
      {
        res_acc_local+=pos_res[frame_local][i]=rsq;
        INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
        INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history[0]-t_history[1]);
        double alpha_it=alpha_INTpos_res[frame_local-1][i]=alpha_comp(INTpos_res[i], t_history[0], t_history[2]);

        if (alpha_it > alpha_tol) // if we just had a collision
        {
          kill_frames[i] = frame_local-1;
          alpha_kill[i] = alpha_it;
        }
        else res_acc_local+=rsq;

        // write test outputs, being absolutely sure we are writing exactly what is being calculated
        TEST_pos_res[i+foffset]=pos_res[frame_local][i]; TEST_INTpos_res[i+foffset]=INTpos_res[i][0];
        TEST_alpha_INTpos_res[i+(foffset-n)]=alpha_INTpos_res[frame_local-1][i];
      }
      else
      {
        // write test outputs anyway.
        TEST_pos_res[i+foffset]=rsq;
        double INTp1=TEST_INTpos_res[i+foffset]=TEST_INTpos_res[i+(foffset-n)]+0.5*(rsq+TEST_pos_res[i+(foffset-n)])*(t_history[0]-t_history[1]);
        TEST_alpha_INTpos_res[i+(foffset-n)]=alpha_comp(INTp1,TEST_INTpos_res[i+(foffset-2*n)], t_history[0], t_history[2]);
      }
    }
    if ((frame_local+1)==end_)
    {
      for (int i=0; i < n; i++)
      {
        if (!(kill_frames[i]))
        {
          kill_frames[i] = end_-1;
          alpha_kill[i] = alpha_INTpos_res[frame_local-2][i];
        }
      }
      break;
    }
  } while(true);
  pos_res_acc=res_acc_local;
  earliest_event=Frames; latest_event=0;
  for (int i = 0; i < n; i++)
  {
    if (kill_frames[i]<earliest_event) earliest_event=kill_frames[i];
    if (kill_frames[i]>latest_event) latest_event=kill_frames[i];
    event_frame_count[i][kill_frames[i]]++;
  }
  return res_acc_local_full;
}


void medic::write_test_run(double accresfull_, int i_)
{
  int int_len=2, double_len=2;
  int header[] = {int_len, double_len};
  int int_params[] = {rt->earliest_event, rt->latest_event};
  double double_params[] = {accresfull_, rt->pos_res_acc};

  sprintf(test_directory+buf_end, "par%d.redat", i_);
  FILE * test_file = fopen(test_directory, "wb");

  fwrite(header, sizeof(int), 2, test_file);
  fwrite(int_params, sizeof(int), int_len, test_file);
  fwrite(double_params, sizeof(double), double_len, test_file);
  fwrite(rt->TEST_sim_pos, sizeof(double), 2*beads*Frame_end, test_file);
  fwrite(rt->TEST_pos, sizeof(double), 2*beads*Frame_end, test_file);
  fwrite(rt->TEST_pos_res, sizeof(double), beads*Frame_end, test_file);
  fwrite(rt->TEST_alpha_INTpos_res, sizeof(double), beads*Frame_end, test_file);
  fwrite(rt->TEST_INTpos_res, sizeof(double), beads*Frame_end, test_file);
  fwrite(rt->pos_res, sizeof(double), beads*rt->Frames, test_file);
  fwrite(rt->alpha_INTpos_res, sizeof(double), beads*rt->Frames, test_file);
  fwrite(rt->kill_frames, sizeof(int), beads, test_file);
  fwrite(rt->alpha_kill, sizeof(double), beads, test_file);
  fclose(test_file);
}

void runner::init_test_run(int Frame_end_)
{
  TEST_sim_pos=new double[2*n*Frame_end_]; TEST_pos=new double[2*n*Frame_end_];
  TEST_pos_res=new double[n*Frame_end_];
  TEST_alpha_INTpos_res=new double[n*Frame_end_]; TEST_INTpos_res=new double[n*Frame_end_];
}

void runner::free_test_buffers()
{
  delete [] TEST_sim_pos; delete [] TEST_pos;
  delete [] TEST_pos_res;
  delete [] TEST_alpha_INTpos_res; delete [] TEST_INTpos_res;
}
