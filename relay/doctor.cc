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


void doctor::init_test(int Frame_end_, int test_id_)
{
  Frame_end = Frame_end;

  char * test_file_name = new char[strlen(rep->out_buf) + 100];
  sprintf(test_file_name, "%s%s.re%d_test%d.redat", rep->out_buf, rep->file_name, rep->relay_id, test_id_);
  FILE * test_file = fopen(test_file_name, "r");
  printf("reading matlab test file: %s\n", test_file_name);

  int header[2];
  fread_safe(header, sizeof(int), 2, test_file);
  double *raw_params = new double[header[0]*header[1]];
  fread_safe(param_chunk[nlead], sizeof(double), header[0]*header[1], test_file);
  fclose(test_file);

  delete [] test_file_name;
}

void doctor::run_test(int Frame_end_)
{

#pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    rt->clear_event_data();
    rt->init_test_run(Frame_end_);
#pragma omp for
    for (int i = 0; i < npool; i++)
    {
      rt->test_run(pool[i], 0, Frame_end_);
      rt->write_test_run(Frame_end_);
    }
    rt->conclude_test_run(Frame_end_);
  }
}

void runner::init_test_run(int Frame_end_)
{
  TEST_pos = new double[2*n*Frame_end_];
  TEST_ref_pos = new double[2*n*Frame_end_];
  TEST_pos_res = new double[n*Frame_end_];
  TEST_alpha_INTpos_res = new double[n*Frame_end_];
  TEST_INTpos_res = new double[n*Frame_end_];
}

int runner::start_test_detection(int start_)
{
  int frame_local = start_+1;
  double dur, *f = xs+(2*n*frame_local), res_acc_local=0.0;
  advance(dur=(ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);

  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_now = (q[i].x-cx)*cl_im + cx_im, y_now = (q[i].y-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt, rsq_old=TEST_pos_res[i+(frame_local-1)*n];

    // debugging code:
    TEST_pos[j+frame_local*n]=x_now; TEST_pos[j+1+frame_local*n]=y_now;
    TEST_ref_pos[j+frame_local*n]=x_ref; TEST_ref_pos[j+1+frame_local*n]=y_ref;
    TEST_pos_res[i+frame_local*n]=rsq; TEST_alpha_INTpos_res[i+frame_local*n]=0.0; TEST_INTpos_res[i+frame_local*n]=0.5*(rsq_old+rsq)*dur;
    // end debugging code

    pos_res[start_][i]=INTpos_res[start_][i]=0.0; // zero out the first frame
    res_acc_local += pos_res[frame_local][i]=rsq;
    INTpos_res[i][2] = 0.5*(pos_res[start_][i]+rsq)*dur;
    // for now, still assuming that we start on the dot
  }
  frame_local++;
  advance(dur=(ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  f = xs+(2*n*frame_local);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_now = (q[i].x-cx)*cl_im + cx_im, y_now = (q[i].y-cy)*cl_im + cy_im,
            x_ref=f[j], y_ref=f[j+1],
            xt=x_now-x_data,yt=y_now-y_ref, rsq=xt*xt+yt*yt, rsq_old=TEST_pos_res[i+(frame_local-1)*n];

    // debugging code:
    TEST_pos[j+frame_local*n]=x_now; TEST_pos[j+1+frame_local*n]=y_now;
    TEST_ref_pos[j+frame_local*n]=x_ref; TEST_ref_pos[j+1+frame_local*n]=y_ref;
    TEST_pos_res[i+frame_local*n]=rsq; TEST_alpha_INTpos_res[i+frame_local*n]=0.0; TEST_INTpos_res[i+frame_local*n]=0.5*rsq*dur;
    // end debugging code

    res_acc_local += pos_res[frame_local][i] = rsq;
    INTpos_res[i][1] = INTpos_res[i][2]+ 0.5*(rsq_old+rsq)*dur;
  }
  pos_res_acc+=res_acc_local;
  return frame_local+1;
}

void runner::test_run(record * rec_, int start_, int end_)
{
  // debugging code: trying to ensure that we are computing everything we think we are computing
  bool not_printed = true;
  // end debugging code

  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  for (int i = 0; i < n; i++) kill_frames[i] = 0;

  // debugging code:
  int frame_local=start_;
  double *f_local = xs+(2*n*frame_local);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  x_now = (q[i].x-cx)*cl_im + cx_im, y_now = (q[i].y-cy)*cl_im + cy_im,
            x_ref=f_local[j], y_ref=f_local[j+1],
            xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt;
    TEST_pos[j+frame_local*n]=x_now; TEST_pos[j+1+frame_local*n]=y_now;
    TEST_ref_pos[j+frame_local*n]=x_ref; TEST_ref_pos[j+1+frame_local*n]=y_ref;
    TEST_pos_res[i+frame_local*n]=rsq; TEST_alpha_INTpos_res[i+frame_local*n]=0.0; TEST_INTpos_res[i+frame_local*n]=0.0;
  }
  // end debugging code

  int frame_back = start_test_detection(start_);

  double t_i=ts[frame-1], t_m1=ts[frame-2], res_acc_local=pos_res_acc;
  for ( frame = frame_back; frame < end_; frame++)
  {
    double t_p1 = ts[frame], dur;
    advance(dur=(t_p1-t_i)/t_phys, d_ang[frame-1], comega_s[frame], dt_sim);
    double * f = xs+(2*n*frame);

    bool all_dead=true;
    for (int i=0, j=0; i < n; i++, j+=2)
    {
      // debugging code:
      double  x_now = (q[i].x-cx)*cl_im + cx_im, y_now = (q[i].y-cy)*cl_im + cy_im,
              x_ref=f_local[j], y_ref=f_local[j+1],
              xt=x_now-x_ref,yt=y_now-y_ref, rsq=xt*xt+yt*yt, rsq_old=TEST_pos_res[i+(frame_local-1)*n];
      TEST_pos[j+frame_local*n]=x_now; TEST_pos[j+1+frame_local*n]=y_now;
      TEST_ref_pos[j+frame_local*n]=x_ref; TEST_ref_pos[j+1+frame_local*n]=y_ref;
      TEST_pos_res[i+frame_local*n]=rsq; TEST_alpha_INTpos_res[i+frame_local*n]=rsq; TEST_INTpos_res[i+frame_local*n]=TEST_INTpos_res[i+frame_local*n] + (+rsq);
      // end debugging code

      if (!(kill_frames[i])) // if this bead does not yet have a kill frame, continue search
      {
        all_dead=false; // then we still have at least one
        pos_res[frame][i] = rsq;
        INTpos_res[i][0] = INTpos_res[i][1] + 0.5*(rsq)*dur;
        double alpha_it=alpha_INTpos_res[frame-1][i]=alpha_comp(INTpos_res[i], t_m1, t_p1);
        if (alpha_it > alpha_tol) // if we just had a collision
        {
          kill_frames[i] = frame-1;
          alpha_kill[i] = alpha_it;
        }
        else res_acc_local+=rsq;
      }
    }
    if (all_dead&&not_printed)
    {
      printf("(thread %d): particle %d all dead at frame %d. \n", thread_id, rec_->global_index, frame);
      not_printed = false;
    }
    else if (++frame==end_)
    {
      for (int i=0; i < n; i++) if (!(kill_frames[i]))
      {
        kill_frames[i] = end_-1;
        alpha_kill[i] = alpha_INTpos_res[end_-2][i];
      }
    }
    t_m1 = t_i; t_i = t_p1;
  }
  for (int i = 0; i < n; i++) event_frame_count[i][kill_frames[i]]++;

}
