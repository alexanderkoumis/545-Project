/*============================================================================
==============================================================================
                      
                              balance_task.c
 
==============================================================================
Remarks:

      sekeleton to create the sample task

============================================================================*/

// system headers
#include "SL_system_headers.h"

#include <fstream>
/* SL includes */
#include "SL.h"
#include "SL_user.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_kinematics.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_man.h"

// defines

// local variables
static double      start_time = 0.0;
static SL_DJstate  target[N_DOFS+1];
static SL_Cstate   cog_target;
static SL_Cstate   cog_traj;
static SL_Cstate   cog_ref;
static double      delta_t = 0.01;
static double      duration = 4.0;
static double      time_to_go;
static int         which_step;

// possible states of a state machine
enum Steps {
  COG_RIGHT,
  LEFT_LEG_UP,
  LEFT_LEG_DOWN,
  COG_LEFT,
  RIGHT_LEG_UP,
  RIGHT_LEG_DOWN,
  DONE,
};
static int num_steps = DONE;

#define IS_COG_STEP(STEP)  \
  (STEP == COG_RIGHT) || \
  (STEP == COG_LEFT)

// variables for COG control
static iMatrix     stat;
static Matrix      Jccogp;
static Matrix      NJccog;
static Matrix      fc;
static Vector      thd;
static Vector      cog_ref_xd_mat;

// global functions 
extern "C" void
add_balance_task( void );
// 
// local functions
static int  init_balance_task(void);
static int  run_balance_task(void);
static int  change_balance_task(void);

static int  min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
        double t_togo, double dt,
        double *x_next, double *xd_next, double *xdd_next);

/*** Added functions for the project ***/
static void min_jerk_joints();
static void min_jerk_cog();
static void move_cog();
static void next_step();
static void assign_next_target();


void add_balance_task( void )
{
  int i, j;
  
  addTask("Balance Task", init_balance_task, 
    run_balance_task, change_balance_task);

}

static int init_balance_task(void)
{
  int j, i;
  int ans;
  static int firsttime = TRUE;
  
  if (firsttime) {
    firsttime = FALSE;

    // allocate memory
    stat   = my_imatrix(1,N_ENDEFFS,1,2*N_CART);
    Jccogp = my_matrix(1,N_DOFS,1,N_CART);
    NJccog = my_matrix(1,N_DOFS,1,N_DOFS+2*N_CART);
    fc     = my_matrix(1,N_ENDEFFS,1,2*N_CART);

    thd = my_vector(1, N_DOFS);
    cog_ref_xd_mat = my_vector(1, N_CART);

  }

  // this is an indicator which Cartesian components of the endeffectors are constraints
  // i.e., both feet are on the ground and cannot move in position or orientation
  for (i=1; i<=2*N_CART; ++i) {
    stat[RIGHT_FOOT][i] = TRUE;
    stat[LEFT_FOOT][i] = TRUE;
  }

  // prepare going to the default posture
  bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
  for (i=1; i<=N_DOFS; i++) {
    target[i] = joint_default_state[i];
  }

  // Initialize cog_target and cog_traj to 0s
  bzero((void *)&cog_target,sizeof(cog_target));
  bzero((void *)&cog_traj,sizeof(cog_traj));

  // go to the target using inverse dynamics (ID)
  if (!go_target_wait_ID(target))
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
    start_time, task_servo_time);

  // start data collection
  scd();

  // next_step() will take us to step 0 and initialize properly
  which_step = -1;
  next_step();

  return TRUE;
}

/*
  Computes desired joint state from target joint state using min jerk movement.
*/
static void min_jerk_joints()
{
  int i;
  for (i=1; i<=N_DOFS; ++i) {
    min_jerk_next_step(joint_des_state[i].th,
     joint_des_state[i].thd,
     joint_des_state[i].thdd,
     target[i].th,
     target[i].thd,
     target[i].thdd,
     time_to_go,
     delta_t,
     &(joint_des_state[i].th),
     &(joint_des_state[i].thd),
     &(joint_des_state[i].thdd));
  }
}

/*
  Computes desired trajectory of the center of gravity using min jerk movement.
*/
static void min_jerk_cog()
{
  int i;
  for (i=1; i<=N_CART; ++i) {
    min_jerk_next_step(
      cog_traj.x[i],
      cog_traj.xd[i],
      cog_traj.xdd[i],
      cog_target.x[i],
      cog_target.xd[i],
      cog_target.xdd[i],
      time_to_go,
      delta_t,
      &(cog_traj.x[i]),
      &(cog_traj.xd[i]),
      &(cog_traj.xdd[i]));
  }
}

/*
  Moves the center of gravity to the cartesian coordinates (x, y, z).
*/
static void move_cog()
{
  double kp = 0.1;

  compute_cog_kinematics(stat, TRUE, FALSE, TRUE, Jccogp, NJccog);

  min_jerk_cog();

  for (int i=1; i<=N_CART; ++i){
    cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i];
    cog_ref_xd_mat[i] = cog_ref.xd[i];    
  }
  
  mat_vec_mult(Jccogp, cog_ref_xd_mat, thd);

  // compute the joint_des_state[i].th and joint_des_state[i].thd  
  for (int i=1; i<=N_DOFS; ++i) {
    joint_des_state[i].th = thd[i] * delta_t + joint_des_state[i].th;
    joint_des_state[i].thd = thd[i];
    joint_des_state[i].thdd = 0;
    joint_des_state[i].uff  = 0;
  }

  inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
}

/*
  Go to the next step and prepare for it.
*/
static void next_step()
{
  // Go to the next step, repeat cycle if done
  which_step = (which_step + 1) % num_steps;
  time_to_go = duration;

  // TODO: insert freeze condition.

  // Assign the target for next step
  assign_next_target();
}

static void assign_next_target()
{
  int i;

  // Initialize variables holding the target state
  for (int i=1; i<=2*N_CART; ++i) {
      cog_traj.x[i] = cog_des.x[i];
  }
  for (int i=1; i<=N_DOFS; ++i) {
      target[i] = joint_des_state[i];
  }

  /*
  Modifiy next target state manually for each step to simulate a
  walking movement (or whatever you want).

  The following values will perform a single step using the left leg.
  */
  switch (which_step) {

    case COG_RIGHT:
      /* Move center of gravity to the right */
      cog_target.x[_X_] = cart_state[RIGHT_FOOT].x[_X_];
      break;

    case LEFT_LEG_UP:
      /* Lift left leg up */
      target[L_HFE].th += 0.8;
      target[L_KFE].th += 0.8;
      target[L_AFE].th += 0.2;
      break;

    case LEFT_LEG_DOWN:
      /* Bend right leg to lower body */
      target[R_HFE].th += 0.1;
      target[R_KFE].th += 0.6;
      target[R_AFE].th += 0.5;

      target[L_AAA].th -= 0.15;
      break;

    case COG_LEFT:
      cog_target.x[_X_] = cart_state[LEFT_FOOT].x[_X_] + 0.01;  // less left
      cog_target.x[_Y_] = cart_state[LEFT_FOOT].x[_Y_] + 0.01;  // more forward
      break;

    case RIGHT_LEG_UP:
      /* Straighten left leg */
      target[L_HFE].th -= 0.5;
      target[L_KFE].th -= 0.2;
      target[L_AFE].th -= 0.1;

      /* Carry right leg (bend a little) */
      target[R_HFE].th += 0.5;
      target[R_KFE].th += 0.5;
      target[R_AFE].th += 0.2;
      /* I don't know lol */
      target[R_FB].th  += 0.1;
      break;

    case RIGHT_LEG_DOWN:
      /* Straighten out (go to default state) */
      for (i=1; i<=N_DOFS; i++) {
        target[i] = joint_default_state[i];
      }
      break;

    case DONE:
      // Do nothing
      time_to_go = 0;
      break;

    default:
      // Do nothing
      time_to_go = 0;
      break;
  }
}

static int run_balance_task(void)
{

  // XXX: remove this when running on real nao!!
  printf("%d, %f\n", which_step, time_to_go);

  /* Initialize target/trajectory for next step */
  if (IS_COG_STEP(which_step)) {
    move_cog();
  } else {
    min_jerk_joints();
  }

  // time_step
  time_to_go -= delta_t;
  if (time_to_go <= 0) {
    next_step();
  }

  return TRUE;
}

static int change_balance_task(void)
{
  int    ivar = 0;
  double dvar = 0;

  get_int("This is how to enter an integer variable",ivar,&ivar);
  get_double("This is how to enter a double variable",dvar,&dvar);

  return TRUE;

}

static int min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
        double t_togo, double dt,
        double *x_next, double *xd_next, double *xdd_next)
{
  double t1,t2,t3,t4,t5;
  double tau,tau1,tau2,tau3,tau4,tau5;
  int    i,j;

  // a safety check
  if (dt > t_togo || dt <= 0) {
    return FALSE;
  }

  t1 = dt;
  t2 = t1 * dt;
  t3 = t2 * dt;
  t4 = t3 * dt;
  t5 = t4 * dt;

  tau = tau1 = t_togo;
  tau2 = tau1 * tau;
  tau3 = tau2 * tau;
  tau4 = tau3 * tau;
  tau5 = tau4 * tau;

  // calculate the constants
  const double dist   = t - x;
  const double p1     = t;
  const double p0     = x;
  const double a1t2   = tdd;
  const double a0t2   = xdd;
  const double v1t1   = td;
  const double v0t1   = xd;
  
  const double c1 = 6.*dist/tau5 + (a1t2 - a0t2)/(2.*tau3) - 
    3.*(v0t1 + v1t1)/tau4;
  const double c2 = -15.*dist/tau4 + (3.*a0t2 - 2.*a1t2)/(2.*tau2) +
    (8.*v0t1 + 7.*v1t1)/tau3; 
  const double c3 = 10.*dist/tau3+ (a1t2 - 3.*a0t2)/(2.*tau) -
    (6.*v0t1 + 4.*v1t1)/tau2; 
  const double c4 = xdd/2.;
  const double c5 = xd;
  const double c6 = x;
  
  *x_next   = c1*t5 + c2*t4 + c3*t3 + c4*t2 + c5*t1 + c6;
  *xd_next  = 5.*c1*t4 + 4*c2*t3 + 3*c3*t2 + 2*c4*t1 + c5;
  *xdd_next = 20.*c1*t3 + 12.*c2*t2 + 6.*c3*t1 + 2.*c4;
  
  return TRUE;
}
