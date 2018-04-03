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
static double      duration = 5.0;
// static double      duration = 10.0;
static double      time_to_go;
static int         which_step;

// possible states of a state machine
enum Steps {
  ASSIGN_COG_TARGET_RIGHT     = 0,
  MOVE_COG_TARGET_RIGHT       = 1,
  ASSIGN_JOINT_LEFT_LEG_UP    = 2,
  MOVE_JOINT_LEFT_LEG_UP      = 3,
  ASSIGN_JOINT_LEFT_LEG_DOWN  = 4,
  MOVE_JOINT_LEFT_LEG_DOWN    = 5,
  ASSIGN_COG_TARGET_LEFT      = 6,
  MOVE_COG_TARGET_LEFT        = 7,
  ASSIGN_JOINT_RIGHT_LEG_UP   = 8,
  MOVE_JOINT_RIGHT_LEG_UP     = 9,
  ASSIGN_JOINT_RIGHT_LEG_DOWN = 10,
  MOVE_JOINT_RIGHT_LEG_DOWN   = 11
};

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
  
  if (firsttime){
    firsttime = FALSE;

    // allocate memory
    stat   = my_imatrix(1,N_ENDEFFS,1,2*N_CART);
    Jccogp = my_matrix(1,N_DOFS,1,N_CART);
    NJccog = my_matrix(1,N_DOFS,1,N_DOFS+2*N_CART);
    fc     = my_matrix(1,N_ENDEFFS,1,2*N_CART);

    thd = my_vector(1, N_DOFS);
    cog_ref_xd_mat = my_vector(1, N_CART);

    // this is an indicator which Cartesian components of the endeffectors are constraints
    // i.e., both feet are on the ground and cannot move in position or orientation
    stat[RIGHT_FOOT][1] = TRUE;
    stat[RIGHT_FOOT][2] = TRUE;
    stat[RIGHT_FOOT][3] = TRUE;
    stat[RIGHT_FOOT][4] = TRUE;
    stat[RIGHT_FOOT][5] = TRUE;
    stat[RIGHT_FOOT][6] = TRUE;

    stat[LEFT_FOOT][1] = TRUE;
    stat[LEFT_FOOT][2] = TRUE;
    stat[LEFT_FOOT][3] = TRUE;
    stat[LEFT_FOOT][4] = TRUE;
    stat[LEFT_FOOT][5] = TRUE;
    stat[LEFT_FOOT][6] = TRUE;

  }

  // prepare going to the default posture
  bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
  for (i=1; i<=N_DOFS; i++)
    target[i] = joint_default_state[i];

  // go to the target using inverse dynamics (ID)
  if (!go_target_wait_ID(target)) 
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
   start_time, task_servo_time);

  // start data collection
  scd();

  // state machine starts at ASSIGN_COG_TARGET
  which_step = ASSIGN_COG_TARGET_RIGHT;

  return TRUE;
}


static void min_jerk_cog()
{
  for (int i=1; i<=N_CART; ++i) {
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

static void min_jerk_joints()
{
  for (int i = 1; i <= N_DOFS; ++i) {
    SL_DJstate &joint_ = joint_des_state[i];
    SL_DJstate &target_ = target[i];
    min_jerk_next_step(
        joint_.th,
        joint_.thd,
        joint_.thdd,
        target_.th,
        target_.thd,
        target_.thdd,
        time_to_go,
        delta_t,
        &joint_.th,
        &joint_.thd,
        &joint_.thdd
    );
  }
}


static Steps time_step(Steps current_step, Steps next_step, bool should_freeze=false)
{
  time_to_go -= delta_t;
  if (time_to_go <= 0) {
    if (should_freeze) {
      freeze();
    } else {
      return next_step;
    }
  }
  return current_step;
}


static void cog_assignment(float x, float y, float z)
{
  bzero((void *)&cog_target,sizeof(cog_target));

  cog_target.x[_X_] = x;
  cog_target.x[_Y_] = y;
  cog_target.x[_Z_] = z;

  for (int i=1; i<=6; ++i) {
    stat[RIGHT_FOOT][i] = TRUE;
    stat[LEFT_FOOT][i] = TRUE;
  }

  // the structure cog_des has the current position of the COG computed from the
  // joint_des_state of the robot. cog_des should track cog_traj
  bzero((void *)&cog_traj,sizeof(cog_traj));
  for (int i=1; i<=N_CART; ++i) {
    cog_traj.x[i] = cog_des.x[i];
  }
}


static void cog_computation()
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

    // intialize to zero
    joint_des_state[i].thd  = 0;
    joint_des_state[i].thdd = 0;
    joint_des_state[i].uff  = 0;

    joint_des_state[i].th = thd[i] * delta_t + joint_des_state[i].th;
    joint_des_state[i].thd = thd[i];
    joint_des_state[i].thdd = 0;
    joint_des_state[i].uff  = 0;
  }
}


static Steps state_assign_cog_target_right()
{
  cog_assignment(
    cart_des_state[RIGHT_FOOT].x[_X_] * 1.2,
    cart_des_state[RIGHT_FOOT].x[_Y_],
    cart_des_state[RIGHT_FOOT].x[_Z_]
  );

  // time to go
  time_to_go = 10.0;
  // time_to_go = 5.0;

  // switch to next step of state machine
  return MOVE_COG_TARGET_RIGHT;
}


static Steps state_move_cog_target_right()
{
  cog_computation();
  // this is a special inverse dynamics computation for a free standing robot
  inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
  return time_step(MOVE_COG_TARGET_RIGHT, ASSIGN_JOINT_LEFT_LEG_UP);
}


static Steps state_assign_joint_left_leg_up()
{
  // initialize the target structure from the joint_des_state
  for (int i=1; i<=N_DOFS; ++i) {
    target[i] = joint_des_state[i];
  }

  target[L_HFE].th =  0.6;
  target[L_HAA].th = -0.4;

  // time_to_go = 10.0;
  time_to_go = 5.0;

  return MOVE_JOINT_LEFT_LEG_UP;
}

static Steps state_move_joint_left_leg_up()
{

  min_jerk_joints();
  return time_step(MOVE_JOINT_LEFT_LEG_UP, ASSIGN_JOINT_LEFT_LEG_DOWN);
}


static Steps state_assign_joint_left_leg_down()
{
  for (int i=1; i<=N_DOFS; ++i) {
    target[i] = joint_des_state[i];
  }

  target[L_HAA].th = -0.25;
  target[R_HFE].th = 0.2;
  target[R_AFE].th = 0.05;

  time_to_go = 5.0;

  return MOVE_JOINT_LEFT_LEG_DOWN;
}

static Steps state_move_joint_left_leg_down()
{
  min_jerk_joints();
  return time_step(MOVE_JOINT_LEFT_LEG_DOWN, ASSIGN_COG_TARGET_LEFT);
}

static Steps state_assign_cog_target_left()
{
  cog_assignment(
      cart_des_state[LEFT_FOOT].x[_X_],
      cart_des_state[LEFT_FOOT].x[_Y_],
      cart_des_state[LEFT_FOOT].x[_Z_]
  );
  time_to_go = 10.0;
  return MOVE_COG_TARGET_LEFT;
}


static Steps state_move_cog_target_left()
{
  cog_computation();
  // this is a special inverse dynamics computation for a free standing robot
  inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
  return time_step(MOVE_COG_TARGET_LEFT, ASSIGN_JOINT_RIGHT_LEG_UP);
}


static Steps state_assign_joint_right_leg_up()
{
  // initialize the target structure from the joint_des_state
  for (int i=1; i<=N_DOFS; ++i) {
    target[i] = joint_des_state[i];
  }

  target[L_HFE].th =  0.6;
  target[R_HFE].th =  0.6;
  target[R_HAA].th = -0.4;

  time_to_go = 10.0;

  return MOVE_JOINT_RIGHT_LEG_UP;
}

static Steps state_move_joint_right_leg_up()
{

  min_jerk_joints();
  return time_step(MOVE_JOINT_RIGHT_LEG_UP, ASSIGN_JOINT_RIGHT_LEG_DOWN);
}


static Steps state_assign_joint_right_leg_down()
{
  for (int i=1; i<=N_DOFS; ++i) {
    target[i] = joint_des_state[i];
  }

  target[R_HAA].th = -0.25;

  time_to_go = 10.0;

  return MOVE_JOINT_RIGHT_LEG_DOWN;
}

static Steps state_move_joint_right_leg_down()
{
  min_jerk_joints();
  return time_step(MOVE_JOINT_RIGHT_LEG_DOWN, ASSIGN_COG_TARGET_RIGHT);
}




static int run_balance_task(void)
{

  printf("%d, %f\n", which_step, time_to_go);
  switch (which_step) {
    case ASSIGN_COG_TARGET_RIGHT     : which_step = state_assign_cog_target_right();    break;
    case MOVE_COG_TARGET_RIGHT       : which_step = state_move_cog_target_right();      break;
    case ASSIGN_JOINT_LEFT_LEG_UP    : which_step = state_assign_joint_left_leg_up();   break;
    case MOVE_JOINT_LEFT_LEG_UP      : which_step = state_move_joint_left_leg_up();     break;
    case ASSIGN_JOINT_LEFT_LEG_DOWN  : which_step = state_assign_joint_left_leg_down(); break;
    case MOVE_JOINT_LEFT_LEG_DOWN    : which_step = state_move_joint_left_leg_down();   break;
    case ASSIGN_COG_TARGET_LEFT      : which_step = state_assign_cog_target_left();     break;
    case MOVE_COG_TARGET_LEFT        : which_step = state_move_cog_target_left();       break;
    case ASSIGN_JOINT_RIGHT_LEG_UP   : which_step = state_assign_joint_right_leg_up();   break;
    case MOVE_JOINT_RIGHT_LEG_UP     : which_step = state_move_joint_right_leg_up();     break;
    case ASSIGN_JOINT_RIGHT_LEG_DOWN : which_step = state_assign_joint_right_leg_down(); break;
    case MOVE_JOINT_RIGHT_LEG_DOWN   : which_step = state_move_joint_right_leg_down();   break;
  }

  return TRUE;
}

static int change_balance_task(void)
{
  int    ivar;
  double dvar;

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
