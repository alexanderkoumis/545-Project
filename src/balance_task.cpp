/*============================================================================
==============================================================================
                      
                              balance_task.c
 
==============================================================================
Remarks:

      sekeleton to create the sample task

============================================================================*/

// system headers
#include "SL_system_headers.h"

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
static double      duration = 10.0;
static double      time_to_go;
static int         which_step;

// possible states of a state machine
enum Steps {
  ASSIGN_COG_TARGET,
  MOVE_TO_COG_TARGET,
  ASSIGN_JOINT_TARGET_LIFT_UP,
  MOVE_JOINT_TARGET_LIFT_UP
};

// variables for COG control
static iMatrix     stat;
static Matrix      Jccogp;
static Matrix      NJccog;
static Matrix      fc;

// global functions 
extern "C" void
add_balance_task( void );

// local functions
static int  init_balance_task(void);
static int  run_balance_task(void);
static int  change_balance_task(void);

static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next);


/*****************************************************************************
******************************************************************************
Function Name	: add_balance_task
Date		: Feb 1999
Remarks:

adds the task to the task menu

******************************************************************************
Paramters:  (i/o = input/output)

none

*****************************************************************************/
void
add_balance_task( void )
{
  int i, j;
  
  addTask("Balance Task", init_balance_task, 
	  run_balance_task, change_balance_task);

}    

/*****************************************************************************
******************************************************************************
  Function Name	: init_balance_task
  Date		: Dec. 1997

  Remarks:

  initialization for task

******************************************************************************
  Paramters:  (i/o = input/output)

       none

 *****************************************************************************/
static int 
init_balance_task(void)
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

  // ready to go
  ans = 999;
  while (ans == 999) {
    if (!get_int("Enter 1 to start or anthing else to abort ...",ans,&ans))
      return FALSE;
  }
  
  // only go when user really types the right thing
  if (ans != 1) 
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
	 start_time, task_servo_time);

  // start data collection
  scd();

  // state machine starts at ASSIGN_COG_TARGET
  which_step = ASSIGN_COG_TARGET;

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: run_balance_task
  Date		: Dec. 1997

  Remarks:

  run the task from the task servo: REAL TIME requirements!

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
run_balance_task(void)
{
  int j, i, n;

  double task_time;
  double kp = 0.1;

  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  task_time = task_servo_time - start_time;

  // the following code computes the contraint COG Jacobian 
  // Jccogp is an N_DOFS x N_CART matrix
  // NJccog is an N_DOFS x N_DOF+2*N_CART matrix -- most likely this is not needed

  compute_cog_kinematics(stat, TRUE, FALSE, TRUE, Jccogp, NJccog);

  // switch according to the current state of the state machine
  switch (which_step) {
    
  case ASSIGN_COG_TARGET:

    // what is the target for the COG?
    bzero((void *)&cog_target,sizeof(cog_target));
    cog_target.x[_X_] = ??;
    cog_target.x[_Y_] = ??;
    cog_target.x[_Z_] = ??;

    // the structure cog_des has the current position of the COG computed from the
    // joint_des_state of the robot. cog_des should track cog_traj
    bzero((void *)&cog_traj,sizeof(cog_traj));
    for (i=1; i<=N_CART; ++i)
      cog_traj.x[i] = cog_des.x[i];

    // time to go
    time_to_go = duration;

    // switch to next step of state machine
    which_step = MOVE_TO_COG_TARGET;

    break;

  case MOVE_TO_COG_TARGET: // this is for inverse kinematics control

    // plan the next step of cog with min jerk
    for (i=1; i<=N_CART; ++i) {
      min_jerk_next_step(cog_traj.x[i],
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

    // inverse kinematics: we use a P controller to correct for tracking erros
    for (i=1; i<=N_CART; ++i)
      cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i];


    // compute the joint_des_state[i].th and joint_des_state[i].thd  
    for (i=1; i<=N_DOFS; ++i) {

      // intialize to zero
      joint_des_state[i].thd  = 0;
      joint_des_state[i].thdd = 0;
      joint_des_state[i].uff  = 0;

      joint_des_state[i].thd = ?;
      joint_des_state[i].thd = ?;

    }

    // decrement time to go
    time_to_go -= delta_t;
    if (time_to_go <= 0) {
      which_step = ASSIGN_JOINT_TARGET_LIFT_UP;
    }

    break;

  case ASSIGN_JOINT_TARGET_LIFT_UP:

    // initialize the target structure from the joint_des_state
    for (i=1; i<=N_DOFS; ++i)
      target[i] = joint_des_state[i];

    time_to_go = duration;  // this may be too fast -- maybe a slower movement is better

    which_step = MOVE_JOINT_TARGET_LIFT_UP;

    break;

  case MOVE_JOINT_TARGET_LIFT_UP:

    // compute the update for the desired states
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

    // decrement time to go
    time_to_go -= delta_t;

    if (time_to_go <= 0)
      freeze();

    }

    break;

  }

  // this is a special inverse dynamics computation for a free standing robot
  inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);


  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: change_balance_task
  Date		: Dec. 1997

  Remarks:

  changes the task parameters

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
change_balance_task(void)
{
  int    ivar;
  double dvar;

  get_int("This is how to enter an integer variable",ivar,&ivar);
  get_double("This is how to enter a double variable",dvar,&dvar);

  return TRUE;

}


/*!*****************************************************************************
 *******************************************************************************
\note  min_jerk_next_step
\date  April 2014
   
\remarks 

Given the time to go, the current state is updated to the next state
using min jerk splines

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]          x,xd,xdd : the current state, vel, acceleration
 \param[in]          t,td,tdd : the target state, vel, acceleration
 \param[in]          t_togo   : time to go until target is reached
 \param[in]          dt       : time increment
 \param[in]          x_next,xd_next,xdd_next : the next state after dt

 ******************************************************************************/
static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next)

{

  // do something here

  ??
  
  return TRUE;
}

