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
static SL_Cstate   cog_desired;
static SL_Cstate   cog_ref;
static double      delta_t = 0.01;
static double      duration = 10.0;
static double      time_to_go;
static int         which_step;

enum Steps {
  ASSIGN_COG_TARGET,
  MOVE_TO_COG_TARGET,
  ASSIGN_JOINT_TARGET_LIFT_UP,
  MOVE_JOINT_TARGET_LIFT_UP
};

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

    stat   = my_imatrix(1,N_ENDEFFS,1,2*N_CART);
    Jccogp = my_matrix(1,N_DOFS,1,N_CART);
    NJccog = my_matrix(1,N_DOFS,1,N_DOFS+2*N_CART);
    fc     = my_matrix(1,N_ENDEFFS,1,2*N_CART);


  }

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


  // prepare going to the default posture
  bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
  for (i=1; i<=N_DOFS; i++)
    target[i] = joint_default_state[i];

  ans = 0;
  get_int("Go to one leg balance posture?",ans,&ans);
  if (ans) {
    target[R_SAA].th = -1.5;
    target[R_HAA].th = 0.25;
    target[R_AAA].th = 0.25;
    target[L_HAA].th = -0.35;
    target[L_AAA].th = 0.25;
    target[L_HFE].th = 0.3;
    target[L_KFE].th = 0.8;
    target[L_AFE].th = 0.5;
  }

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


  //  target[R_SFE].th  = -1.83;
  //target[L_SFE].th  = -1.83;
  target[R_HAA].th  =  0.22;
  target[R_AAA].th  =  0.22;
  target[L_HAA].th  = -0.22;
  target[L_AAA].th  =  0.26;
  time_to_go = duration;

  //which_step = MOVE_JOINT_TARGET_LIFT_UP;

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
  double aux;
  static int firsttime = TRUE;
  static int wait_ticks;
  double s = 0.1;


  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  task_time = task_servo_time - start_time;

  // the following code computes the contraint COG Jacobian 
  // Jccogp is an N_DOFS x N_CART matrix
  // NJccog is an N_DOFS x N_DOF+2*N_CART matrix

  compute_cog_kinematics(stat, TRUE, FALSE, FALSE, Jccogp, NJccog);
  
  /*
  if (firsttime) {
    FILE *fp;
    fp = fopen("mist.mat","r");
    fread_mat(fp,Jccogp);
    fclose(fp);
    firsttime = FALSE;
    print_mat("Jccogp",Jccogp);
  }
  */

  //print_mat("Jccogp",Jccogp);
  //print_mat("Jcog",Jcog);

  //print_mat("Jccogp",Jccogp);
  //getchar();

  switch (which_step) {
    
  case ASSIGN_COG_TARGET:

    bzero((void *)&cog_target,sizeof(cog_target));
    cog_target.x[_X_] = 0.05;
    cog_target.x[_Y_] = 0.00;
    cog_target.x[_Z_] = -.13;

    bzero((void *)&cog_desired,sizeof(cog_desired));
    for (i=1; i<=N_CART; ++i)
      cog_desired.x[i] = cog_des.x[i];

    // time to go
    time_to_go = duration;

    // switch to next step of state machine
    which_step = MOVE_TO_COG_TARGET;

    break;

  case MOVE_TO_COG_TARGET:

    // plan the next step of cog with min jerk
    for (i=1; i<=N_CART; ++i) {
      min_jerk_next_step(cog_desired.x[i],
			 cog_desired.xd[i],
			 cog_desired.xdd[i],
			 cog_target.x[i],
			 cog_target.xd[i],
			 cog_target.xdd[i],
			 time_to_go,
			 delta_t,
			 &(cog_desired.x[i]),
			 &(cog_desired.xd[i]),
			 &(cog_desired.xdd[i]));
    }

    // inverse kinematics
    for (i=1; i<=N_CART; ++i)
      cog_ref.xd[i] = kp*(cog_desired.x[i] - cog_des.x[i]) + cog_desired.xd[i];

    //printf("cog_ref = %f %f %f\n",cog_ref.xd[1],cog_ref.xd[2],cog_ref.xd[3]);

    // initialize
    for (i=1; i<=N_DOFS; ++i) {
      joint_des_state[i].thd  = 0;
      joint_des_state[i].thdd = 0;
    }

    for (i=1; i<=N_DOFS; ++i) {
      for (j=1; j<=N_CART; ++j) {
	joint_des_state[i].thd += cog_ref.xd[j] * Jccogp[i][j];
      }
    }

    // integrate the desired state
    for (i=1; i<=N_DOFS; ++i) 
      joint_des_state[i].th += joint_des_state[i].thd * delta_t;


    // decrement time to go
    time_to_go -= delta_t;
    if (time_to_go <= 0) {
      which_step = ASSIGN_JOINT_TARGET_LIFT_UP;
    }

    break;

  case ASSIGN_JOINT_TARGET_LIFT_UP:
    for (i=1; i<=N_DOFS; ++i)
      target[i] = joint_des_state[i];

    target[L_HFE].th += s;
    target[L_KFE].th += 2*s;
    target[L_AFE].th += s;
    target[L_HAA].th -= 0.3;
    target[L_AAA].th  = 0.38;

    target[R_HAA].th -= 0.3;
    //target[R_AAA].th += 0.1;
    time_to_go = duration*4;

    for (n=1; n<=6; ++n)
      stat[LEFT_FOOT][n] = FALSE;

    wait_ticks = task_servo_rate;

    which_step = MOVE_JOINT_TARGET_LIFT_UP;

    break;

  case MOVE_JOINT_TARGET_LIFT_UP:

    if (--wait_ticks > 0)
      break;

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
    if (time_to_go <= 0) {
      freeze();
    }

    break;

  }


  if (which_step > ASSIGN_JOINT_TARGET_LIFT_UP) {

    ;


  } else {


    for (i=1; i<=N_DOFS; ++i)
      joint_des_state[i].uff = 0.0;
    
    inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);
    //double xdd_ref[3+1] = {0,0,0,0};
    //double add_ref[3+1] = {0,0,0,0};
    //inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, xdd_ref, add_ref, fc);
  }



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

