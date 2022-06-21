

#include<stdbool.h> 
#include <math.h>

#include "mujoco.h"
#include "glfw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

double simend = 5;

FILE *fid;
int loop_index = 0;
const int data_frequency = 10;


char path[] = "../Projects/Double_pendulum/";
char xmlfile[] = "pendulum.xml";


char datafile[] = "data.csv";


mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
mjvCamera cam;                      // abstract camera
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     


bool button_left = false;
bool button_middle = false;
bool button_right =  false;
double lastx = 0;
double lasty = 0;


mjtNum position_history = 0;
mjtNum previous_time = 0;

float_t ctrl_update_freq = 100;
mjtNum last_update = 0.0;
mjtNum ctrl;

void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods)
{
    
    if( act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE )
    {
        mj_resetData(m, d);
        mj_forward(m, d);
    }
}


void mouse_button(GLFWwindow* window, int button, int act, int mods)
{
    
    button_left =   (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
    button_right =  (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

    
    glfwGetCursorPos(window, &lastx, &lasty);
}


void mouse_move(GLFWwindow* window, double xpos, double ypos)
{
   
    if( !button_left && !button_middle && !button_right )
        return;

    
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    
    int width, height;
    glfwGetWindowSize(window, &width, &height);

   
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

    
    mjtMouse action;
    if( button_right )
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    else if( button_left )
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    else
        action = mjMOUSE_ZOOM;

   
    mjv_moveCamera(m, action, dx/height, dy/height, &scn, &cam);
}


void scroll(GLFWwindow* window, double xoffset, double yoffset)
{
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scn, &cam);
}


//****************************
void init_save_data()
{
  
   fprintf(fid,"t, ");
   fprintf(fid,"PE, KE, TE, ");
   fprintf(fid,"\n");
}

//***************************

void save_data(const mjModel* m, mjData* d)
{
  fprintf(fid,"%f, ",d->time);
  fprintf(fid,"%f, %f, %f",d->energy[0],d->energy[1],d->energy[0]+d->energy[1]);
  fprintf(fid,"\n");
}

//**************************
void mycontroller(const mjModel* m, mjData* d)
{
   mj_energyPos(m, d);
   mj_energyVel(m, d);
   
  //equation of motion
   
   double dense_M[4]={0};
   mj_fullM(m, dense_M, d->qM); 
   double M[2][2]={0};
   M[0][0]=dense_M[0]; 
   M[0][1]=dense_M[1]; 
   M[1][0]=dense_M[2]; 
   M[1][1]=dense_M[3]; 

   double qdd[2]={0};
   qdd[0]=d->qacc[0];
   qdd[1]=d->qacc[1];

   double b[2]={0};
   b[0]=d->qfrc_bias[0];
   b[1]=d->qfrc_bias[1];

   double lhs[2]={0};
   mju_mulMatVec(lhs,dense_M,qdd,2,2);
   lhs[0]=lhs[0]+b[0];
   lhs[1]=lhs[1]+b[1];

   d->qfrc_applied[0] = 0.1*b[0];
   d->qfrc_applied[1] = 0.5*b[1];

   double rhs[2]={0};
   rhs[0]=d->qfrc_applied[0];
   rhs[1]=d->qfrc_applied[1];

   //control (feedback linearization)
   double kp1=100, kp2=100, kv1=10, kv2=10;
   double qref1=0, qref2=0;

   double tau[2]={0};
   tau[0]=-kp1*(d->qpos[0]-qref1)-kv1*d->qvel[0];
   tau[1]=-kp1*(d->qpos[1]-qref2)-kv1*d->qvel[01];
   mju_mulMatVec(tau,dense_M,tau,2,2);
   tau[0]=tau[0]+b[0];
   tau[1]=tau[1]+b[1];

   d->qfrc_applied[0] = tau[0];
   d->qfrc_applied[1] = tau[1];

  
  if ( loop_index%data_frequency==0)
    {
      save_data(m,d);
    }
  loop_index = loop_index + 1;
}


//************************
// main function
int main(int argc, const char** argv)
{

    mj_activate("mjkey.txt");

    char xmlpath[100]={};
    char datapath[100]={};

    strcat(xmlpath,path);
    strcat(xmlpath,xmlfile);

    strcat(datapath,path);
    strcat(datapath,datafile);


    char error[1000] = "Could not load binary model";

    if( argc<2 )
        m = mj_loadXML(xmlpath, 0, error, 1000);

    else
        if( strlen(argv[1])>4 && !strcmp(argv[1]+strlen(argv[1])-4, ".mjb") )
            m = mj_loadModel(argv[1], 0);
        else
            m = mj_loadXML(argv[1], 0, error, 1000);
    if( !m )
        mju_error_s("Load model error: %s", error);


    d = mj_makeData(m);


    
    if( !glfwInit() )
        mju_error("Could not initialize GLFW");

    
    GLFWwindow* window = glfwCreateWindow(1244, 700, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);
    mjv_makeScene(m, &scn, 2000);                
    mjr_makeContext(m, &con, mjFONTSCALE_150);   

    
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    double arr_view[] = {90, -5, 6, 0, -1, 2};
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];

    
    mjcb_control = mycontroller;

    fid = fopen(datapath,"w");
    init_save_data();

    d->qpos[0] = -2.2;
    d->qpos[1] = 3.5;
    
    while( !glfwWindowShouldClose(window))
    {
        mjtNum simstart = d->time;
        while( d->time - simstart < 1.0/60.0 )
        {
            mj_step(m, d);
        }

        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

         
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);
        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    mjv_freeScene(&scn);
    mjr_freeContext(&con);

    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    #if defined(__APPLE__) || defined(_WIN32)
        glfwTerminate();
    #endif

    return 1;
}
