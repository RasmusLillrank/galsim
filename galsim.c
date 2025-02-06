#include "graphics/graphics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
typedef struct vec vec_t;
typedef struct body body_t;

struct body{  
    double posX;
    double posY;
    double mass;
    double velX;
    double velY;
    double brightness;
};

struct vec {
  double x;
  double y;
}


const float circleRadius=0.005, circleColor=0;
const int windowWidth=800;

const double epsilon = 0.001;
body_t* read_file(int N, char* input_name) {
    FILE* file = fopen(input_name, "r");
    if (!file){
      printf("Error opening file!");
      return NULL;
    }

    body_t *data = (body_t*)malloc(N*6*sizeof(double));  // check it later- 
    if (!data){
      printf("File allocation error!");
      fclose(file);
      return NULL;
    }
    
    // size_t item_raad= fread(data,sizeof)
    if (fread(data, sizeof(body_t), N, file) != N) {
        printf("Incorrect number of particles!");
        free(data);
        fclose(file);
        return NULL;
    } else {
        fclose(file);
        return data;
    }
}

/*
calculate force
for j
 sum += -gravity * mass[i] * mass[j] * rel_dist(i,j) / (dist(i,j) + epsilon) ** 3
*/

vec_t force(int* r, int j, body_t* bodies, int N){
    
        double G = 100/N;
        double sum[2] = {0, 0};
     for(int i = 0; i < N; i++){

        double dist_x = bodies[j].posX - bodies[i].posX;
        double dist_y = bodies[j].posY - bodies[i].posY;
        double relative = sqrt((dist_x*dist_x)+(dist_y*dist_y));

        double dist_ux = (bodies[j].posX - bodies[i].posX);
        double dist_uy = (bodies[j].posY - bodies[i].posY);
        
        double distu[2] = {dist_ux, dist_uy};

        sum[0] += (bodies[i].mass)/(relative-epsilon)*distu[0];
        sum[1] += (bodies[i].mass)/(relative-epsilon)*distu[1];

     }

     
    vec_t results = {-G * bodies[j].mass * sum[0], -G * bodies[j].mass * sum[1]};
    return results;
}

void update_accelerate(body_t* bodies, int N){
  for (int i =1 ; i <N,i++){
    double sum[1];
    force(force_result)
  }
  
}

void update_velocity(){
  //previous the vecility
  // and then force/mass - new vector of the speed
  //update the velocity 
}


int main(int argc, char *argv[]) {
  if (argc != 6) {
    printf("Usage: %s N filename nsteps delta_t graphics", argv[0]);
  }
  const int N = atoi(argv[1]);
  body_t* bodies = read_file(N, argv[2]);
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  const int graphics = atoi(argv[5]);
  if(!bodies) { return; }

  float L=1, W=1;

  InitializeGraphics(argv[0],windowWidth,windowWidth);
  SetCAxes(0,1);
  //int N = 50;
  //body_t* bodies = read_file(N, "input_data/ellipse_N_00050.gal");

  printf("Hit q to quit.\n");
  while(!CheckForQuit()) {

    /* Call graphics routines. */
    ClearScreen();
    for(int i = 0; i < N; i++) {
        DrawCircle(bodies[i].posX, bodies[i].posY, L, W, bodies[i].mass*0.005, circleColor);
    }
    Refresh();
    /* Sleep a short while to avoid screen flickering. */
    usleep(3000);
  }
  FlushDisplay();
  CloseDisplay();
  return 0;
}