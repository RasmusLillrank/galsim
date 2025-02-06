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
};


const float circleColor=0;
const int windowWidth=800;

const double epsilon = 0.001;
body_t* read_file(int N, char* input_name) {
    FILE* file = fopen(input_name, "r");
    if (!file){
      printf("Error opening file!\n");
      return NULL;
    }

    body_t *data = (body_t*)malloc(N*6*sizeof(double));  // check it later- 
    if (!data){
      printf("File allocation error!\n");
      fclose(file);
      return NULL;
    }
    
    // size_t item_raad= fread(data,sizeof)
    if (fread(data, sizeof(body_t), N, file) != N) {
        printf("Incorrect number of particles!\n");
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

vec_t force(int j, body_t* bodies, int N){
    
  double G = 100/N;
  vec_t sum = {0, 0};
  for(int i = 0; i < N; i++){
    if (i == j) {continue;}

    double dist_x = bodies[j].posX - bodies[i].posX;
    double dist_y = bodies[j].posY - bodies[i].posY;
    double relative = sqrt((dist_x*dist_x)+(dist_y*dist_y)); 

    double dist_eps = (relative+epsilon) * (relative+epsilon) * (relative+epsilon);
    sum.x += (bodies[i].mass * relative)/dist_eps*dist_x;
    sum.y += (bodies[i].mass * relative)/dist_eps*dist_y;

  }
  vec_t results = {-G * bodies[j].mass * sum.x, -G * bodies[j].mass * sum.y};
  return results;
}

vec_t acceleration(int i, body_t* bodies, vec_t force) {
  return (vec_t) {force.x/bodies[i].mass, force.y/bodies[i].mass};
}

void update_bodies(body_t* bodies, int N, double dt){
  for(int i = 0; i < N; i++) {
    vec_t f = force(i, bodies, N);
    vec_t a = (i, bodies, f);
    bodies[i].velX += dt*a.x;
    bodies[i].velY += dt*a.y;
    bodies[i].posX += dt*bodies[i].velX;
    bodies[i].posY += dt*bodies[i].velY;
  }
}


int main(int argc, char *argv[]) {
  if (argc != 6) {
    printf("Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
    return 0;
  }
  const int N = atoi(argv[1]);
  printf("N: %i \n", N);
  body_t* bodies = read_file(N, argv[2]);
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  const int graphics = atoi(argv[5]);
  if(!bodies) {return -1;} // Something went wrong reading the input file

  float L=1, W=1;

  InitializeGraphics(argv[0],windowWidth,windowWidth);
  SetCAxes(0,1);
  printf("Hit q to quit.\n");
  int steps = 0;
  while(steps < nsteps && !CheckForQuit()) {
    steps++;
    /* Call graphics routines. */
    ClearScreen();
    update_bodies(bodies, N, dt);
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