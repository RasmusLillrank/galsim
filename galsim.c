#include "graphics/graphics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

typedef struct vec vec_t;
typedef struct body body_t;

// Struct used for reading and saving
struct body{  
    double posX;
    double posY;
    double mass;
    double velX;
    double velY;
    double brightness;
};


const float circleColor=0;
const int windowWidth=800;
const double epsilon = 0.001;

// Gets wall time
// Taken from Lab 7, Task 6
double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

// Reading the input .gal file
body_t* read_file(int N, char* input_name) {
    // If file is not found or readable, etc
    FILE* file = fopen(input_name, "r");
    if (!file){
      printf("Error opening file!\n");
      return NULL;
    }

    // Allocating data for file, if out of memory exit
    body_t *data = (body_t*)malloc(N*6*sizeof(double)); 
    if (!data){
      printf("File allocation error!\n");
      fclose(file);
      return NULL;
    }
    
    // Check if N > particles in file, exit
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

// Writing the result to result.gal
int write_file(int N, body_t* data) {
  // Create/open result.gal
  FILE* file = fopen("result.gal", "w");
  if(!file){
    printf("Error opening file!\n");
    return 0;
  }
  // Write data to result.gal
  fwrite((void*)data, sizeof(body_t), N, file);
  fclose(file);
  // Freeing the data
  free(data);
  return 0;
}


void update_bodies(int N, double dt, double* restrict posX, double* restrict posY, double* restrict mass, double* restrict velX, double* restrict velY){
  // Set the gravity constant
  const double G = (double)100/N;
  //
  for(int i = 0; i < N; i++) {
    double accX = 0;
    double accY = 0;
    for(int j = i+1; j < N; j++){
      // Distance related calculations
      const double dist_x = posX[i] - posX[j];
      const double dist_y = posY[i] - posY[j];
      const double relative = sqrt(dist_x*dist_x + dist_y*dist_y);
      const double rel_eps = relative+epsilon;
      const double dist_eps = rel_eps * rel_eps * rel_eps; //  relative+epsilon can be pre-computed

      // Calculate force and force vectors between i and j
      const double force = -G / dist_eps; // -G, mass can be simplified
      const double forceX = force * dist_x;
      const double forceY = force * dist_y;

      // Update velocities for all j particles
      velX[j] -= dt * (forceX * mass[i]);
      velY[j] -= dt * (forceY * mass[i]);

      // Update force sum for i particle
      accX += forceX * mass[j];
      accY += forceY * mass[j];
    }
    // Update velocities for i particle
    velX[i] += dt * accX;
    velY[i] += dt * accY;
  }
  // Update position for all particles
  for(int i = 0; i < N; i++) {
    posX[i]+= dt*velX[i];
    posY[i]+= dt*velY[i];
  }
}

// After loading data store it in separate arrays
void split_bodies(int N, body_t* bodies, double* posX, double* posY, double* mass, double* velX, double* velY) {
  for(int i = 0; i < N; i++) {
    posX[i] = bodies[i].posX;
    posY[i] = bodies[i].posY;
    mass[i] = bodies[i].mass;
    velX[i] = bodies[i].velX;
    velY[i] = bodies[i].velY;
  }
}

// After running the simulation store the updated data back into the body_t data array
void join_bodies(int N, body_t* bodies, double* posX, double* posY, double* mass, double* velX, double* velY) {
  for(int i = 0; i < N; i++) {
    bodies[i].posX = posX[i];
    bodies[i].posY = posY[i];
    bodies[i].velX = velX[i];
    bodies[i].velY = velY[i];
  }
  // Free the individual data arrays
  free(posX);
  free(posY);
  free(mass);
  free(velX);
  free(velY);
}

int main(int argc, char *argv[]) {
  if (argc != 6) {
    printf("Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
    return 0;
  }
  const int N = atoi(argv[1]);
  body_t* bodies = read_file(N, argv[2]);
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  const int graphics = atoi(argv[5]);
  if(!bodies) {return -1;} // Something went wrong reading the input file
  
  double* posX = (double *)malloc(sizeof(double)*N);
  double* posY = (double *)malloc(sizeof(double)*N);
  double* mass = (double *)malloc(sizeof(double)*N);
  double* velX = (double *)malloc(sizeof(double)*N);
  double* velY = (double *)malloc(sizeof(double)*N);
  split_bodies(N, bodies, posX, posY, mass, velX, velY);

  double time = get_wall_seconds();
  // Running without graphics
  if(!graphics) {
    for(int i = 0; i < nsteps; i++) {
      update_bodies(N, dt, posX, posY, mass, velX, velY);
    }
    time = get_wall_seconds()-time;
    printf("No graphics, time = %lf\n", time);
    join_bodies(N, bodies, posX, posY, mass, velX, velY);
    write_file(N, bodies);
    return 0;
  }

  float L=1, W=1;
  InitializeGraphics(argv[0],windowWidth,windowWidth);
  SetCAxes(0,1);
  printf("Hit q to quit.\n");
  int steps = 0;
  while(steps < nsteps && !CheckForQuit()) {
    steps++;
    /* Call graphics routines. */
    ClearScreen();
    update_bodies(N, dt, posX, posY, mass, velX, velY);
    for(int i = 0; i < N; i++) {
        DrawCircle(posX[i], posY[i], L, W, mass[i]*0.005, circleColor);
    }
    Refresh();
    /* Sleep a short while to avoid screen flickering. */
    usleep(3000);
  }
  FlushDisplay();
  CloseDisplay();
  time = get_wall_seconds() - time;
  printf("Graphics, time = %lf\n", time);
  join_bodies(N, bodies, posX, posY, mass, velX, velY);
  write_file(N, bodies);
  return 0;
}
