#include "graphics/graphics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <string.h>

typedef struct vec vec_t;
typedef struct body body_t;
typedef struct args args_t;

// Initializing global pointers
double* restrict posX;
double* restrict posY;
double* restrict mass;
double* restrict velX;
double* restrict velY;

// Struct used for reading and saving
struct body{  
    double posX;
    double posY;
    double mass;
    double velX;
    double velY;
    double brightness;
};

struct args{
  int N;
  int nthreads;
  int id;
  double dt;
  double* dvelX;
  double* dvelY;
}; 


const float circleColor=0;
const int windowWidth=800;
const double epsilon = 0.001;

pthread_mutex_t lock1;
pthread_mutex_t lock2;

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


void* update_bodies(void* args){

  args_t* args2 = (args_t*)args;

  const int N = args2->N;
  const int nthreads = args2->nthreads;
  const int id = args2->id;
  const double dt = args2->dt;
  double* restrict dvelX = args2->dvelX;
  double* restrict dvelY = args2->dvelY;
  // Initialize delta velocity to 0
  memset(args2->dvelX, 0, N*sizeof(double));
  memset(args2->dvelY, 0, N*sizeof(double));

  
  // Set the gravity constant
  const double G = (double)100/N;
  #pragma omp simd
  for(int i = id; i < N-1; i+= nthreads) {
    double accX = 0;
    double accY = 0;
    const double px = posX[i];
    const double py = posY[i];
  #pragma omp simd
    for(int j = i+1; j < N; j++){
      // Distance related calculations
      const double dist_x = px - posX[j];
      const double dist_y = py - posY[j];
      const double relative = sqrt(dist_x*dist_x + dist_y*dist_y);
      const double rel_eps = relative+epsilon;
      const double dist_eps = rel_eps * rel_eps * rel_eps; //  relative+epsilon can be pre-computed

      // Calculate force and force vectors between i and j
      const double force = -G / dist_eps; // -G, mass can be simplified
      const double forceX = force * dist_x;
      const double forceY = force * dist_y;

      // Update velocities for all j particles
      dvelX[j] -= dt * (forceX * mass[i]);
      dvelY[j] -= dt * (forceY * mass[i]);

      // Update force sum for i particle
      accX += forceX * mass[j];
      accY += forceY * mass[j];
    }
    // Update velocities for i particle
    dvelX[i] += dt * accX;
    dvelY[i] += dt * accY;
  }
  pthread_mutex_lock(&lock1);
  for(int k = 0; k < N; k++) {
      velX[k] += dvelX[k];
  }
  pthread_mutex_unlock(&lock1);

  pthread_mutex_lock(&lock2);
  for(int k = 0; k < N; k++) {
      velY[k] += dvelY[k];
  }
  pthread_mutex_unlock(&lock2);
  // Update position for all particles 
  return NULL;
}

void* update_positions(void *args) {
  args_t* args2 = (args_t*)args;

  const int N = args2->N;
  const int nthreads = args2->nthreads;
  const int id = args2->id;
  const double dt = args2->dt;
  double* restrict dvelX = args2->dvelX;
  double* restrict dvelY = args2->dvelY;

  for(int i = id; i < N; i+=nthreads) {
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
  if (argc != 7) {
    printf("Usage: %s N filename nsteps delta_t graphics nthreads\n", argv[0]);
    return 0;
  }
  const int N = atoi(argv[1]);
  body_t* bodies = read_file(N, argv[2]);
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  const int graphics = atoi(argv[5]);

  const int nthreads = atoi(argv[6]);
  printf("nthreads = %d\n", nthreads);

  pthread_mutex_init(&lock1, NULL);
  pthread_mutex_init(&lock2, NULL);
  
  if(!bodies) {return -1;} // Something went wrong reading the input file
  
  posX = (double *)malloc(sizeof(double)*N);
  posY = (double *)malloc(sizeof(double)*N);
  mass = (double *)malloc(sizeof(double)*N);
  velX = (double *)malloc(sizeof(double)*N);
  velY = (double *)malloc(sizeof(double)*N);
  split_bodies(N, bodies, posX, posY, mass, velX, velY);

  // Create threads
  pthread_t threads[nthreads];
  args_t args[nthreads];
  double* dvelX[nthreads];
  double* dvelY[nthreads];

  // Initialize arguments
  for(int i = 0; i < nthreads; i++) {
    args[i].dt = dt;
    args[i].N = N;
    args[i].id = i;
    args[i].nthreads = nthreads;
    args[i].dvelX = (double*)malloc(sizeof(double)*N);
    args[i].dvelY = (double*)malloc(sizeof(double)*N);
  }

  double time = get_wall_seconds();
  // Running without graphics
  if(!graphics) {
    for(int i = 0; i < nsteps; i++) {
      for(int j = 0; j < nthreads; j++) { 
        // Launch threads
        pthread_create(&threads[j], NULL, update_bodies, (void*)&args[j]); 
      }
      // Join threads
      for(int j = 0; j < nthreads; j++) {
        pthread_join(threads[j], NULL);
      }
      // Update positions
      // Parallellize
      for(int j = 0; j < nthreads; j++) { 
        // Launch threads
        pthread_create(&threads[j], NULL, update_positions, (void*)&args[j]); 
      }
      for(int j = 0; j < nthreads; j++) {
        pthread_join(threads[j], NULL);
      }
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

  // TODO: Fix graphics
  int steps = 0;
  while(steps < nsteps && !CheckForQuit()) {
    steps++;
    /* Call graphics routines. */
    ClearScreen();
    
    //update_bodies(N, dt, posX, posY, mass, velX, velY);
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
