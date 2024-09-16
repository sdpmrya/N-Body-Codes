#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <omp.h>
#include "gnuplot-iostream.h"

// Constants
const double G =0.01;       // Gravitational constant
const int N = 512;          // Number of particles
const double M = 1.0;         // Mass of one particle

struct Particle 
{
    double x, y, z;
    double vx, vy, vz;
    double mass;
};

double plummer(double r) {
    double b = 1.0, pi=M_PI;
    return 3/(4*pi*pow(b,3.0))*pow(1 + (r*r)/(b*b), -5/2);
}

void randomDist(Particle* particles) {
    using namespace std;

    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i=0; i<N; i++) {
    particles[i].x = distribution(generator);
    particles[i].y = distribution(generator);
    particles[i].z = distribution(generator);
    particles[i].vx = 0.0;
    particles[i].vy = 0.0;
    particles[i].vz = 0.0;
    particles[i].mass = M;
    }
}

void initializeparticles(Particle* particles) {
    using namespace std;
    //random number generator
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, 4.0);

    for (int i=0; i<N; i++) {
        double theta=distribution(generator)*M_PI, phi=distribution(generator)*2*M_PI, r=distribution(generator);
        double r1 = plummer(r);
        particles[i].x = r1*sin(theta)*cos(phi) + 0.5;
        particles[i].y = r1*sin(theta)*sin(phi) + 0.5;
        particles[i].z = r1*cos(theta) + 0.5;

        double x = particles[i].x - 0.5;
        double y = particles[i].y - 0.5;
        double z = particles[i].z - 0.5;

        double radius = std::sqrt(x*x + y*y + z*z);

        double vel = 2.0;
        particles[i].vx = radius*vel*distribution(generator);
        particles[i].vy = radius*vel*distribution(generator);
        particles[i].vz = radius*vel*distribution(generator);
        particles[i].mass = M;
    }
}

void filesave(Particle* particles) {
    std::ofstream outFile("ParticleDist.txt");
    if (outFile.is_open()) {
        for (int i=0; i<N; ++i){
            for (int j=0; j<3; ++j){
                if (j==0) outFile << particles[i].x << " ";
                if (j==1) outFile << particles[i].y << " ";
                if (j==2) outFile << particles[i].z << " ";
            }
            outFile << "\n";
        }
        outFile.close();
    }
}

void computeDirectForces(Particle* particles, double dt) {
    double eps = 0.001; //softening
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            // initializing force
            double force = 0.0;
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;

            // Apply periodic boundary conditions
            if (dx > 0.5) dx -= 1.0;
            if (dx < -0.5) dx += 1.0;
            if (dy > 0.5) dy -= 1.0;
            if (dy < -0.5) dy += 1.0;
            if (dz > 0.5) dz -= 1.0;
            if (dz < -0.5) dz += 1.0;

            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            force = (G * particles[i].mass * particles[j].mass) / (dist*std::sqrt(dist) + eps);
            
            double fx = force * dx;
            double fy = force * dy;
            double fz = force * dz;

            //std::cout << dx << " , " << dy << " , " << dz <<  "\n" ;
            //std::cout << force << "\n" ;

            //if (i == j) fx = 0.0, fy = 0.0, fz = 0.0;

            //Velocity update  
            particles[i].vx += fx * dt;
            particles[i].vy += fy * dt;
            particles[i].vz += fz * dt;

            particles[j].vx -= fx * dt;
            particles[j].vy -= fy * dt;
            particles[j].vz -= fz * dt;
        }
    }
}

void positionupdate(Particle* particles, double dt) {
    for (int i=0; i<N; i++) {
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;

        //Periodic boundary consditions
        particles[i].x = fmod(particles[i].x + 1.0, 1.0);
        particles[i].y = fmod(particles[i].y + 1.0, 1.0);
        particles[i].z = fmod(particles[i].z + 1.0, 1.0);
    }
}

int main() {
    Particle* particles = new Particle[N];
    initializeparticles(particles);

    //for saving the data:
    //std::ofstream outFile("Simul.txt")

    double TotTime = 1.0; 
    const double dt = 0.0001;     // Time step
    int num_steps = (int) TotTime/dt;

    std::ofstream outFile("Simul.txt");

    for (int i=0; i< num_steps; i++) {
        computeDirectForces(particles, dt);
        positionupdate(particles, dt);
        for (int j=0; j<N; j++) {
            outFile << particles[j].x << " " << particles[j].y << " " << particles[j].z << "\n";
            }
            outFile << "\n";
        }
    outFile.close();
    return 0;
}