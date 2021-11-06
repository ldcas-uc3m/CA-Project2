using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>
#include <iomanip>


class Object{
    public:
        Object(double x, double y, double z, double mass){
            px = x;
            py = y;
            pz = z;
            vx = 0;
            vy = 0;
            vz = 0;
            m = mass;
            fx = 0;
            fy = 0;
            fz = 0;
        }
        double px;
        double py;
        double pz;
        double vx;
        double vy;
        double vz;
        double m;
        double fx;
        double fy;
        double fz;    
};

// constants
const double g = 6.674e-11;
const double COL_DISTANCE = 1;  // minimum colision distance

Object * populateWorld( int num_objects, int size_enclosure, int random_seed,  const char ** argcv){
      // distribution generation 
    random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{1e21, 1e15};
    
    gen64.seed(random_seed);  // introduce seed

    // memory alloc
    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);
    
    // init file
    ofstream inFile("init_config.txt", ofstream::out);  // open file
    
    inFile << fixed << setprecision(3) << argcv[4] << " " << argcv[5] << " " << argcv[1];
    inFile << endl;

    // populate

    double x, y, z, m;
    for (int i = 0; i < num_objects; i++){
        
        x = dis(gen64);
        y = dis(gen64);
        z = dis(gen64);
        m = d(gen64);
        universe[i] = Object(x, y, z, m);
        // write to file
        inFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m << endl;
    }

    inFile.close();

    return universe;
}
void updatePosition( Object *a, double time_step){

            // velocity calculation
            a->vx += (a->fx/a->m) * time_step;
            a->vy += (a->fy/a->m) * time_step;
            a->vz += (a->fz/a->m) * time_step;
            
            //reset force 
            a->fx = 0;
            a->fy = 0;
            a->fz = 0;
            
            // position calculation
            a->px += a->vx * time_step;
            a->py += a->vy * time_step;
            a->pz += a->vz * time_step;  

}
void mergeObjects(Object * a, Object * b){
/* ---OBJECT COLLISION--- */
// merge objects into a
a->m = a->m + b->m;
a->vx = a->vx + b->vx;
a->vy = a->vy + b->vy;
a->vz = a->vz + b->vz;
}
void forceBetweenTwoObjects(Object *a, Object *b, double distance, double dx, double dy, double dz){
                    
                    double dfx = (g * a->m * b->m * dx) / (distance*distance*distance);
                    double dfy = (g * a->m * b->m * dy) / (distance*distance*distance);
                    double dfz = (g * a->m * b->m * dz) / (distance*distance*distance);

                    // a forces
                    a->fx += dfx;
                    a->fy += dfy;
                    a->fz += dfz;

                    // b forces
                    b->fx -= dfx;
                    b->fy -= dfy;
                    b->fz -= dfz;
}

void reboundEffect(Object * a, int size_enclosure){


            if(a->px <= 0){
                a->px = 0;
                a->vx = - a->vx;
            } else if(a->px >= size_enclosure){
                a->px = size_enclosure;
                a->vx = - a->vx;
            }

            if(a->py <= 0){
                a->py = 0;
                a->vy = - a->vy;
            } else if(a->py >= size_enclosure){
                a->py = size_enclosure;
                a->vy = - a->vy;
            }

            if(a->pz <= 0){
                a->pz = 0;
                a->vz = - a->vz;
            } else if(a->pz >= size_enclosure){
                a->pz = size_enclosure;
                a->vz = - a->vz;
            }
            
}


void kernel(Object * universe,int num_objects, int size_enclosure, int num_iterations,double time_step, const char ** argcv){
     /* ---
    OUTPUT
    --- */
    
    ofstream outFile("final_config.txt");

    outFile << fixed << setprecision(3) << argcv[4] << " " << argcv[5] << " " << argcv[1] << endl;


    // extra vars
    int curr_objects = num_objects;
    bool *deleted = (bool *)calloc(num_objects, sizeof(bool)); // bytemap of objects -> if true, object is deleted

    /* ---
    KERNEL
    --- */
 
    for(int iteration = 0; iteration < num_iterations; iteration++){
        for(int i = 0; i < num_objects; i++){
            if(deleted[i]) continue;
            Object *a = &universe[i];

            for(int j = i + 1; j < num_objects; j++){
                if(deleted[j]) continue;
                
                Object *b = &universe[j];

                /* ---
                FORCE COMPUTATION
                --- */
                
                // distance
                double dx = b->px - a->px;
                double dy = b->py - a->py;
                double dz = b->pz - a->pz;
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

                if(distance <= COL_DISTANCE){
                    /* ---
                    OBJECT COLLISION
                    --- */

                    // merge objects into a
                    mergeObjects(a, b);

                    // del b
                    curr_objects--;
                    deleted[j] = true;

                    // force between a & b is 0
                } else{
                
               forceBetweenTwoObjects(a, b, distance, dx, dy, dz);
                }
            }
            /* ---
            UPDATE POSITION
            --- */
            updatePosition(a, time_step);
            /* ---
            REBOUND EFFECT
            --- */
            reboundEffect(a, size_enclosure);
            // print to output
            
            if((iteration == num_iterations - 1) || curr_objects == 1){  // final positions

            outFile << fixed << setprecision(3) << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
            << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
            << " " << universe[i].m << endl;
            }
        }
    }
    outFile.close();
}

int main(int argc, const char ** argcv){

    /* ---
    PARAMETERS
    --- */
     if(argc < 6){
        switch(argc){
            case 5:
                cerr  << "sim-soa invoked with " << argc << " parameters."
                 << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                 << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                 << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
                 << " time_step: ?"<< endl ;
                break;
            case 4:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                     << argcv[3] << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
            case 3:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
            case 2:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
            case 1:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: ?"
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
        }
        return -1;
    }else if(argc > 6){
        cerr << "sim-soa invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
             << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
             << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
             << " time_step: "<< argcv[5] << endl ;
        return -1;
    }


    // parameters init & casting
    const int num_objects = atoi(argcv[1]);
    const int num_iterations = atoi(argcv[2]);
    const int random_seed = atoi(argcv[3]);
    const double size_enclosure = atof(argcv[4]);
    const double time_step = atof(argcv[5]);

    // check correct parameters
    if(num_objects <= 0) {
        cerr << "Invalid number of object "<<endl << "sim-aos invoked with " << argc << " parameters."
              << endl << "Arguments: "<< endl << " num_objects: " << num_objects
              << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
              << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
              << " time_step: "<< time_step << endl ;
        return -2;
    }
    if(num_iterations <= 0){
        cerr << "Invalid number of iterations "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    if(random_seed<= 0){
        cerr << "Invalid seed "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    if(size_enclosure<= 0){
        cerr << "Invalid box size "<< endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    Object * universe =  populateWorld(num_objects, size_enclosure, random_seed, argcv);
    kernel(universe, num_objects , size_enclosure, num_iterations, time_step, argcv);
   
    return 0;
}
