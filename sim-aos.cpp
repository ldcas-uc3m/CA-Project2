using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <omp.h>

struct Force{
    double fx;
    double fy;
    double fz;
};

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

/* ---
CONSTANTS
--- */
const double g = 6.674e-11;
const double COL_DISTANCE = 1;  // minimum colision distance
const int THREAD_NUM = 4 ;


int checkArguments(int argc, const char ** argcv){
    /*
    Checks the arguments of the program and returns the
    error codes.
    */

    // check number of parameters
    if(argc < 6){
        switch(argc){
            case 5:
                cerr  << "sim-soa invoked with " << argc << " parameters."
                 << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                 << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                 << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
                 << " time_step: ?"<< endl;
                return -1;
            case 4:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                     << argcv[3] << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                return -1;
            case 3:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                return -1;
            case 2:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                return -1;
            case 1:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: ?"
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                return -1;
        }
    } else if(argc > 6){
        cerr << "sim-soa invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
             << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
             << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
             << " time_step: "<< argcv[5] << endl;
        return -1;
    }
    const int num_objects = atoi(argcv[1]);
    const int num_iterations = atoi(argcv[2]);
    const int random_seed = atoi(argcv[3]);
    const double size_enclosure = atof(argcv[4]);
    const double time_step = atof(argcv[5]);

    // check correct parameters
    if(num_objects <= 0){ // num_objects
        cerr << "Invalid number of object "<<endl << "sim-aos invoked with " << argc << " parameters."
              << endl << "Arguments: "<< endl << " num_objects: " << num_objects
              << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
              << random_seed << endl << " size_enclosure: " << size_enclosure << endl
              << " time_step: "<< time_step << endl;
        return -2;
    }
    if(num_iterations <= 0){ // num_iterations
        cerr << "Invalid number of iterations "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " << size_enclosure << endl
             << " time_step: "<< time_step << endl;
        return -2;
    }
    if(random_seed <= 0){ // random_seed
        cerr << "Invalid seed "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " << size_enclosure << endl
             << " time_step: "<< time_step << endl;
        return -2;
    }
    if(size_enclosure <= 0){ // size_enclosure
        cerr << "Invalid box size "<< endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " << size_enclosure << endl
             << " time_step: "<< time_step << endl;
        return -2;
    }
    return 0;
}



Object * bigBang(int num_objects, int size_enclosure, int random_seed, int time_step){
    /*
    Creates the universe and populates it. Also fills init_config.txt.
    */
    
    // distribution generation 
    random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{1e21, 1e15};
    gen64.seed(random_seed);

    // memory alloc
    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);
    
    // init file
    ofstream inFile("init_config.txt", ofstream::out);  // open file
    inFile << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    // populate
    for(int i = 0; i < num_objects; i++){
        
        double x = dis(gen64);
        double y = dis(gen64);
        double z = dis(gen64);
        double m = d(gen64);
        
        universe[i] = Object(x, y, z, m);

        // write to file
        inFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m << endl;
    }

    inFile.close();

    return universe;
}


void updatePosition(Object *a, double time_step){
    /*
    Updates the position of an object according to the time step
    and its current force.
    */

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


void mergeObjects(Object *a, Object *b, int curr_objects, bool deleted[], int j){
    /*
    Merges two objects into one.
    */
    #pragma omp atomic
    // merge objects into a
    a->m = a->m + b->m;
    a->vx = a->vx + b->vx;
    a->vy = a->vy + b->vy;
    a->vz = a->vz + b->vz;

    //delete(&b);
    curr_objects--;
    deleted[j] = true;
}


Force forceComputation(Object *a, Object *b, int curr_objects, bool *deleted, int j){
    /*
    Calculates the force between two objects and updates it.
    Returns false if there was a colision, true else.
    */

    // distance
    const double dx = b->px - a->px; 
    const double dy = b->py - a->py;
    const double dz = b->pz - a->pz;
    const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    /* note these vars are constant bc they will be out of scope outside the function */

    struct Force force{0, 0, 0};

    if(distance <= COL_DISTANCE){
        // Object colision

        mergeObjects(a, b, curr_objects, deleted, j);

        // force between a & b is 0

    } else{
        const double dfx = (g * a->m * b->m * dx) / (distance*distance*distance);
        const double dfy = (g * a->m * b->m * dy) / (distance*distance*distance);
        const double dfz = (g * a->m * b->m * dz) / (distance*distance*distance);

        force.fx = dfx;
        force.fy = dfy;
        force.fz = dfz;
        
        // b force
        b->fx -= dfx;
        b->fy -= dfy;
        b->fz -= dfz;
    }
    return force;
}


void reboundEffect(Object *a, int size_enclosure){
    /*
    Checks for possible rebounds and updates positions accordingly.
    */
    
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


int main(int argc, const char ** argcv){

    /* ---
    PARAMETERS
    --- */

    const int aerror = checkArguments(argc, argcv);

    if(aerror < 0){
        return aerror;
    }

    // parameters init & casting
    const int num_objects = atoi(argcv[1]);
    const int num_iterations = atoi(argcv[2]);
    const int random_seed = atoi(argcv[3]);
    const double size_enclosure = atof(argcv[4]);
    const double time_step = atof(argcv[5]);


    /* ---
    INITIALIZATION
    --- */

    Object * universe = bigBang(num_objects, size_enclosure, random_seed, time_step);
    
    // prepare output
    ofstream outFile("final_config.txt");
    outFile << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    // extra vars
    int curr_objects = num_objects;
    bool *deleted = (bool *)calloc(num_objects, sizeof(bool)); // bytemap of objects -> if true, object is deleted


    /* ---
    KERNEL
    --- */
    for(int iteration = 0; iteration < num_iterations; iteration++){

        if(curr_objects != 1){
            for(int i = 0; i < num_objects; i++){
                double dfx = 0; double dfy = 0; double dfz = 0;
                if(!deleted[i]){
                    Object *a = &universe[i];
                    #pragma omp parallel for reduction(+:dfx,dfy,dfz)
                    for(int j = i + 1; j < num_objects; j++){
                        
                        if(!deleted[j]){
                            Object *b = &universe[j];
                            struct Force force = forceComputation(a, b, curr_objects, deleted, j);
                            dfx += force.fx;
                            dfy += force.fy;
                            dfz += force.fz;
                        }
                    }
                    // take into account the pre-saved forces on a
                    a->fx += dfx;
                    a->fy += dfy;
                    a->fz += dfz;
                    updatePosition(a, time_step);
                    reboundEffect(a, size_enclosure);
                    

                    if((iteration == num_iterations - 1) || curr_objects == 1){  // final positions
                        // print to output
                        outFile << fixed << setprecision(3) << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
                        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
                        << " " << universe[i].m << endl;
                    }
                }                                
            }       
        }           
    }
    outFile.close();
    return 0;
}
