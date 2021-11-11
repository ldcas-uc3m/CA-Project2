using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>
#include <iomanip>


class Universe{
    public:
        Universe(int num_objects, int size_enclosure){
            px = (double *)malloc(sizeof(double) * num_objects);
            py = (double *)malloc(sizeof(double) * num_objects);
            pz = (double *)malloc(sizeof(double) * num_objects);
            vx = (double *)calloc(num_objects, sizeof(double));
            vy = (double *)calloc(num_objects, sizeof(double));
            vz = (double *)calloc(num_objects, sizeof(double));
            m = (double *)malloc(sizeof(double) * num_objects);
            fx = (double *)calloc(num_objects, sizeof(double));
            fy = (double *)calloc(num_objects, sizeof(double));
            fz = (double *)calloc(num_objects, sizeof(double));
            objects = num_objects;
            size = size_enclosure;
        }
        double * px;
        double * py;
        double * pz;
        double * vx;
        double * vy;
        double * vz;
        double * m;
        double * fx;
        double * fy;
        double * fz;
        int objects;
        int size;
};


/* ---
CONSTANTS
--- */
const double G = 6.674e-11;
const double COL_DISTANCE = 1;  // minimum colision distance

Universe bigBang(double size_enclosure, int random_seed, int num_objects, int time_step){
    
    // distribution generation
    random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{1e21, 1e15};

    gen64.seed(random_seed);  // introduce seed

    // big bang
    Universe universe(num_objects, size_enclosure);
    
    // open input
    ofstream inFile("init_config.txt", ofstream::out);  // open file
    inFile << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;
    
    // populate
    for(int i = 0; i < num_objects; i++){
        universe.px[i] = dis(gen64);
        universe.py[i] = dis(gen64);
        universe.pz[i] = dis(gen64);
        universe.m[i] = d(gen64);

        inFile << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i]
        << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
        << " " << universe.m[i] << endl;
    }

    inFile.close();
    return universe;
}


 void mergeObjects(Universe universe, int i, int j){
    /*
    Merges two objects.
    */
    universe.vx[i] = universe.vx[i] + universe.vx[j];
    universe.vy[i] = universe.vy[i] + universe.vy[j];
    universe.vz[i] = universe.vz[i] + universe.vz[j];
    universe.m[i] = universe.m[i] + universe.m[j];
}


 void updatePosition(Universe universe, int i, double time_step){
    /*
    Updates the position of an object according to the time step
    and its current force.
    */

    // velocity calculation
    universe.vx[i] += (universe.fx[i]/universe.m[i]) * time_step;
    universe.vy[i] += (universe.fy[i]/universe.m[i]) * time_step;
    universe.vz[i] += (universe.fz[i]/universe.m[i]) * time_step;
            
    //reset force 
    universe.fx[i] = 0;
    universe.fy[i] = 0;
    universe.fz[i] = 0;

    // position calculation
    universe.px[i] += universe.vx[i] * time_step;
    universe.py[i] += universe.vy[i] * time_step;
    universe.pz[i] += universe.vz[i] * time_step;
}


void reboundEffect(Universe universe, int i, int size_enclosure){
    /*
    Checks for possible rebounds and updates positions accordingly.
    */

    if(universe.px[i] <= 0){
        universe.px[i] = 0;
        universe.vx[i] = - universe.vx[i];
    } else if(universe.px[i] >= size_enclosure){
        universe.px[i] = size_enclosure;
        universe.vx[i] = - universe.vx[i];
    }

    if(universe.py[i] <= 0){
        universe.py[i] = 0;
        universe.vy[i] = - universe.vy[i];
    } else if(universe.py[i] >= size_enclosure){
        universe.py[i] = size_enclosure;
        universe.vy[i] = - universe.vy[i];
    }

    if(universe.pz[i] <= 0){
        universe.pz[i] = 0;
        universe.vz[i] = - universe.vz[i];
    } else if(universe.pz[i] >= size_enclosure){
        universe.pz[i] = size_enclosure;
        universe.vz[i] = - universe.vz[i];
    }
}


bool forceComputation(Universe universe, int i, int j){
    /*
    Calculates the force between two objects and updates it.
    Returns true if there was a colision, false else.
    */

    // distance
    double dx = universe.px[j] - universe.px[i];
    double dy = universe.py[j] - universe.py[i];
    double dz = universe.pz[j] - universe.pz[i];
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

    if(distance <= COL_DISTANCE){
        // Object colision

        mergeObjects(universe, i, j);
        // force between a & b is 0

        return true;  // flag colision

    } else{
    
        double dfx = (G * universe.m[i] * universe.m[j] * dx) / (distance*distance*distance);
        double dfy = (G * universe.m[i] * universe.m[j] * dy) / (distance*distance*distance);
        double dfz = (G * universe.m[i] * universe.m[j] * dz) / (distance*distance*distance);

        // a forces
        universe.fx[i] += dfx;
        universe.fy[i] += dfy;
        universe.fz[i] += dfz;

        // b forces
        universe.fx[j] -= dfx;
        universe.fy[j] -= dfy;
        universe.fz[j] -= dfz;

        return false;
    }
}


int main(int argc, const char ** argcv){

    /* ---
    PARAMETERS
    --- */

    // check argc
    if(argc < 6){
        switch(argc){
            case 5:
                cerr  << "sim-soa invoked with " << argc << " parameters."
                 << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                 << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                 << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
                 << " time_step: ?"<< endl;
                break;
            case 4:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                     << argcv[3] << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                break;
            case 3:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                break;
            case 2:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                break;
            case 1:
                cerr << "sim-soa invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: ?"
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl;
                break;
        }
        return -1;
    } else if(argc > 6){
        cerr << "sim-soa invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
             << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
             << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
             << " time_step: "<< argcv[5] << endl;
        return -1;
    }

    // parameters init & casting
    const int num_objects = atoi(argcv[1]);
    const int num_iterations = atoi(argcv[2]);
    const int random_seed = atoi(argcv[3]);
    const double size_enclosure = atof(argcv[4]);
    const double time_step = atof(argcv[5]);

    // chech correct parameters
    if(num_objects <= 0){
        cerr << "Invalid number of object "<<endl << "sim-aos invoked with " << argc << " parameters."
              << endl << "Arguments: "<< endl << " num_objects: " << num_objects
              << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
              << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
              << " time_step: "<< time_step << endl;
        return -2;
    }
    if(num_iterations <= 0){
        cerr << "Invalid number of iterations "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl;
        return -2;
    }
    if(random_seed <= 0){
        cerr << "Invalid seed "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl;
        return -2;
    }
    if(size_enclosure <= 0){
        cerr << "Invalid box size "<< endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl;
        return -2;
    }

    Universe universe = bigBang(size_enclosure, random_seed, num_objects, time_step);

    // open output
    ofstream outFile("final_config.txt");
    outFile << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    // extra vars
    int curr_objects = num_objects;
    bool *deleted = (bool *)calloc(num_objects, sizeof(bool));  // bytemap of objects -> if true, object is deleted


    /* ---
    KERNEL
    --- */
    
    for(int iteration = 0; iteration < num_iterations; iteration++){
        for(int i = 0; i < num_objects; i++){
            if(deleted[i]) continue;

            for(int j = i + 1; j < num_objects; j++){
                if(deleted[j]) continue;
                
                if(forceComputation(universe, i, j)){
                    // delete b
                    curr_objects--;
                    deleted[j] = true;
                }
            }

            updatePosition(universe, i, time_step);
            reboundEffect(universe, i, size_enclosure);
    
            if((iteration == num_iterations - 1) ||  curr_objects == 1){  // final positions
                outFile << fixed << setprecision(3)  << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i] 
                << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
                << " " << universe.m[i] << endl;
            }
        }
    }

    outFile.close();
    return 0;
}
