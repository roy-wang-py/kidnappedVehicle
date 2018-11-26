/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

#define DEFAULT_PARTICLES_NUM 100
#define DBL_MAX 1.7976931348623158e+308 /* max value */
default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);    
    
    num_particles  = DEFAULT_PARTICLES_NUM;

    if(is_initialized)
    {
        //inited alreay 
        return;
    }

    for(int i=0; i<num_particles; i++)
    {
        Particle particle_tmp;
        //wmemset(&particle_tmp, 0, sizeof(particle_tmp));
        particle_tmp.id = i;
        particle_tmp.x = dist_x(gen);
        particle_tmp.y = dist_y(gen);
        particle_tmp.theta = dist_theta(gen);
        particle_tmp.weight = 1.0;
        
        particles.push_back(particle_tmp);
    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);


    for(int i=0; i<num_particles; i++)
    {
        double theta = particles[i].theta;

        if (fabs(yaw_rate) > 0.00001) 
        { 
            //theta changed
            double theta_dt = yaw_rate*delta_t;
            particles[i].x += velocity / yaw_rate * (sin(theta + theta_dt) - sin(theta));
            particles[i].y += velocity / yaw_rate * (cos(theta) - cos(theta + theta_dt));
            particles[i].theta += theta_dt;

        } 
        else 
        {
            //theta was the same
            particles[i].x += velocity * delta_t * cos(theta);
            particles[i].y += velocity * delta_t * sin(theta);

        }

        // add noise.
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
     }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    //num of observations
    unsigned int ob_num = observations.size();
    //num of predicteds
    unsigned int pre_num = predicted.size();

    for(unsigned int i=0; i < ob_num; i++)
    {
        double min_distance = DBL_MAX;
        double distance = 0.0;
        int land_id = -1;
        

        for(unsigned int j=0; j<pre_num; j++)
        {
            //find the instace whitch is the closest to  observations[i]
            /*
            if(observations[i].id)
            {
                continue;
            }
            */
                
            //double x_distance = observations[i].x - predicted[j].x;
            //double y_distance = observations[i].y - predicted[j].y;

            //distance = (x_distance*x_distance) + (y_distance*y_distance);
            distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

            if(distance < min_distance)
            {
                min_distance = distance;
                land_id = predicted[j].id;
            }
        }
        
        observations[i].id = land_id;
    }
    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    for(int i=0; i <num_particles ; i++)
    {
        
        double p_x = particles[i].x;
        double p_y = particles[i].y;;
        double p_theta = particles[i].theta;
        std::vector<LandmarkObs> predicted;

        //predict: find the landmarks those are in range
        for(unsigned int j=0; j<map_landmarks.landmark_list.size(); j++)
        {

            double m_x = map_landmarks.landmark_list[j].x_f;
            double m_y = map_landmarks.landmark_list[j].y_f;
            int m_id = map_landmarks.landmark_list[j].id_i;

            double distance = dist(p_x, p_y, m_x, m_y);
            if(distance > sensor_range)
            {
                //could not find this landmark,continue
                continue;
            }
                        
            LandmarkObs find_landmark;
            find_landmark.x = m_x;
            find_landmark.y = m_y;
            find_landmark.id = m_id;
                
            
            predicted.push_back(find_landmark);
            
        }

        //from vehicle's corodinate to map's coordinate system
        std::vector<LandmarkObs> observations_map;
        for(unsigned int j=0; j<observations.size() ; j++)
        {
            LandmarkObs ob_map;
            ob_map.x = cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y + p_x;
            ob_map.y = sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y + p_y;
            ob_map.id = observations[j].id;
            
            observations_map.push_back(ob_map);
        }

        //association
        dataAssociation(predicted, observations_map);

        //calc weights
        double sig_x = std_landmark[0];
        double sig_y = std_landmark[1];
        double x_obs,y_obs,mu_x,mu_y;


        //reset weight
        particles[i].weight = 1.0;   

        for(unsigned int j=0; j<observations_map.size(); j++)
        {
            unsigned int k=0;
                        
            x_obs = observations_map[j].x;
            y_obs = observations_map[j].y;
            
            for(k=0; k<map_landmarks.landmark_list.size(); k++)
            {
                if(observations_map[j].id == map_landmarks.landmark_list[k].id_i)
                {
                    break;
                }
            }

            if(k >= map_landmarks.landmark_list.size())
            {
                //map land mark not found
                continue;
            }

            mu_x = map_landmarks.landmark_list[k].x_f;
            mu_y = map_landmarks.landmark_list[k].y_f;

            //calculate normalization term
            double gauss_norm= (1.0/(2.0 * M_PI* sig_x * sig_y));
            double exponent= pow(x_obs-mu_x, 2.0)/(2 * pow(sig_x, 2.0)) + pow(y_obs-mu_y, 2.0)/(2 * pow(sig_y, 2.0));
            //calculate weight using normalization terms and exponent
            double weight= gauss_norm * exp(-exponent);       
            particles[i].weight *= weight;
        }

        //weights[i] = particles[i].weight;
        //weights.push_back(particles[i].weight);
    }
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> new_particles;

    weights.clear();
    for(int i=0; i <num_particles ; i++)
    {
        weights.push_back(particles[i].weight);
    }

    //double max_weight = *max_element(weights.begin(),weights.end());  

    // create distributions
    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::discrete_distribution<int> d_idx(weights.begin(), weights.end());

    // Generating index.
    int index;

    // the wheel
    for(int i = 0; i < num_particles; i++) 
    {
        index =  d_idx(gen);
        new_particles.push_back(particles[index]);
    }

    particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
