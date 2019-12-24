/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::cout;
using std::endl;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  std::default_random_engine gen;
  Particle p;
  num_particles = 100;  // TODO: Set the number of particles
  is_initialized = 1;
  
  // This line creates a normal (Gaussian) distribution
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i = 0; i < num_particles; ++i) { 
    // TODO: Sample from these normal distributions like this: 
    // where "gen" is the random engine initialized earlier.
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1;
    particles.push_back(p);
  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  double angle, x, y, theta;
  if(yaw_rate)
  {
    for (int i = 0; i < num_particles; ++i) { 
      
      angle = particles[i].theta + (yaw_rate * delta_t); 
     
      x = velocity * (sin(angle) - sin(particles[i].theta)) / yaw_rate;
      y = velocity * (cos(particles[i].theta) - cos(angle)) / yaw_rate;
      theta = yaw_rate * delta_t;

      normal_distribution<double> dist_x(x, std_pos[0]);
      normal_distribution<double> dist_y(y, std_pos[1]);
      normal_distribution<double> dist_theta(theta, std_pos[2]);
      particles[i].x += dist_x(gen);
      particles[i].y += dist_y(gen);
      particles[i].theta += dist_theta(gen);
    }
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     LandmarkObs& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    
  double dist_min, _dist;
  int m_idx;

  dist_min=1000;

  for(int j=0;j<predicted.size();j++)
  {
    _dist=dist(predicted[j].x, predicted[j].y, observations.x, observations.y);
    if(_dist < dist_min)
    {
      dist_min = _dist;
      m_idx = j;
    }
  }
  observations.x=predicted[m_idx].x;
  observations.y=predicted[m_idx].y;
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */ 

      
  	// calculate normalization term
    double gauss_norm;
    gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  
  	std::vector <LandmarkObs> map_m;
  	for(int k=0;k<map_landmarks.landmark_list.size();k++)
  	{
    	LandmarkObs _map_m;
    	_map_m.id = map_landmarks.landmark_list[k].id_i;
    	_map_m.x = map_landmarks.landmark_list[k].x_f;
    	_map_m.y = map_landmarks.landmark_list[k].y_f;
    	map_m.push_back(_map_m);
  	}

   	for(int i=0;i<particles.size();i++)
   	{
        LandmarkObs obs_m, obs_n;  

        double exponent, weight=1;
      
   		for(int j=0;j<observations.size();j++)
   		{
     		if(dist(0, 0, observations[j].x, observations[j].y) <= sensor_range)
     		{
       			obs_m.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
     			obs_m.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
				obs_n = obs_m;

                dataAssociation(map_m, obs_m);
              
              	exponent = (pow(obs_m.x - obs_n.x, 2) / (2 * pow(std_landmark[0], 2))) + (pow(obs_m.y - obs_n.y, 2) / (2 * pow(std_landmark[1], 2)));
          		// calculate weight using normalization terms and exponent
        		weight *= gauss_norm * exp(-exponent);
            }
        }
      	particles[i].weight=weight;
   	}
  	
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   std::vector <double> dd;
   for(int i=0;i<particles.size();i++)
     dd.push_back(particles[i].weight);
   
   std::vector <Particle> P;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::discrete_distribution<> d(dd.begin(),dd.end());
   
   for(int i=0; i<particles.size(); ++i) {
       P.push_back(particles[d(gen)]);
   }
   particles = P;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}