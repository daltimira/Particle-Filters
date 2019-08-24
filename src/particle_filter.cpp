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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int iParticle = 0; iParticle < num_particles; iParticle++) {
    Particle particle;
    particle.id = iParticle;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);
  }
  is_initialized = true;
}

std::default_random_engine gen;

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

  for (int iParticle = 0; iParticle < num_particles; iParticle++) {
    Particle & particle = particles[iParticle];

    if (fabs(yaw_rate) < TOLERANCE) {
      particle.x += velocity * delta_t * cos(particle.theta) + dist_x(gen);
		  particle.y += velocity * delta_t * sin(particle.theta) + dist_y(gen);
    } else {
      particle.x += (velocity / yaw_rate) * (sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta)) + dist_x(gen);
      particle.y += (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate*delta_t)) + dist_y(gen);
      particle.theta += yaw_rate * delta_t + dist_theta(gen);
    }

  }
   
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for (size_t i = 0; i<observations.size(); i++) {
    LandmarkObs & landMarkObs = observations[i];
    // now find the predicted measurement that is closest to this observation measurement
    double dist = 100000000;
    for (size_t j = 0; j<predicted.size(); j++) {
      const LandmarkObs & landMarkPred = predicted[j];
      double distPredObsSquared = (landMarkPred.x-landMarkObs.x)*(landMarkPred.x-landMarkObs.x) + (landMarkPred.y-landMarkObs.y)*(landMarkPred.y-landMarkObs.y);
      if (distPredObsSquared < dist) {
        landMarkObs.id = landMarkPred.id;
        dist = distPredObsSquared;
      }
    }
  }

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

  for (int iParticle = 0; iParticle < num_particles; iParticle++) {
    Particle & particle = particles[iParticle];
    particle.weight = 1.0;

    // build the predicted observations. These are the landmarks that we might expect to sense if we have the vehicle positioned in this particle and the sensor range
    vector<LandmarkObs> predicted;
    for (size_t iLandMark = 0; iLandMark < map_landmarks.landmark_list.size(); iLandMark++) {
      const Map::single_landmark_s & landMarkMap = map_landmarks.landmark_list[iLandMark];
      double distSquared = (landMarkMap.x_f-particle.x)*(landMarkMap.x_f-particle.x) + (landMarkMap.y_f-particle.y)*(landMarkMap.y_f-particle.y);
      if (distSquared <= sensor_range*sensor_range) {
        LandmarkObs landmark;
        landmark.id = landMarkMap.id_i;
        landmark.x = landMarkMap.x_f;
        landmark.y = landMarkMap.y_f;
        predicted.push_back(landmark);
      }
    }

    vector<LandmarkObs> observationsMapSpace;
    for (int iObservation = 0; iObservation<(int)observations.size(); iObservation++) {
      LandmarkObs landmark;
      landmark.id = observations[iObservation].id;
      landmark.x = particle.x + (cos(particle.theta)*observations[iObservation].x) - (sin(particle.theta)*observations[iObservation].y);
      landmark.y = particle.y + (sin(particle.theta)*observations[iObservation].x) + (cos(particle.theta)*observations[iObservation].y);
      observationsMapSpace.push_back(landmark);
    }

    dataAssociation(predicted, observationsMapSpace);

    for (int iObservation = 0; iObservation<(int)observationsMapSpace.size(); iObservation++) {
      double x_pred = 0, y_pred = 0;
      for (size_t iLandMark = 0; iLandMark < predicted.size(); iLandMark++) {
        const LandmarkObs & landMarkMap = predicted[iLandMark];
        if (landMarkMap.id == observationsMapSpace[iObservation].id) {
          x_pred = landMarkMap.x;
          y_pred = landMarkMap.y;
          break;
        }
      }
      double prob = multivariateGaussianProbability(observationsMapSpace[iObservation].x, observationsMapSpace[iObservation].y, x_pred, y_pred, std_landmark[0], std_landmark[1]);
      if (prob == 0.0) {
        prob = 0.00000000001;
      }
      particle.weight *= prob;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  // First let's guess a index uniformally [1,N]
  std::uniform_int_distribution<int>uniformDistributionInteger(0, num_particles - 1);
  std::uniform_real_distribution<double>uniformDistributionReal(0, 1);
  int index = uniformDistributionInteger(gen);
  double beta = 0;
  double maxW = 0;
  // get the max weight
  for (int iParticle = 0; iParticle < num_particles; iParticle++) {
    if (particles[iParticle].weight > maxW) {
      maxW = particles[iParticle].weight;
    }
  }

  std::vector<Particle> newSetParticles;
  for (int i = 0; i<num_particles; i++) {
    beta += uniformDistributionReal(gen)*2.0*maxW;
    while (beta > particles[index].weight) {
      beta -= particles[index].weight;
      index = (index + 1) % num_particles;
    }
    newSetParticles.push_back(particles[index]);
  }

  particles = newSetParticles;
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