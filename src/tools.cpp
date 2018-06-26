#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse = VectorXd(4);
	rmse.fill(0);

	//  estimations vector size should not be zero
	//  and should equal ground truth vector size
	if ((estimations.size() <= 0) || (estimations.size() != ground_truth.size())) {
	    cout << "Error: Invalid vector size.\n";
	    return rmse;
	}

	//accumulate squared residuals
	for(size_t i=0; i < estimations.size(); ++i){
        
        VectorXd residuals = estimations[i] - ground_truth[i];

	    VectorXd mean = residuals.array() * residuals.array();
    
        rmse += mean;
	}
		
	//calculate the squared root
	rmse = rmse/estimations.size();
	
	rmse = rmse.array().sqrt();
	
	return rmse;
}