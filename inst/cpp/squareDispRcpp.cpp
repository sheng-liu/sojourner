/*
squareDispRcpp.cpp

Wu Lab, Johns Hopkins University
(Referenced original codebase squareDisp.R algorithm from Sheng Liu)

Author: Sun Jay Yoo
Date: May 22, 2017
*/

#include <armadillo>
#include <vector>
#include <string>
#include <RcppArmadillo.h>
#include <stdexcept>
using namespace Rcpp;

//Required comment headers for Rcpp Armadillo (DO NOT DELETE FOLLOWING COMMENTS)
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export()]]
List squareDispRcpp(arma::mat track, int dt = 1, double resolution = 0.107){

	//Throw error if dt is greater than the track length - 1
	if (dt >= track.n_rows){
	  std::string error = 
	    "\ntrack length:\t" + std::to_string(track.n_rows) + "\n" + 
	    "dt:     \t" + std::to_string(dt) + "\n" +
	    "Time interval (dt) greater than track length-1";
		throw std::invalid_argument(error);
	}

	//Create a vector of matrices with algorithmically accurate dimensions according to dt time steps
	std::vector<arma::mat> trackOut;
	trackOut.push_back(arma::mat(track.n_rows/dt + track.n_rows % dt, 7));
	for (int i = 1; i < dt; i++){
		trackOut.push_back(arma::mat(track.n_rows / dt, 7));
	}

	//Fill matrices in the vector with square displacement calculations for each index of the input track
	int c = 0; //Counter variable for each matrix (ex. counter alternates index 0 and 1 if dt = 2)
	int j = 0; //Index for each matrix

	for (int i = 0; i < track.n_rows; i++){

			//Copying the original coordinates and index per track input index
			trackOut[c](j, 0) = track (i, 0);
			trackOut[c](j, 1) = track (i, 1);
			trackOut[c](j, 2) = track (i, 2);
			trackOut[c](j, 3) = i + 1;

			//Displacement data null if coordinate at previous time step doesn/t exist
			if (i < dt){
				trackOut[c](j, 4) = arma::datum::nan;
				trackOut[c](j, 5) = arma::datum::nan;
				trackOut[c](j, 6) = arma::datum::nan;

			//Calculate displacement and square displacement
			} else {
				double dx = (track (i, 0)-track (i - dt, 0))*resolution;
				double dy = (track (i, 1)-track (i - dt, 1))*resolution;
				trackOut[c](j, 4) = dx * dx + dy * dy;
				trackOut[c](j, 5) = dx;
				trackOut[c](j, 6) = dy;

			}

			//Alternate through indexes
			c++;
			if (c == dt){
				c = 0;
				j++;
			}
	}

	//Create Rcpp::List type.
	//(List type not used in algorithm as .push_back() for List is extremely memory ineffecient
	//and element access operations difficult in C++ for non-native List type)
	List tracklist(dt);
	
	//Cast back each matrix from the vector into a R Numeric Matrix and set appropriate dimension names.
	//Place ordered Numeric Matrices into List.
	for (int i = 0; i < dt; i++){
	  NumericMatrix m = wrap(trackOut[i]);
	  colnames(m) = CharacterVector::create("x", "y", "z", "index", "square.disp", "dx", "dy");
	  arma::colvec r = trackOut[i].col(3);
	  CharacterVector v = as<CharacterVector>(wrap(r));
	  rownames(m) = v;
	  
		tracklist[i] = m;
	}
	return tracklist;
}
