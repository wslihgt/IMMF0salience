/*
  copyright (C) 2010 Jean-Louis Durrieu
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Remember to use a different guard symbol in each header!
#ifndef _F0_SALIENCE_H_
#define _F0_SALIENCE_H_

#include <vamp-sdk/Plugin.h>
#include <complex>
#include <fftw3.h>
#include <cstdlib>
#include <cmath>
#include "cppblas.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

#define log2(x) (log((double) x)/log(2.0))
#define max(x,y) ((x<y)?y:x)

using namespace std;
using std::string;
using std::endl;
using std::cerr;

class F0Salience : public Vamp::Plugin
{
 public:
  F0Salience(float inputSampleRate);
  virtual ~F0Salience();
  
  string getIdentifier() const;
  string getName() const;
  string getDescription() const;
  string getMaker() const;
  int getPluginVersion() const;
  string getCopyright() const;
  
  InputDomain getInputDomain() const;
  size_t getPreferredBlockSize() const;
  size_t getPreferredStepSize() const;
  size_t getMinChannelCount() const;
  size_t getMaxChannelCount() const;
  
  ParameterList getParameterDescriptors() const;
  float getParameter(string identifier) const;
  void setParameter(string identifier, float value);
  
  ProgramList getPrograms() const;
  string getCurrentProgram() const;
  void selectProgram(string name);
  
  OutputList getOutputDescriptors() const;
  
  bool initialise(size_t channels, size_t stepSize, size_t blockSize);
  void reset();
  
  FeatureSet process(const float *const *inputBuffers,
		     Vamp::RealTime timestamp);
  
  FeatureSet getRemainingFeatures();
  
 protected:
  void genHannBasis();
  void hanning(size_t length, float *window); 
  size_t genBaseWF0(float Fm, float FM, float inputSamplingRate, 
		    size_t NFT, size_t stepnote, float opencoef, 
		    float *WF0, float *f0Table); 
  void genODGDspec(float f0, float inputSamplingRate, 
		   size_t blockSize, size_t NFT, float opencoef, 
		   complex<double> *in, complex<double> *odgdspec, 
		   fftw_plan p);
  void updateHmatrices();
  
  
  string m_currentProgram;
  size_t m_stepSize;
  size_t m_blockSize;
  float m_inputSampleRate;
  size_t m_nbChannels; // added here, but with min/max = 1 maybe dealt in SV directly
  size_t m_numberBins;
  size_t m_nbAccumulatedFrames; // Activates the computation as soon as there are enough frames in the matrix
  size_t m_numberFrames; // for later implementation: stacking the frames, working on matrices, instead of vector.
  
  // some variables for the design of the basis matrix for the source WF0
  float m_fmin;
  float m_fmax;
  size_t m_stepnote;
  float m_opencoef;
  size_t m_numberOfF0;
  size_t m_numberOfGamma;
  float m_overlapGamma;
  size_t m_maxIterations;
  float m_thresholdEnergy;
  FeatureSet m_featureSet;
  
  // basis matrices
  float *m_WF0;    // for the source part
  float *m_HF0;
  float *m_f0Table;
  float *m_WGAMMA; // for the filter part
  float *m_HGAMMA;
  float *m_SX; // matrix or vector of input to be processed
  int m_computeBasis;
};



#endif
