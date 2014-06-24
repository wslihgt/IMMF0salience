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

#include "f0Salience.h"

#define db(x) (10.0 * log(x) / log(10.0))
#define invdb(y) (exp((y) / 10.0 * log(10.0)))

// including some pre-computed matrices
// the following included files are added, each with a 
// precomputed matrix (using python/numpy version of this 
// program). 
//
//     wf0_fmin_fmax_inputSampleRate_blockSize_stepnote_opencoef
//
// TODO: some differences between these versions and matrices directly 
// computed may be seen - this should be fixed, if one finds the time
// to do so...
#include "wf0_50_600_16000_512_16_05.h"
#include "wf0_50_600_44100_2048_16_05.h" 
#include "wf0_100_800_44100_2048_8_05.h" 
#include "wf0_100_800_16000_2048_16_05.h"
#include "wf0_100_800_44100_2048_16_05.h"
#include "wf0_50_2500_.h" // for Fs=44100, stepnote=20 and Ot=0.5

F0Salience::F0Salience(float inputSampleRate) :
  Plugin(inputSampleRate),
  m_currentProgram("default"),
  m_stepSize(0),
  m_blockSize(0),
  m_nbAccumulatedFrames(0),
  m_numberFrames(10),
  m_fmin(100.0),
  m_fmax(800.0),
  m_stepnote(8.0),
  m_opencoef(0.5),
  m_numberOfF0(0),
  m_numberOfGamma(50),
  m_overlapGamma(.75),
  m_maxIterations(50),
  m_thresholdEnergy(30.0),
  m_featureSet(FeatureSet()),
  m_computeBasis(0)
{
  m_inputSampleRate = inputSampleRate;
}

F0Salience::~F0Salience()
{
  if(m_computeBasis==1)
    {
      m_computeBasis = 0;
      delete [] m_WF0;
      delete [] m_HF0;
      delete [] m_f0Table;
      delete [] m_WGAMMA;
      delete [] m_HGAMMA;
    }
}

string
F0Salience::getIdentifier() const
{
  return "f0salience";
}

string
F0Salience::getName() const
{
  return "F0 Salience";
}

string
F0Salience::getDescription() const
{
  return "This plug-in produces a representation exhibiting the salience of the fundamental frequency (F0). Based on Durrieu et al, ``A Musically Motivated Representation For Pitch Estimation And Musical Source Separation'', 2010. Plug-in may freeze during initialization: this is a normal behavior, the program needs to compute big matrices before starting to process the data. Apologies for this (necessary) computational issue! ";
}

string
F0Salience::getMaker() const
{
  return "J.-L. Durrieu";
}

int
F0Salience::getPluginVersion() const
{
  // Increment this each time you release a version that behaves
  // differently from the previous one
  
  // version 1: initial version, implementing IMM 
  // version 2: adding an output descriptor = spectral envelope.
  // version 3: adding pre-computed matrices and programs.
  // version 4: (rev 265) adding control on dynamics of output
  return 4;
}

string
F0Salience::getCopyright() const
{
  // This function is not ideally named.  It does not necessarily
  // need to say who made the plugin -- getMaker does that -- but it
  // should indicate the terms under which it is distributed.  For
  // example, "Copyright (year). All Rights Reserved", or "GPL"
  return "GPL";
}

F0Salience::InputDomain
F0Salience::getInputDomain() const
{
  return FrequencyDomain;
}

size_t
F0Salience::getPreferredBlockSize() const
{
  return 2048; // 0 means "I can handle any block size"
}

size_t 
F0Salience::getPreferredStepSize() const
{
  return 1024; // 0 means "anything sensible"; in practice this
              // means the same as the block size for TimeDomain
              // plugins, or half of it for FrequencyDomain plugins
}

size_t
F0Salience::getMinChannelCount() const
{
  return 1;
}

size_t
F0Salience::getMaxChannelCount() const
{
  return 1;
}

F0Salience::ParameterList
F0Salience::getParameterDescriptors() const
{
  ParameterList list;
  
  // If the plugin has no adjustable parameters, return an empty
  // list here (and there's no need to provide implementations of
  // getParameter and setParameter in that case either).

  ParameterDescriptor numberFrames;
  numberFrames.identifier = "numberProcFrames";
  numberFrames.name = "Number of processed frames";
  numberFrames.description = "Number of frames to be processed as one chunk";
  numberFrames.unit = "Frame(s)";
  numberFrames.minValue = 1; 
  numberFrames.maxValue = 5000;
  numberFrames.defaultValue = 1;
  numberFrames.isQuantized = true;
  numberFrames.quantizeStep = 1.0;
  list.push_back(numberFrames);
  
  ParameterDescriptor fmin;
  fmin.identifier = "Fmin";
  fmin.name = "Minimum F0";
  fmin.description = "Minimum F0";
  fmin.unit = "Hz";
  fmin.minValue = 50.0;
  fmin.maxValue = 2000.0;
  fmin.defaultValue = 100.0;
  fmin.isQuantized = false;
  list.push_back(fmin);
  
  ParameterDescriptor fmax;
  fmax.identifier = "FMax";
  fmax.name = "Maximum F0";
  fmax.description = "Maximum F0";
  fmax.unit = "Hz";
  fmax.minValue = 50.0;
  fmax.maxValue = 6000.0;
  fmax.defaultValue = 800.0;
  fmax.isQuantized = false;
  list.push_back(fmax);
  
  ParameterDescriptor stepnote;
  stepnote.identifier = "stepnote";
  stepnote.name = "Step between two notes";
  stepnote.description = "This is the step between two consecutive F0s, such that there are ``stepnote'' F0s per semitone";
  stepnote.unit = "'per semitone'";
  stepnote.minValue = 1.0;
  stepnote.maxValue = 200.0;
  stepnote.defaultValue = 8.0;
  stepnote.isQuantized = true;
  stepnote.quantizeStep = 1.0;
  list.push_back(stepnote);
  
  ParameterDescriptor opencoef;
  opencoef.identifier = "Ot";
  opencoef.name = "Open Coefficient";
  opencoef.description = "Open coefficient for the KLGLOTT88 glottal source model";
  opencoef.unit = "";
  opencoef.minValue = 0.0;
  opencoef.maxValue = 1.0;
  opencoef.defaultValue = 0.5;
  opencoef.isQuantized = false;
  list.push_back(opencoef);
  
  ParameterDescriptor nGamma;
  nGamma.identifier = "numberOfGamma";
  nGamma.name = "Number of elementary filters";
  nGamma.description = "The number of elementary filters for the filter decomposition.";
  nGamma.unit = "";
  nGamma.minValue = 4.0;
  nGamma.maxValue = 1000.0;
  nGamma.defaultValue = 50.0;
  nGamma.isQuantized = true;
  nGamma.quantizeStep = 1.0;
  list.push_back(nGamma);
  
  ParameterDescriptor overlapGamma;
  overlapGamma.identifier = "overlapGamma";
  overlapGamma.name = "Overlap of elementary filters";
  overlapGamma.description = "Overlapping rate of the elementary filters.";
  overlapGamma.unit = "";
  overlapGamma.minValue = .25;
  overlapGamma.maxValue = .98;
  overlapGamma.defaultValue = 0.75;
  overlapGamma.isQuantized = false;
  list.push_back(overlapGamma);
  
  ParameterDescriptor maxIter;
  maxIter.identifier = "maximumIterations";
  maxIter.name = "Max. Iterations";
  maxIter.description = "The maximum number of iterations for the estimation algorithm.";
  maxIter.unit = "";
  maxIter.minValue = 5;
  maxIter.maxValue = 1000;
  maxIter.defaultValue = 50;
  maxIter.isQuantized = true;
  maxIter.quantizeStep = 1;
  list.push_back(maxIter);

  ParameterDescriptor thresDisp;
  thresDisp.identifier = "thresholdEnergyDisplay";
  thresDisp.name = "threshold energy display";
  thresDisp.description = "threshold for energy display (in db)";
  thresDisp.unit = "dB";
  thresDisp.minValue = 0.0;
  thresDisp.maxValue = 500.0;
  thresDisp.defaultValue = 30.0;
  thresDisp.isQuantized = false;
  list.push_back(thresDisp);

  return list;
}

float
F0Salience::getParameter(string identifier) const
{
  if (identifier == "Fmin") {
    return m_fmin;
  } else if (identifier == "FMax") {
    return m_fmax;
  } else if (identifier == "stepnote"){
    return m_stepnote;
  } else if (identifier == "Ot"){
    return m_opencoef;
  } else if (identifier == "numberOfGamma"){
    return m_numberOfGamma;
  } else if (identifier == "overlapGamma"){
    return m_overlapGamma;
  } else if (identifier == "maximumIterations"){
    return m_maxIterations;
  } else if (identifier == "thresholdEnergyDisplay"){
    return m_thresholdEnergy;
  } else if (identifier == "numberProcFrames"){
    return m_numberFrames;
  } 
    
  return 0;
}

void
F0Salience::setParameter(string identifier, float value) 
{
  if (identifier == "Fmin") {
    m_fmin = value;
  } else if (identifier == "FMax") {
    m_fmax = value;
  } else if (identifier == "stepnote"){
    m_stepnote = value;
  } else if (identifier == "Ot"){
    m_opencoef = value;
  } else if (identifier == "numberOfGamma"){
    m_numberOfGamma = value;
  } else if (identifier == "overlapGamma"){
    m_overlapGamma = value;
  } else if (identifier == "maximumIterations"){
    m_maxIterations = value;
  } else if (identifier == "thresholdEnergyDisplay"){
    m_thresholdEnergy = value;
  } else if (identifier == "numberProcFrames"){
    m_numberFrames = value;
  } 
}

F0Salience::ProgramList
F0Salience::getPrograms() const
{
  ProgramList list;

  // If you have no programs, return an empty list (or simply don't
  // implement this function or getCurrentProgram/selectProgram)
  
  list.push_back("default"); 
  list.push_back("low and high pitched, fine scale"); 
  list.push_back("speech, 16kHz"); 
  list.push_back("speech, 44.1kHz"); 
  list.push_back("average popular music, coarse scale"); 
  list.push_back("average popular music, finer scale"); 

  return list;
}

string
F0Salience::getCurrentProgram() const
{
  return m_currentProgram; // no programs
}

void
F0Salience::selectProgram(string name)
{
  m_currentProgram = name;
  if (name == "default")
    {
      m_fmin = 100;
      m_fmax = 800;
      m_stepnote = 8;
      m_opencoef = 0.5;
    }
  if (name == "speech, 16kHz")
    {
      m_fmin = 50;
      m_fmax = 600;
      m_stepnote = 16;
      m_opencoef = 0.5;
      m_blockSize = 512;
    }
  if (name == "speech, 44.1kHz")
    {
      m_fmin = 50;
      m_fmax = 600;
      m_stepnote = 16;
      m_opencoef = 0.5;
    }
  if (name == "low and high pitched, fine scale")
    {
      m_fmin = 50;
      m_fmax = 2500;
      m_stepnote = 20;
      m_opencoef = 0.5;
    }
  if (name == "average popular music, coarse scale")
    {
      m_fmin = 100;
      m_fmax = 800;
      m_stepnote = 8;
      m_opencoef = 0.5;
    }
  if (name == "average popular music, finer scale")
    {
      m_fmin = 100;
      m_fmax = 800;
      m_stepnote = 16;
      m_opencoef = 0.5;
    }
}

F0Salience::OutputList
F0Salience::getOutputDescriptors() const
{
  OutputList list;

  // See OutputDescriptor documentation for the possibilities here.
  // Every plugin must have at least one output.

  vector<string> f0values;
  size_t f;
  char f0String[10];
  for (f=0; f<m_numberOfF0; f++)
    {
      sprintf(f0String, "%d Hz", (int)m_f0Table[f]);
      f0values.push_back(f0String);
    }

  OutputDescriptor d;
  d.identifier = "f0salience";
  d.name = "Salience of F0s.";
  d.description = "This representation should emphasize the salience of the different F0s in the signal.";
  d.unit = "";
  d.hasFixedBinCount = true;
  if (m_blockSize == 0) {
    // Just so as not to return "1". This is the bin count that
    // would result from a block size of 1024, which is a likely
    // default -- but the host should always set the block size
    // before querying the bin count for certain.
    d.binCount = 513;
  } else {
    d.binCount = m_numberOfF0; 
  } 
  d.binNames = f0values;
  d.hasKnownExtents = false;
  d.isQuantized = false;
  d.sampleType = OutputDescriptor::OneSamplePerStep; // VariableSampleRate;// OneSamplePerStep;
  d.hasDuration = false;
  list.push_back(d);
  
  // TODO: add an output descriptor for the spectral envelope
  OutputDescriptor specEnv;
  float my_inputSampleRate;
  specEnv.identifier = "spectralEnvelope";
  specEnv.name = "Spectral Envelope";
  specEnv.description = "The computed spectral envelope is drawn.";
  specEnv.unit = "";
  specEnv.hasFixedBinCount = true;
  if (m_blockSize == 0) {
    // Just so as not to return "1". This is the bin count that
    // would result from a block size of 1024, which is a likely
    // default -- but the host should always set the block size
    // before querying the bin count for certain.
    specEnv.binCount = 513;
    my_inputSampleRate = 44100.0;
  } else {
    specEnv.binCount = m_blockSize/2+1; 
    my_inputSampleRate = m_inputSampleRate;
  } 
  vector<string> freqNames;
  char fString[100];
  for (f=0; f<m_blockSize/2+1; f++)
    {
      sprintf(fString, "%d Hz", (int)(f / ((specEnv.binCount-1.0)*2.0) * my_inputSampleRate));
      freqNames.push_back(fString);
    }
  specEnv.binNames = freqNames;
  specEnv.hasKnownExtents = false;
  specEnv.isQuantized = false;
  specEnv.sampleType = OutputDescriptor::OneSamplePerStep;//VariableSampleRate;//
  specEnv.hasDuration = false;
  list.push_back(specEnv);
  
  return list;
}

bool
F0Salience::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
  cout << "in initialise" << endl;
  if (channels < getMinChannelCount() ||
      channels > getMaxChannelCount()) return false;
  m_nbChannels = channels;
  m_blockSize = blockSize;
  m_stepSize = stepSize;
  m_numberBins = m_blockSize / 2 + 1;
  m_numberOfF0 = (int) floor(12*m_stepnote*log2(m_fmax/m_fmin));
  
  if(m_computeBasis==1)
    {
      m_computeBasis = 0;
      delete [] m_WF0;
      delete [] m_HF0;
      delete [] m_f0Table;
      delete [] m_WGAMMA;
      delete [] m_HGAMMA;
    }
  
  m_WF0 = new float [m_numberOfF0 * m_numberBins];
  m_f0Table = new float [m_numberOfF0];
  m_WGAMMA = new float [m_numberOfGamma * m_numberBins];
  m_HF0 = new float [m_numberFrames * m_numberOfF0];
  m_HGAMMA = new float [m_numberFrames * m_numberOfGamma];
 
  cout << "generate WGAMMA" << endl;
  genHannBasis();
  cout << "generate WF0" << endl;
  m_numberOfF0 = genBaseWF0(m_fmin, m_fmax, m_inputSampleRate, 
			    m_blockSize, m_stepnote, m_opencoef, 
			    m_WF0, m_f0Table);
  cout << "leaving initiate" << endl;
  /*
    int n;
    for ( n=0; n<m_numberOfF0; n++)
    {cout << m_f0Table[n] << " ";} cout << endl;
  */
  m_computeBasis = 1;
  
  return true;
}

void
F0Salience::reset()
{
  cout << "--> reset()" << endl;
  delete [] m_WF0;
  delete [] m_HF0;
  delete [] m_f0Table;
  delete [] m_WGAMMA;
  delete [] m_HGAMMA;
  m_computeBasis = 0;
}

void 
F0Salience::hanning(size_t length, float *window)
{
  size_t i=0;
  for(i=0;i<length;i++)
    {
      window[i] = 0.5 - 0.5 * cos(2 * M_PI * i / (float)length);
    }
}

void
F0Salience::genHannBasis()
{
  size_t f, n; // variables for the loops
  size_t lengthGamma = (int) floor(m_numberBins / 
                                ((1.0 - m_overlapGamma) * 
				 (m_numberOfGamma - 1) + 
				 1 - 2.0 * m_overlapGamma));
  // making an (oversized) prototype window, to make the assignment easier.
  float *bigWindow;
  bigWindow = new float [3*m_numberBins];
  
  // the center for each basis function is one of the bins of the
  // frequency range
  int center;// size_t center; // change this to int center; //?
  
  for(f=0; f<3*m_numberBins; f++)
    {
      bigWindow[f] = 0.0;
    }
  // hann function for the elementary window, 
  // TODO: could be parameterized to adapt the user's need... 
  hanning(lengthGamma, &bigWindow[(int) (1.5*m_numberBins-lengthGamma/2)]);
  
  for(n=0; n<m_numberOfGamma; n++)
    {
      // The centers are uniformly and linearly spaced:
      center = (- 1 / (1 - m_overlapGamma) + 1 + n) * 
	(1 - m_overlapGamma) * lengthGamma + lengthGamma / 2;
      // filling in the elementary functions.
      for(f=0; f<m_numberBins; f++)
	{
	  m_WGAMMA[n*m_numberBins+f] = 
	    bigWindow[(int)(1.5*m_numberBins-center+f)];
	}
    }
  // freeing memory
  delete [] bigWindow;
}

size_t
F0Salience::genBaseWF0(float Fm, float FM, float inputSamplingRate, 
		       size_t NFT, size_t stepnote, float opencoef, 
		       float *WF0, float *f0Table)
/*
  Generates the WF0 source spectra basis. 

  It is somehow an independent function, in which the arguments are
  all given, instead of relying on the class attributes.

  If possible, tries to load the data from internally saved and 
  precalculated bases.
  
 */
{
  size_t N_notes = (int) floor(12 * stepnote * log2(FM / Fm));
  size_t ff=0, F=0, Nf=0;
  float eps=1.0e-10;
  Nf = m_numberBins; // renaming, for ease of use
  
  // If the parameters are "known", directly load the matrix
  if (Fm==50 && FM==600 && inputSamplingRate==16000 && 
      NFT==512 && stepnote==16 && opencoef==0.5)
    {
      // N_notes should be 690 - it actually is less, but that s fine.
      // cout << "N_notes should be 690, is " << N_notes << endl;
      // cout << "Total number in WF0 should be 177330, is ";
      // cout << N_notes*Nf << endl;
      cout << "Loading from pre-computed Matrices" << endl;
      for(ff=0;ff<N_notes-1;ff++)
	{
	  f0Table[ff] = F0Table_50_600_16000_512_16_05[ff];
	  for(F=0;F<Nf;F++)
	    {
	      WF0[ff*Nf+F] = WF0_50_600_16000_512_16_05[ff*Nf+F] + eps;
	    }
	}
      
      // last element actually noise:
      for(F=0;F<Nf;F++)
	{
	  WF0[(N_notes-1)*Nf+F] = 1.0;
	}
      
      return N_notes;
    }
  if (Fm==50 && FM==600 && inputSamplingRate==44100 && 
      NFT==2048 && stepnote==16 && opencoef==0.5)
    {
      // N_notes should be 690
      // cout << "N_notes should be 690, is " << N_notes << endl;
      // cout << "Total number in WF0 should be 707250, is ";
      // cout << N_notes*Nf << endl;
      cout << "Loading from pre-computed Matrices" << endl;
      for(ff=0;ff<N_notes-1;ff++)
	{
	  f0Table[ff] = F0Table_50_600_44100_2048_16_05[ff];
	  for(F=0;F<Nf;F++)
	    {
	      WF0[ff*Nf+F] = WF0_50_600_44100_2048_16_05[ff*Nf+F] + eps;
	    }
	}
      
      // last element actually noise:
      for(F=0;F<Nf;F++)
	{
	  WF0[(N_notes-1)*Nf+F] = 1.0;
	}
      
      return N_notes;
    }
  if (Fm==100 && FM==800 && inputSamplingRate==16000 && 
      NFT==2048 && stepnote==16 && opencoef==0.5)
    {
      // N_notes should be 577
      // cout << "N_notes should be 577, is " << N_notes << endl;
      // cout << "Total number in WF0 should be 591425, is ";
      // cout << N_notes*Nf << endl;
      cout << "Loading from pre-computed Matrices" << endl;
      for(ff=0;ff<N_notes-1;ff++)
	{
	  f0Table[ff] = F0Table_100_800_16000_2048_16_05[ff];
	  for(F=0;F<Nf;F++)
	    {
	      WF0[ff*Nf+F] = WF0_100_800_16000_2048_16_05[ff*Nf+F] + eps;
	    }
	}
      
      // last element actually noise:
      for(F=0;F<Nf;F++)
	{
	  WF0[(N_notes-1)*Nf+F] = 1.0;
	}
      
      return N_notes;
    }
  if (Fm==100 && FM==800 && inputSamplingRate==44100 && 
      NFT==2048 && stepnote==8 && opencoef==0.5)
    {
      // N_notes should be 289
      // cout << "N_notes should be 289, is " << N_notes << endl;
      // cout << "Total number in WF0 should be 296225, is ";
      // cout << N_notes*Nf << endl;
      cout << "Loading from pre-computed Matrices" << endl;
      for(ff=0;ff<N_notes-1;ff++)
	{
	  f0Table[ff] = F0Table_100_800_44100_2048_8_05[ff];
	  for(F=0;F<Nf;F++)
	    {
	      WF0[ff*Nf+F] = WF0_100_800_44100_2048_8_05[ff*Nf+F] + eps;
	    }
	}
      
      // last element actually noise:
      for(F=0;F<Nf;F++)
	{
	  WF0[(N_notes-1)*Nf+F] = 1.0;
	}
      
      return N_notes;
    }
  if (Fm==100 && FM==800 && inputSamplingRate==44100 && 
      NFT==2048 && stepnote==16 && opencoef==0.5)
    {
      // N_notes should be 577
      // cout << "N_notes should be 577, is " << N_notes << endl;
      // cout << "Total number in WF0 should be 591425, is ";
      // cout << N_notes*Nf << endl;
      cout << "Loading from pre-computed Matrices" << endl;
      for(ff=0;ff<N_notes-1;ff++)
	{
	  f0Table[ff] = F0Table_100_800_44100_2048_16_05[ff];
	  for(F=0;F<Nf;F++)
	    {
	      WF0[ff*Nf+F] = WF0_100_800_44100_2048_16_05[ff*Nf+F] + eps;
	    }
	}
      
      // last element actually noise:
      for(F=0;F<Nf;F++)
	{
	  WF0[(N_notes-1)*Nf+F] = 1.0;
	}
      
      return N_notes;
    }
  if (Fm==50 && FM==2500 && inputSamplingRate==44100 && 
      NFT==2048 && stepnote==20 && opencoef==0.5)
    {
      // N_notes should be 1356
      // cout << "N_notes should be 1356, is " << N_notes << endl;
      // cout << "Total number in WF0 should be 1389900, is ";
      // cout << N_notes*Nf << endl; 
      cout << "Loading from pre-computed Matrices" << endl;
      for(ff=0;ff<N_notes-1;ff++)
	{
	  f0Table[ff] = F0Table_50_2500_44100_2048_20_05[ff];
	  for(F=0;F<Nf;F++)
	    {
	      WF0[ff*Nf+F] = WF0_50_2500_44100_2048_20_05[ff*Nf+F] + eps;
	    }
	}
      
      // last element actually noise:
      for(F=0;F<Nf;F++)
	{
	  WF0[(N_notes-1)*Nf+F] = 1.0;
	}
      return N_notes;
    }
  cout << "Computing directly the matrices..." << endl;
  float norm;
  complex<double> *in, *odgdspec; 
  // from my former C implementation: used to be fftw_complex 
  fftw_plan p;
  
  in = new complex<double> [NFT];
  odgdspec = new complex<double> [NFT];
  
  // calling the fftw function, not sure it does what I expect...
  // the info comes from: 
  // http://www.fftw.org/fftw3_doc/Complex-numbers.html
  p = fftw_plan_dft_1d(NFT, reinterpret_cast<fftw_complex*>(in), 
		       reinterpret_cast<fftw_complex*>(odgdspec), 
		       FFTW_FORWARD, FFTW_MEASURE);
  
  f0Table[0] = Fm;
  for(ff=1;ff<N_notes;ff++)
    {
      f0Table[ff] = f0Table[ff-1] * exp(log(2)/(12*stepnote));
    }
  
  for(ff=0;ff<N_notes;ff++)
    {
      // cout << "    computing ODGD " << ff << "..." << endl;
      norm = 0;
      genODGDspec(f0Table[ff], inputSamplingRate, NFT, 
		  NFT, opencoef, in, odgdspec, p);
      // cout << "    additional computations (abs, norm)..." << endl;
      for(F=0;F<Nf;F++)
	{
	  WF0[ff*Nf+F] = abs(odgdspec[F])*abs(odgdspec[F]);
	  if(norm<WF0[ff*Nf+F])
	    {
	      norm = WF0[ff*Nf+F];
	    }
	}
      // cout << "    normalizing..." << endl;
      for(F=0;F<Nf;F++)
	{
	  WF0[ff*Nf+F] = WF0[ff*Nf+F]/norm + eps;
	  /* if (ff == 100){cout << setiosflags(ios::fixed) << setprecision(10)<<WF0[ff*Nf+F] << ", ";} */
	}
    }
  // last element actually noise:
  for(F=0;F<Nf;F++)
    {
      WF0[(N_notes-1)*Nf+F] = 1.0;
    }
  
  
  fftw_destroy_plan(p);
  //fftw_free(in); fftw_free(odgdspec);
  delete [] in;
  delete [] odgdspec;
  return N_notes;
}

void
F0Salience::genODGDspec(float f0, float inputSamplingRate, size_t blockSize, 
			size_t NFT, float opencoef,  complex<double> *in, 
			 complex<double> *odgdspec, fftw_plan p)
/*
  Note that in, odgdspec and p should be linked _beforehand_, before 
  calling this method, in order to be able to run 
  fftw_execute(p);
  Typical call to fftw should look like that:
  p = fftw_plan_dft_1d(NFT, in, odgdspec, FFTW_FORWARD, FFTW_MEASURE);
*/
{
  size_t ff=0,t=0;
  complex<double> ctemp(0.0, 0.0), *comp; // temporary quantities
  complex<double> I(0.0, 1.0); 
  float *W, *odgd; // analysis window, glottal waveform 
  
  size_t L = blockSize; 
  if(L>NFT){ L = NFT; } // L is the window length.
  
  W = new float [L];
  odgd = new float [L];
  hanning(L, W); // using Hann window for analysis. 
  /* 
     note that one should use the same as the one 
     chosen within the program.
  */
  
  size_t PartielMax = (size_t) floor((inputSamplingRate/2.0)/f0);
  
  complex<double> tempVal = (complex<double>)2*I*
    (complex<double>)M_PI/((complex<double>)inputSamplingRate);
  //float *f_partiels;
  //f_partiels = new float [PartielMax];
  complex<double> *f_partiels;
  f_partiels = new complex<double> [PartielMax];
  
  // cout << "        computing comp" << endl;
  comp = new complex<double> [PartielMax];
  for(ff=0;ff<PartielMax;ff++)
    {
      f_partiels[ff] = tempVal * complex<double> ((ff+1) * f0);
      comp[ff] = (ff+1) * M_PI * opencoef * 2 * I;
    }

  // cout << "        computing amplitudes" << endl;
  complex<double> *amplitudes;
  amplitudes = new complex<double> [PartielMax];
  for(ff=0;ff<PartielMax;ff++)
    {
      // TODO: the type casting in this equation is a problem...
      // hopefully not much impact on the actual result?
      amplitudes[ff] = (complex<double>)(f0 * 27 / 4.0) * 
	(exp(-comp[ff]) + (((complex<double>)2*
			    ((complex<double>)1+
			     (complex<double>)2*exp(-comp[ff])))/
			   (comp[ff]))-(((complex<double>)6 * 
					 ((complex<double>)1 - 
					  exp(-comp[ff])))/
					(comp[ff]*comp[ff]))) / comp[ff];
    }
  // cout << "        computing odgd" << endl;
  for(t=0;t<L;t++)
    {
      // cout << "            odgd: time n" << t << endl;
      ctemp=(complex<double>)0*I;
      for(ff=0;ff<PartielMax;ff++)
	{
	  // cout << "                odgd: partial n" << ff << endl;
	  ctemp = ctemp + amplitudes[ff]*
	    exp((complex<double>)f_partiels[ff]*
		(complex<double>)t);
	}
      odgd[t] = W[t] * 2 * (float) real(ctemp);
    }
  
  for(t=0;t<L;t++)
    {
      in[t] = ((complex<double>)odgd[t]);
    }
  for(t=L;t<NFT;t++)
    {
      in[t] = ((complex<double>)0.0);
    }
  // cout << "        computing the FFT" << endl;
  fftw_execute(p);
  // cout << "        end of FFT, freeing mem" << endl;
  // freeing memory
  delete [] odgd;
  delete [] W;
  delete [] comp;
  delete [] amplitudes;
  delete [] f_partiels;
}

void
F0Salience::updateHmatrices()
{
  // just for debugging... 
  // cout << "we are here updateHmat" << endl;
  float eps = 1.0e-20;// 0.00000000000000000000000001;
  size_t niter, n, f;
  // all the following matrices are useful for use with cblas
  float *sF0Hat, *sPhiHat, *tmpNum, *tmpDen, 
    *numeraHF0, *denomiHF0, *numeraHG, *denomiHG, *maxGamma;

  size_t nbMus=20;
  float *WM, *HM, *sM; // accompaniment part
  float *numHM, *denHM, *numWM, *denWM;
  
  float *sXHat;

  sF0Hat = new float [m_numberFrames * m_numberBins];
  sPhiHat = new float [m_numberFrames * m_numberBins];
  tmpNum = new float [m_numberFrames * m_numberBins];
  tmpDen = new float [m_numberFrames * m_numberBins];
  numeraHF0 = new float [m_numberFrames * m_numberOfF0];
  denomiHF0 = new float [m_numberFrames * m_numberOfF0];
  numeraHG = new float [m_numberFrames * m_numberOfGamma];
  denomiHG = new float [m_numberFrames * m_numberOfGamma];
  maxGamma = new float [m_numberFrames];

  WM = new float [nbMus * m_numberBins];
  HM = new float [m_numberFrames * nbMus];
  sM = new float [m_numberFrames * m_numberBins];
  // initialize the accompaniment matrices
  for(f=0; f<nbMus; f++)
    {
      for (n=0; n<m_numberFrames; n++)
	{
	  HM[n*nbMus+f] = 0.5;
	}
      for (n=0; n<m_numberBins; n++)
	{
	  WM[f*m_numberBins+n] = 0.5;
	}
    }
  numHM = new float [m_numberFrames * nbMus];
  denHM = new float [m_numberFrames * nbMus];
  numWM = new float [m_numberBins * nbMus];
  denWM = new float [m_numberBins * nbMus];
  
  sXHat = new float [m_numberFrames * m_numberBins];
  
  for (niter=0; niter<m_maxIterations; niter++)
    {
      // computing the approximated frame:
      /* sF0Hat = m_HF0 m_WF0 (matrix product) */
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberFrames, m_numberBins, m_numberOfF0, 1.0,
		  m_HF0, m_numberOfF0, m_WF0, m_numberBins, 0.0, 
		  sF0Hat, m_numberBins);
      cout << "ok1"<<endl;
      /* sPhiHat = m_HGAMMA m_WGAMMA*/
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberFrames, m_numberBins, m_numberOfGamma, 1.0,
		  m_HGAMMA, m_numberOfGamma, m_WGAMMA, m_numberBins, 0.0,
		  sPhiHat, m_numberBins);
      cout << "ok2"<<endl;
      /* sF0Hat = m_HF0 m_WF0 (matrix product) */
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberFrames, m_numberBins, nbMus, 1.0,
		  HM, nbMus, WM, m_numberBins, 0.0, 
		  sM, m_numberBins);
      cout << "ok3"<<endl;
      // update HF0:
      for (n=0; n<m_numberFrames; n++)
	{
	  for(f=0; f<m_numberBins; f++)
	    {
	      // note: tmpM matrices are F x N 
	      // and not N x F (like m_SX):
	      // this helps for sgemm. A trick to go faster?
	      
	      sXHat[n*m_numberBins+f] = (sPhiHat[n*m_numberBins+f] * 
					 sF0Hat[n*m_numberBins+f] + 
					 sM[n*m_numberBins+f]);
	      
	      /* 
	                            m_SX
	         tmpNum^T = --------------------  
		             sPhiHat * sF0Hat^2

	                            m_SX * sPhiHat
	         tmpNum^T = --------------------  
		             (sPhiHat * sF0Hat + sM )**2
	      */
	      tmpNum[f*m_numberFrames+n] = 
		m_SX[n*m_numberBins+f] * sPhiHat[n*m_numberBins+f] / 
		max(sXHat[n*m_numberBins+f] * sXHat[n*m_numberBins+f], 
		    eps);
	      
	      /* 
	                       1
	         tmpDen^T = --------
		             sF0Hat 

	                       sPhiHat
	         tmpDen^T = --------
		             sF0Hat *sPhiHat + sM
	      */
	      tmpDen[f*m_numberFrames+n] = 
		 sPhiHat[n*m_numberBins+f] / max(sXHat[n*m_numberBins+f], eps);
	    }
	}
      /* numeraHF0 = m_WF0 tmpNum */
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  m_numberOfF0, m_numberFrames, m_numberBins, 1.0,
		  m_WF0, m_numberBins, tmpNum, m_numberFrames, 0.0,
		  numeraHF0, m_numberFrames);
      cout << "ok4"<<endl;
      /* denomiHF0 = m_WF0 tmpDen */
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberOfF0, m_numberFrames, m_numberBins, 1.0,
		  m_WF0, m_numberBins, tmpDen, m_numberFrames, 0.0,
		  denomiHF0, m_numberFrames);
      cout << "ok5"<<endl;
      for (n=0; n<m_numberFrames; n++)
	{
	  for(f=0; f<m_numberOfF0; f++)
	    {
	      /* m_HF0 <- m_HF0 * numeraHF0^T / denomiHF0^T  */
	      m_HF0[n*m_numberOfF0+f] = m_HF0[n*m_numberOfF0+f] *
		(numeraHF0[f*m_numberFrames+n]) / 
		max(denomiHF0[f*m_numberFrames+n], eps);
	    }
	}
      // computing the approximated frame:
      /* sF0Hat = m_HF0 m_WF0 (matrix product) */
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberFrames, m_numberBins, m_numberOfF0, 1.0,
		  m_HF0, m_numberOfF0, m_WF0, m_numberBins, 0.0, 
		  sF0Hat, m_numberBins);
      cout << "ok6"<<endl;
      
      // update HGAMMA:
      /* equations analogous to above ones */
      for (n=0; n<m_numberFrames; n++)
	{
	  for(f=0; f<m_numberBins; f++)
	    {
	      	      
	      sXHat[n*m_numberBins+f] = (sPhiHat[n*m_numberBins+f] * 
					 sF0Hat[n*m_numberBins+f] + 
					 sM[n*m_numberBins+f]);
	      
	      // note: tmpM matrices are F x N 
	      // and not N x F (like m_SX):
	      // this helps for sgemm.
	      tmpNum[f*m_numberFrames+n] = 
		m_SX[n*m_numberBins+f] * sF0Hat[n*m_numberBins+f] / 
		max(sXHat[n*m_numberBins+f]*sXHat[n*m_numberBins+f], 
		    eps);
	      
	      tmpDen[f*m_numberFrames+n] = 
		sF0Hat[n*m_numberBins+f] / 
		max(sXHat[n*m_numberBins+f], eps);
	    }
	}
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  m_numberOfGamma, m_numberFrames, m_numberBins, 1.0,
		  m_WGAMMA, m_numberBins, tmpNum, m_numberFrames, 0.0,
		  numeraHG, m_numberFrames);
      cout << "ok8"<<endl;
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberOfGamma, m_numberFrames, m_numberBins, 1.0,
		  m_WGAMMA, m_numberBins, tmpDen, m_numberFrames, 0.0,
		  denomiHG, m_numberFrames);
      cout << "ok8"<<endl;
      for (n=0; n<m_numberFrames; n++)
	{
	  maxGamma[n] = 0;
	  for(f=0; f<m_numberOfGamma; f++)
	    {
	      m_HGAMMA[n*m_numberOfGamma+f] = m_HGAMMA[n*m_numberOfGamma+f] *
		(numeraHG[f*m_numberFrames+n]) / 
		max(denomiHG[f*m_numberFrames+n], eps);
	      if (maxGamma[n]<m_HGAMMA[n*m_numberOfGamma+f])
		{
		  maxGamma[n] = m_HGAMMA[n*m_numberOfGamma+f];
		}
	    }
	}

      /* Normalization */
      for (n=0; n<m_numberFrames; n++)
	{
	  if(maxGamma[n]==0)
	    { //reset the values of HGAMMA to uniform distrib, HF0 to eps
	      for (f=0; f<m_numberOfGamma; f++)
		{
		  m_HGAMMA[n*m_numberOfGamma+f] = 1/m_numberOfGamma;
		}
	      for (f=0; f<m_numberOfF0; f++)
		{
		  m_HF0[n*m_numberOfF0+f] = eps;
		}
	    }
	  else
	    {
	      for (f=0; f<m_numberOfGamma; f++)
		{
		  m_HGAMMA[n*m_numberOfGamma+f] = m_HGAMMA[n*m_numberOfGamma+f] / 
		    maxGamma[n];
		}
	      for (f=0; f<m_numberOfF0; f++)
		{
		  m_HF0[n*m_numberOfF0+f] = m_HF0[n*m_numberOfF0+f] * maxGamma[n];
		}
	    }
	}

      /* sPhiHat = m_HGAMMA m_WGAMMA*/
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberFrames, m_numberBins, m_numberOfGamma, 1.0,
		  m_HGAMMA, m_numberOfGamma, m_WGAMMA, m_numberBins, 0.0,
		  sPhiHat, m_numberBins);
      cout << "ok9"<<endl;

      // Update accompaniment elements
      for (n=0; n<m_numberFrames; n++)
	{
	  for(f=0; f<m_numberBins; f++)
	    {
	      
	      sXHat[n*m_numberBins+f] = (sPhiHat[n*m_numberBins+f] * 
					 sF0Hat[n*m_numberBins+f] + 
					 sM[n*m_numberBins+f]);
	      
	      // note: tmpM matrices are F x N 
	      // and not N x F (like m_SX):
	      // this helps for sgemm.
	      tmpNum[f*m_numberFrames+n] = 
		m_SX[n*m_numberBins+f] / 
		max(sXHat[n*m_numberBins+f]*sXHat[n*m_numberBins+f], 
		    eps);
	      
	      tmpDen[f*m_numberFrames+n] = 
		1 / 
		max(sXHat[n*m_numberBins+f], eps);
	    }
	}
      // computing the num and denom for the update
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  nbMus, m_numberFrames, m_numberBins, 1.0,
		  WM, m_numberBins, tmpNum, m_numberFrames, 0.0,
		  numHM, m_numberFrames);
      cout << "ok10"<<endl;
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  nbMus, m_numberFrames, m_numberBins, 1.0,
		  WM, m_numberBins, tmpDen, m_numberFrames, 0.0,
		  denHM, m_numberFrames);
      cout << "ok11"<<endl;
      for (n=0; n<m_numberFrames; n++)
	{
	  for(f=0; f<nbMus; f++)
	    {
	      HM[n*nbMus+f] = HM[n*nbMus+f] *
		(numHM[f*m_numberFrames+n]) / 
		max(denHM[f*m_numberFrames+n], eps);
	    }
	}
      // computing the num and denom for the update
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		  m_numberBins, nbMus, m_numberFrames, 1.0,
		  tmpNum, m_numberFrames, HM, nbMus, 0.0,
		  numWM, nbMus);
      cout << "ok12"<<endl;
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberBins, nbMus, m_numberFrames, 1.0,
		  tmpDen, m_numberFrames, HM, nbMus, 0.0,
		  denWM, nbMus);
      cout << "ok13"<<endl;
      for (n=0; n<m_numberBins; n++)
	{
	  for(f=0; f<nbMus; f++)
	    {
	      WM[f*m_numberBins+n] = WM[f*m_numberBins+n] *
		(numWM[n*nbMus+f]) / 
		max(denWM[n*nbMus+f], eps);
	    }
	}
      
    }

  // freeing memory
  delete[] denomiHF0;
  delete[] numeraHF0;
  delete[] tmpDen ;
  delete[] tmpNum;
  delete[] numeraHG;
  delete[] denomiHG;
  delete[] sF0Hat;
  delete[] sPhiHat;
  delete[] maxGamma;  
  delete[] numHM;
  delete[] denHM;
  delete[] numWM;
  delete[] denWM;
  delete[] WM;
  delete[] HM;
  delete[] sM;
}

F0Salience::FeatureSet
F0Salience::process(const float *const *inputBuffers, 
		    Vamp::RealTime timestamp)
{
  size_t f, nn, nch, nf;
  // FeatureSet fs;
  float sumOfFramePower = 0.0, frameIsNotNull = 0.0;
  
  if(m_blockSize == 0) {
    cerr << "ERROR: F0Salience::process: Not initialised" << endl;
    return FeatureSet();
  }
  
  // size_t n = m_numberBins; // for convenience? but not very readable...
  size_t N_notes = m_numberOfF0;

  if (m_nbAccumulatedFrames==0)
    {
      // initialize the power spectrum matrix
      m_SX = new float [m_numberFrames * m_numberBins];
    }
  
  /*
    although above has m_numberFrames, below works only
    for  m_numberFrames=1
    Not sure if it would be useful to work on matrices.

    EDIT 20140624: extending to m_numberFrames>1
  */
  for (f=0; f<m_numberBins; f++)
    {
      m_SX[m_nbAccumulatedFrames * m_numberBins + f] = 
	(inputBuffers[0][2*f]*inputBuffers[0][2*f]) + 
	(inputBuffers[0][2*f+1]*inputBuffers[0][2*f+1]); 
      // taking first channel power spectrum. 
      for(nch=1;nch<m_nbChannels;nch++)
	{
	  m_SX[m_nbAccumulatedFrames * m_numberBins + f] += 
	    (inputBuffers[nch][2*f]*inputBuffers[nch][2*f]) + 
	    (inputBuffers[nch][2*f+1]*inputBuffers[nch][2*f+1]); 
	}
      // the power summation does not make sense anymore... 
      sumOfFramePower += m_SX[m_nbAccumulatedFrames * m_numberBins + f];
    }
  m_nbAccumulatedFrames += 1; // keep track of how many frames we stored
  
  if (m_computeBasis==0)
    {
      m_WF0 = new float [N_notes * m_numberBins];
      m_f0Table = new float [N_notes];
      m_WGAMMA = new float [m_numberOfGamma * m_numberBins];
      m_HF0 = new float [m_numberFrames * N_notes];
      m_HGAMMA = new float [m_numberFrames * m_numberOfGamma];
      
      cout << "matrices initialized again" << endl;
      
      cout << "generate WGAMMA" << endl;
      genHannBasis();
      cout << "generate WF0" << endl;
      m_numberOfF0 = genBaseWF0(m_fmin, m_fmax, m_inputSampleRate, 
				m_blockSize, m_stepnote, m_opencoef, 
				m_WF0, m_f0Table);
      N_notes = m_numberOfF0;
      cout << "matrices computed" << endl;
      m_computeBasis = 1; // to avoid computing this all the time...
    }
  
  
  if(m_nbAccumulatedFrames == m_numberFrames)
    { 
      // Compute the decomposition when the matrix is full:
      // initializing:
      if (sumOfFramePower > 0.0)
	{
	  frameIsNotNull = 1.0;
	}
      cout << "initializing H matrices:" << endl;
      for(nf=0; nf<m_numberFrames; nf++)
	{
	  for(nn=0; nn<m_numberOfGamma; nn++)
	    {
	      m_HGAMMA[nf*m_numberOfGamma + nn] = 1.0 * frameIsNotNull; // 0 if frame is null
	    }
	  for(nn=0; nn<m_numberOfF0; nn++)
	    {
	      m_HF0[nf*m_numberOfF0 + nn] = 1.0 * frameIsNotNull;
	    }
	}
      
      // make the calculations:
      // cout << "updating H matrices:" << endl;
      if (sumOfFramePower > 0.0) // no need to compute for null frames
	{
	  updateHmatrices();
	}
  
      // fill in the features:
      float *sPhiHat;
      sPhiHat = new float [m_numberFrames * m_numberBins];
      /* sPhiHat = m_HGAMMA m_WGAMMA*/
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  m_numberFrames, m_numberBins, m_numberOfGamma, 1.0,
		  m_HGAMMA, m_numberOfGamma, m_WGAMMA, m_numberBins, 0.0,
		  sPhiHat, m_numberBins);
  
      for (nf=0; nf<m_numberFrames; nf++)
	{
	  Feature f0salienceFeat;
	  f0salienceFeat.hasTimestamp = false;
	  f0salienceFeat.timestamp = timestamp 
	    - Vamp::RealTime::frame2RealTime((m_numberFrames-1-nf)*m_stepSize,
					     lrintf(m_inputSampleRate));
	  f0salienceFeat.values.reserve(m_numberOfF0);
	  
	  float maxInHF0 = *max_element(m_HF0,m_HF0+m_numberOfF0);
	  float targetMinInHF0 = 0.; // DEBUG // invdb((db(maxInHF0)-m_thresholdEnergy));
	  
	  for(nn=0; nn<m_numberOfF0; nn++) 
	    // we do not return the last bin: noise energy 
	    // (TODO: include it as another output descriptor)
	    {
	      f0salienceFeat.values.push_back(max(
                  m_HF0[nf*m_numberOfF0 + nn], 
		  targetMinInHF0));
	    }
	  
	  m_featureSet[0].push_back(f0salienceFeat);
	  
	  Feature envelopeFeat;
	  envelopeFeat.hasTimestamp = false;
	  envelopeFeat.values.reserve(m_numberBins);
	  for(f=0; f<m_numberBins; f++)
	    {
	      envelopeFeat.values.push_back(
                  sPhiHat[nf*m_numberBins + f]);
	    }
	  m_featureSet[1].push_back(envelopeFeat);
	}
      
      // should we free some of these arrays?
      // delete [] m_WF0;
      // delete [] m_HF0;
      // delete [] m_WGAMMA;
      // delete [] m_HGAMMA;
      // delete [] m_f0Table;
  
      // delete [] m_SX;
      m_nbAccumulatedFrames = 0; // resetting buffer pointer to beginning of matrix
      // should we also reset the content of matrix? 
    } // end of if m_nbAccumulatedFrames == m_numberFrames

  return FeatureSet(); // return the featureSet, in any case (even if empty, no prob with this?)
}

F0Salience::FeatureSet
F0Salience::getRemainingFeatures()
{
  return m_featureSet;
}

