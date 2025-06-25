#ifndef __RAT_VertexGen_PhotonRay__
#define __RAT_VertexGen_PhotonRay__

#include <RAT/GLG4VertexGen.hh>
#include <CLHEP/Random/RandGeneral.h>

namespace RAT {

  class VertexGen_PhotonRay : public GLG4VertexGen {
  public:
    VertexGen_PhotonRay(const char *arg_dbname="pRay");
    virtual ~VertexGen_PhotonRay();
    virtual void GeneratePrimaryVertex( G4Event *argEvent,
					G4ThreeVector &dx,
					G4double dt);
    /** State format "num_photons wavelength_nm" */
    virtual void SetState( G4String newValues );
    virtual G4String GetState();
    virtual G4ThreeVector RodriguesRotationFormula(G4ThreeVector point, G4ThreeVector axis, double rotationAngle);

    // Diffuser data from file functions
    virtual G4ThreeVector SampleFromDataDistribution( std::vector<std::vector<double>> diffuserData );
    virtual std::vector<std::vector<double>> ReadDiffuserData(std::string filePath);
    virtual double InterpolateDiffuserData(std::vector<std::vector<double>> diffuserData, double interpolationAngleTheta, double interpolationAnglePhi);
    virtual G4ThreeVector TheoreticalMomentumDirectionGeneration();

  private:
    G4ParticleDefinition *fOpticalPhoton;
    int fNumPhotons;
    double fEnergy;
    CLHEP::RandGeneral* fRndmEnergy;
    double fMinEnergy;
    double fMaxEnergy;
    std::string fMaterial;
    double fExpTime;
	double macroTheta;
	double macroPhi;
  };


} // namespace RAT

#endif
