#include <RAT/VertexGen_PhotonRayWithAngularDist.hh>
#include <RAT/LinearInterp.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <Randomize.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4ThreeVector.hh>
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <TRandom.h>
#include <random>
using namespace std;





namespace RAT {
  VertexGen_PhotonRay::VertexGen_PhotonRay(const char *arg_dbname) // VertexGen_PhotonRay constructor
    : GLG4VertexGen(arg_dbname)
  {
    // std::cout << "Constructing PhotonRay class:\n";
    fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
    SetState("1 400 0.0 0.0"); // one photon per event, 400 nm wavelength, angles both 0.0
    fRndmEnergy = 0;
    fMinEnergy = 0.0;
    fMaxEnergy = 0.0;
    fMaterial = "";
    macroTheta = 0.0;
    macroPhi = 0.0;
  }


  VertexGen_PhotonRay::~VertexGen_PhotonRay() // Destructor
  {
    delete fRndmEnergy;
  }



  G4ThreeVector VertexGen_PhotonRay::RodriguesRotationFormula(G4ThreeVector point, G4ThreeVector axis, double rotationAngle)
  {
    // Rotates in 3D, a point by an angle theta around the given axis
    // Point is a point in Cartesian coordinates (important)
    // Axis should be a unit vector (assumed to be normalised here, could be changed to make axis a unit vector by itself) (also in Cartesian coordinates)
    // Angle is in radians
    return point * cos(rotationAngle) + ( axis.cross(point) ) * sin(rotationAngle) + axis * (axis.dot(point)) * (1-cos(rotationAngle));
  }





  double VertexGen_PhotonRay::InterpolateDiffuserData(std::vector<std::vector<double>> diffuserData, double interpolationAngleTheta, double interpolationAnglePhi)
  {
    //std::cout << "Starting InterpolateDiffuserData\n";
    std::vector<double> diffuserDataThetas; // column index 4 // Now column index 3 because the first column has been removed
    std::vector<double> diffuserDataPhis; // column index 6 // Now column index 5 because the first column has been removed
    std::vector<double> powerProbabilities; // column index 7 // Now column index 6 because the first column has been removed

    // Get all the theta and phi angles, as well as the "probabilities" (powers) for each direction
    for (unsigned long int i=0; i < diffuserData.size(); i++) {
      diffuserDataThetas.push_back(diffuserData[i][3] - 32.68); // -32.68 to correct the angles from the data back to a 0->theta range rather than 32.68->maxAngle range
      diffuserDataPhis.push_back(diffuserData[i][5]);
      powerProbabilities.push_back(diffuserData[i][6]);
    }
    //std::cout << "Got the theta and phi angles\n";

    // For interpolation, need 1) The point to interpolate to in angular coordinates  2) the closest stored point to this
    std::vector<double> angularDistances;
    int minDistanceIndex;
    double minDistanceAngular;
    double distanceAngular;
    //std::cout << "Starting distance checking\n";
    for (unsigned long int i=0; i < diffuserData.size(); i++) {
      distanceAngular = pow(pow(diffuserDataThetas[i] - interpolationAngleTheta, 2) + pow(diffuserDataPhis[i] - interpolationAnglePhi, 2), 0.5);
      angularDistances.push_back( distanceAngular  );
      if (i==0) {
        minDistanceIndex = 0;
        minDistanceAngular = distanceAngular;
      }
      else {
        if (distanceAngular < minDistanceAngular) {
          minDistanceIndex = i;
          minDistanceAngular = distanceAngular;
        }
      }
    }

    // Currently found the closest data point to the requested point
    // For now just use the nearest values' z value as the z value to return
    // Kind of like a flat interpolation. 
    // Currently no handling of really far away zs and this still means a harder cutoff at some point where the boundary box is made, but that can be improved later with a better interpolation
    // The water data goes from -90 to 90 degrees, and with the diffuser construction being the way it is, any photons emitted backwards are going to be extremely rare and unlikely to be relevant.
    //std::cout << "Found min distance index, angle " << diffuserDataThetas[minDistanceIndex] << " " << diffuserDataPhis[minDistanceIndex] << "\n";
    double interpolatedProbValue = powerProbabilities[minDistanceIndex];
    return interpolatedProbValue;

  }


  G4ThreeVector VertexGen_PhotonRay::SampleFromDataDistribution(std::vector<std::vector<double>> diffuserData)
  {

    // For the interpolation, want to:
    // 1) Load the data (already taken care of outside the function)
    // 2) Set the boundary box for the random number sampling (need maximum power value (min is 0), min+max phi and theta)
    // 3) Pick a random point in the (3D) box
    // 4) Find the interpolated value of the "probability" (power) distribution
    // 5) Check if the random z value is lower than the interpolated z value
    // 6) If it is, keep the point, otherwise try again with a new random (3D) point.

    double boundaryMinTheta;
    double boundaryMinPhi;
    double boundaryMaxTheta;
    double boundaryMaxPhi;
    double boundaryMinProb = 0;
    double boundaryMaxProb;

    std::vector<double> diffuserDataThetas; // column index 4 // Now column index 3 because the first column has been removed
    std::vector<double> diffuserDataPhis; // column index 6 // Now column index 5 because the first column has been removed
    std::vector<double> powerProbabilities; // column index 7 // Now column index 6 because the first column has been removed

    //std::cout << "Putting diffuserData into vectors\n";
    // Get all the theta and phi angles, as well as the "probabilities" for each direction
    for (unsigned long int i=0; i < diffuserData.size(); i++) {
      diffuserDataThetas.push_back(diffuserData[i][3]-32.68); // The data has a ~32.68 degree offset in theta, reducing it makes it more straightforward to determine the correct direction of the beam
      diffuserDataPhis.push_back(diffuserData[i][5]); 
      powerProbabilities.push_back(diffuserData[i][6]); 
    }
    //std::cout << "Finished putting diffuserData into vectors\n";

    boundaryMinTheta = *std::min_element( diffuserDataThetas.begin(), diffuserDataThetas.end() );
    boundaryMaxTheta = *std::max_element( diffuserDataThetas.begin(), diffuserDataThetas.end() );
    boundaryMinPhi = *std::min_element( diffuserDataPhis.begin(), diffuserDataPhis.end() );
    boundaryMaxPhi = *std::max_element( diffuserDataPhis.begin(), diffuserDataPhis.end() );
    boundaryMaxProb = *std::max_element( powerProbabilities.begin(), powerProbabilities.end() );

    // now we have the boundaries for the 3D box, pick a random uniformly distributed point within it
    // basic random number generation from cppreference
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> thetaDist(boundaryMinTheta, boundaryMaxTheta);
    std::uniform_real_distribution<> phiDist(boundaryMinPhi, boundaryMaxPhi);
    std::uniform_real_distribution<> probDist(boundaryMinProb, boundaryMaxProb);

    bool findingValidDirection(true); // Start looping until a valid direction is found
    G4ThreeVector randomPoint(0,0,0);

    //std::cout << "Starting interpolation:\n";
    while (findingValidDirection == true) {
      randomPoint.setX(thetaDist(gen));
      randomPoint.setY(phiDist(gen));
      randomPoint.setZ(probDist(gen));
      
      //std::cout << "Test value: " << randomPoint.getX() << " " << randomPoint.getY() << " " << randomPoint.getZ() << "\n";
      //std::cout << "Interpolated value: " << InterpolateDiffuserData(diffuserData, randomPoint.getX(), randomPoint.getY()) << "\n";

      if ( randomPoint.getZ() <= InterpolateDiffuserData(diffuserData, randomPoint.getX(), randomPoint.getY()) ) {
        // Point has been found if it's <= the interpolated probability distribution, this ensures random numbers produced with the correct probability distribution
        findingValidDirection = false;
      }

    }
    //std::cout << "Finished interpolation:\n";
    //std::cout << "Resulting angles: " << randomPoint.getX() << " " << randomPoint.getY() << "\n";

    // Now have a valid theta and phi, need to convert these into a momentum direction in cartesian coordinates

    // Need to do Rodrigues rotation due to the way the coordinate system is laid out
    G4ThreeVector initialDirectionPoint(-1,0,0); // start at unit radius on the -x axis
    G4ThreeVector thetaRotationAxis(0,1,0);
    G4ThreeVector phiRotationAxis(0,0,1);
    G4ThreeVector directionAfterThetaRotation( RodriguesRotationFormula(initialDirectionPoint, thetaRotationAxis, randomPoint.getX()/180*CLHEP::pi) );
    G4ThreeVector directionAfterPhiRotation( RodriguesRotationFormula(directionAfterThetaRotation, phiRotationAxis, randomPoint.getY()/180*CLHEP::pi) );
    //std::cout << "randomPoint.getY(): "  << randomPoint.getY() << "\n";

    //std::cout << "Cartesian Point: " << directionAfterPhiRotation.getX() << " " << directionAfterPhiRotation.getY() << " " << directionAfterPhiRotation.getZ() << "\n";

    return directionAfterPhiRotation;
  
  
  }




  std::vector<std::vector<double>> VertexGen_PhotonRay::ReadDiffuserData(std::string filePath)
  {
    //std::cout << "ReadDiffuserData started\n";
    string DataLine;

    std::vector<std::vector<double>> diffuserDataPoints;
    ifstream DiffuserData(filePath);
    if (DiffuserData.is_open()) {
      std::cout << "Diffuser data opened\n";
    }
    else {
      std::cout << "Diffuser data not opened correctly.\n";
    }
 
    int lineNumber=0;
    std::cout << "Starting diffuserData loading loop:\n";
    while (std::getline (DiffuserData, DataLine)) {
      if (lineNumber==0) {
        lineNumber+=1;
        continue;
      } // Skip the first line because it's a header line
      
      std::stringstream ss(DataLine);
      std::string token;
      double tokenDouble;
      std::vector<double> dataLineToPush;
      //std::cout << "Example data pushing " << DataLine << "\n";
      int tokenCount(0); // Try skipping the first column of data since it's a date that isn't properly converted to a double
      while(getline(ss, token, ','))
      {
        if (tokenCount==0) {
            tokenCount+=1;
            continue;
        }
        tokenDouble = stod(token); 
        
        dataLineToPush.push_back(tokenDouble);
        //std::cout << "Token being pushed: " << tokenDouble << "\n";

        tokenCount += 1;
      }
      diffuserDataPoints.push_back(dataLineToPush);
      //std::cout << "Line being pushed: \n";
      //if (lineNumber != 0) {
        //std::cout << diffuserDataPoints[0][0] << "\n";
      //}
      lineNumber+=1;
      
    }

    // Close the file
    DiffuserData.close();
    std::cout << "Finishing diffuser data load.\n";

    return diffuserDataPoints;
  }





G4ThreeVector VertexGen_PhotonRay::TheoreticalMomentumDirectionGeneration() 
{
  double a_max= 0.698132; // 0.698132rad ~40deg 

  double temp_theta = 0.25; //0.698; / the angle to get 4 pmts 0.25;
  double theta_init=temp_theta; //0.25; //16.47;
  double cosTheta=0.;
  double maxCosTheta = cos(a_max);


  G4ThreeVector momentumDirection;
  /*
  // Alternative method for producing random values for cos theta
  // Want to randomly create number between cos(colldivrad) and 1
  // double colldivrad = 0.698132; // The (half)angle width of the beam, 40 degrees.
  // double range = 1-cos(colldivrad);
  // double theta_init =  acos((range * G4UniformRand()) + (1-range));
  
  
  // double theta_init = (f1->GetRandom())*CLHEP::pi/180;
  // double phi_init = ((r1->Uniform(16))*CLHEP::pi/180);
  */


  // This method produces a random value using std::rand (the maximum value of which varies based on implementation)
  // Limits it to 0->9999, then add 1 and divide by 10000 to get a random number 0->1
  // (with 10,000 steps, is this good enough? Especially given the non-linear nature of the cos function)
  // Rescales 0->1 to the range cos(amax)->1 - This changes the "fine-grained-ness" of the points; i.e. larger cones will have more coarse positioning than smaller cones
  // Then makes a theta from this
  cosTheta = ((std::rand()%10000+1.0)/10000.0 * (1-maxCosTheta) + maxCosTheta);
  theta_init = acos(cosTheta);
  //std::cout << theta_init << "\n";

  // Alternative method with 0->1 cos ranges, maintains positional resolution for different cone widths which may be desired, but needs multiple runs to produce a valid angle.
  /*
  while(not_open==true){
    maxCosTheta = cos(a_max);
    cosTheta = ((std::rand()%10000)+1.0)/10000.0; // //r1->Uniform();
    theta_init = acos(cosTheta); ///2/1.125+100);
    if(theta_init<a_max && theta_init>a_min ){
        not_open=false;}
  }
  not_open = true; // Resets not_open for the next photon
  */



  // Generate a random value for phi in the range 0->2pi.
  double phi_init = 2.0 * G4UniformRand() * CLHEP::pi;
  // double phi_init = (f1->GetRandom())*CLHEP::pi/180;
  
  // theta_init and phi_init describe 2 angles on the cone: theta_init = angle out from centre line of cone, phi_init = angle around cone
  // so if phi=0 is one side of the cone, phi=pi is the other side      


  // For the diffuser's directionality, theta_init and phi_init need to be transformed according to the generator's direction.

  // Need to convert coordinate system from the local theta_init and phi_init given relative to the beam to the global theta and phi actual locations
  // Starting position - this is some point in the beam, relative to the beam, with the beam pointing towards the global theta=0, phi=0 (or z direction in Cartesian)
  double sphericalR = 1;
  double sphericalTheta = theta_init;
  double sphericalPhi = phi_init;
  G4ThreeVector sphericalCoord(sphericalR, sphericalTheta, sphericalPhi);

  //std::cout << "Spherical coord before transform: " << sphericalR << ", " << sphericalTheta << ", " << sphericalPhi << "\n";      

  // Convert this position to cartesian coordinates to use G4ThreeVector's dot and cross product methods.
  // These conversions effectively set the arrangement of the cartesian axes with respect to theta_init and phi_init, as well as with respect to the global coordinates.
  double cartesianX = sphericalR * sin(sphericalTheta) * cos(sphericalPhi);
  double cartesianY = sphericalR * sin(sphericalTheta) * sin(sphericalPhi);
  double cartesianZ = sphericalR * cos(sphericalTheta);
  G4ThreeVector cartesianCoord(cartesianX, cartesianY, cartesianZ);


  return cartesianCoord;
}










  void VertexGen_PhotonRay::GeneratePrimaryVertex(G4Event *event,
                                                  G4ThreeVector &dx,
                                                  G4double dt)
  {
    //double a_min = 0.0;
   
    
    ofstream varfile;
    varfile.open("variable_check.txt"); 
    varfile << "Theta_Init Phi_Init xDir yDir zDir\n"; 

    bool useDiffuserData(true);


    /*
    // Parameters for the (pol7) random number generation, determined from diffuser data(?)
    vector<double> par = {-299.339, 183467, -67597.5, 118294, -147280, 75690, -18345.4, 1736.95};
    //vector<double> par = {3.553944665124477,-0.000513139723363017, -0.000503803173391664, 1.907455214523768e-06, 4.188507416994249e-07, -1.389170400703575e-09, -2.973373805790281e-10};
    /// load in the diffuser profile parameters
    */
    /*
    std::string line;
    ifstream parfile;
    parfile.open("/user/atarrant/hepstore/data/pray/para_dif.txt");

    if(parfile.is_open()){
      while(std::getline(parfile, line)){
	      par.push_back(std::stod(line));
      }
    }
    else{
      cout << "no dif file" << endl;  
    }
    */

    // Alternative setup for generating random numbers on a polynomially-described custom probability distribution
    //    TF1 *f1 = new TF1("f1","pol7(0)",160, 170);
    //TF1 *f1 = new TF1("f1","pol7(0)",0,80); // Make a random distribution called f1, using a 7-degree polynomial with parameters numbered starting at 0 (giving 0-7), with a minimum value of 0 and a max value of 80
    //f1->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7]); // Parameters for the random number generation
    //f1->Draw();
    //TRandom *r1 = new TRandom();


    std::vector<std::vector<double>> diffuserData;
    
    if (useDiffuserData==true) {
      //bool not_open=true;
      std::cout << "Starting diffuser data loading\n";
	  
      std::string filePath("DiffuserWaterDataClean.out");
      //std::vector<std::vector<double>> diffuserData(ReadDiffuserData(filePath));
      diffuserData = ReadDiffuserData(filePath);
      std::cout << "Example diffuser data: " << diffuserData[1][4] << "\n";
    }

    
    //find the photons for the pray
    for (int i = 0; i < fNumPhotons; i++) {
      
	    //theta_init=temp_theta;
      //not_open=true;
      // Choose angles of injection and divergence of beam. Position of beam injection point currently set in macro
      G4ThreeVector mom;
      G4ThreeVector momentumDirection;
      

      if (useDiffuserData == true) {
        // Second is the method based on the diffuser data:

        //std::cout << "Starting sampleFromDataDistribution\n";
        momentumDirection = (SampleFromDataDistribution(diffuserData)); // make a random momentum direction based on random sampling of interpolated probability data 
        //std::cout << "Finished sampleFromDataDistribution\n";
      }

      
      else {
        // First are the old methods with a min/max cos theta and random uniform on cos(theta)
        momentumDirection = (TheoreticalMomentumDirectionGeneration()); // Want momentumDirection to be a unit vector in 3D
      }






      // These momentum directions then need to be rotated by the macro theta/phi for diffuser directionality.

      //std::cout << "Cartesian coord before transform: " << cartesianX << " " << cartesianY << " " << cartesianZ << "\n";

      // Next need to do some rotation transformations with Rodrigues' formula - a general formula for rotations in 3D about an axis of rotation
      // Take as the coordinate system that (theta = pi/2, phi=0) points in x direction, (theta=0) points in z and (theta=pi/2, phi=pi/2) points in y direction.

      // Fist rotate around the y-axis by the beam's given theta direction
      
      G4ThreeVector axis(0,1,0);
      
      G4ThreeVector firstRotationCartesianCoord = ( RodriguesRotationFormula(momentumDirection, axis, macroTheta) );

      // Second is to rotate forwards round the z axis to the beam's phi
      G4ThreeVector secondAxis(0,0,1);
      G4ThreeVector finalRotationCartesianCoord = ( RodriguesRotationFormula(firstRotationCartesianCoord, secondAxis, macroPhi) );

      //std::cout << "Cartesian coord after 2 transform: " << finalRotationCartesianCoord.getX() << " " << finalRotationCartesianCoord.getY() << " " << finalRotationCartesianCoord.getZ() << "\n";

      // This should give a cartesian coordinate for the actual location of the point on the sphere in the global coordinate system.
      
      // At this point, it's easier to use the Cartesian direction directly rather than attempt to convert back to spherical coordinates and have to deal with quadrants, etc.
      



      // Use fixed energy unless spectrum was provided
      double energy;
      if(fRndmEnergy)
	      energy = fMinEnergy + (fMaxEnergy - fMinEnergy) * fRndmEnergy->shoot();
      else
	      energy = fEnergy;
      

      //mom.setRThetaPhi(energy, theta, phi); // Momentum == energy units in GEANT4


      // Set the momentum in Cartesian coordinates rather than spherical to avoid potential quadrant issues when using asin, acos, atan, etc. 
      //mom.setRThetaPhi(energy, phi, theta);
      // Convert direction to momentum; direction has a magnitude of 1, so just need to multiply by energy to get the magnitude of the vector to be energy with the same direction. 
      double finalX(finalRotationCartesianCoord.getX());
      double finalY(finalRotationCartesianCoord.getY());
      double finalZ(finalRotationCartesianCoord.getZ());
      double momentumX = finalX * energy;
      double momentumY = finalY * energy;
      double momentumZ = finalZ * energy;
      mom.set(momentumX, momentumY, momentumZ);
      


      // Distribute times expoenentially, but don't bother picking a
      // random number if there is no time constant
      double expt = 0.0;
      if(fExpTime > 0.0)
	      expt = - fExpTime * log(G4UniformRand());
      G4PrimaryVertex* vertex= new G4PrimaryVertex(dx, dt + expt);

      // Make a photon with the specified momentum
      G4PrimaryParticle* photon = new G4PrimaryParticle(fOpticalPhoton, momentumX, momentumY, momentumZ);

      // Generate random polarization
      double polarisationPhi = (G4UniformRand()*2.0-1.0)*CLHEP::pi;
      G4ThreeVector e1 = mom.orthogonal().unit();
      G4ThreeVector e2 = mom.unit().cross(e1);
      G4ThreeVector pol = e1*cos(polarisationPhi)+e2*sin(polarisationPhi);
      photon->SetPolarization(pol.x(), pol.y(), pol.z());
      photon->SetMass(0.0); // Seems odd, but used in GLG4VertexGen_Gun


      
      // varfile << theta_init << " " << phi_init << " " << " " << finalX << " " << finalY << " " << finalZ << "\n" ;
      //varfile << energy << " " << pol.x() << " " <<  pol.y() << " " << pol.z() << " " << dx << " " << dt << " " << dt +expt << "\n";

      vertex->SetPrimary(photon);
      event->AddPrimaryVertex(vertex);
    } // loop over photons
    varfile.close();
  }

  void VertexGen_PhotonRay::SetState(G4String newValues)
  {
    if (newValues.length() == 0) {
      // print help and current state
      G4cout << "Current state of this VertexGen_PhotonBeam:\n"
             << " \"" << GetState() << "\"\n" << G4endl;
      G4cout << "Format of argument to VertexGen_PhotonBeam::SetState: \n"
        " \"num_photons wavelength_nm angleTheta anglePhi\"\n" << G4endl;
      return;
    }

    std::istringstream is(newValues.c_str());
    int num, wavelength;
    // Going to edit this to add angles, need two to specify rotation
    std::cout << "Check istringstream:\n";
    std::cout << is.str() << "\n";
    double macroArgTheta(0);
    double macroArgPhi(0); // Assume degrees passed in (changed to radians at the end of the function for use in the code)
    std::cout << "Try num:\n";
    is >> num ;
    std::cout << num << "\n";
    std::cout << "Try wavelength:\n";
    is >> wavelength ;
    std::cout << wavelength << "\n";
    std::cout << "Try macroArgTheta:\n";
    is >> macroArgTheta;
    std::cout << macroArgTheta << "\n";
    std::cout << "Try macroArgPhi:\n";
    is >> macroArgPhi;
    std::cout << macroArgPhi << "\n";
    std::cout << "Checking input stream arguments.\n";
    std::cout << "Does is fail?:" << is.fail() << "\n";

    if (is.fail()){
      // check for scintillation wavelength spectrum
      is.str(newValues.c_str());
      is.clear();
      std::string material;
      is >> num >> material >> macroArgTheta >> macroArgPhi;
      if (is.fail())
	      Log::Die("VertexGen_PhotonBeam: Incorrect vertex setting " + newValues);
      fMaterial = material;

      // get the scintillation wavelength spectrum
      DBLinkPtr loptics = DB::Get()->GetLink("OPTICS", material);
      std::vector<double> wlarr = loptics->GetDArray("SCINTILLATION_value1");
      std::vector<double> wlamp = loptics->GetDArray("SCINTILLATION_value2");
      for (unsigned i=0; i<wlarr.size(); i++)
	      wlarr[i] = CLHEP::hbarc * CLHEP::twopi / (wlarr[i] * CLHEP::nm);
      if(wlarr.front() > wlarr.back()){
	      reverse(wlarr.begin(), wlarr.end());
	      reverse(wlamp.begin(), wlamp.end());
      }
      for (unsigned i=1; i<wlarr.size(); i++)
	      if(wlarr[i-1] >= wlarr[i])
	        Log::Die("VertexGen_PhotonBeam: wavelengths out of order");

      // use a linear interpolator to get a uniform sampling with bin
      // size smaller than the smallest bin in order to use RandGeneral
      LinearInterp<double> energyInterp(wlarr, wlamp);
      double step = 1.0e9;
      for(int i=0; i<(int)wlarr.size()-1; i++)
	      step = fmin(step, wlarr[i+1]-wlarr[i]);
      step /= 2;
      int nbins = (int) ((energyInterp.Max()-energyInterp.Min()) / step) + 1;
      step = (energyInterp.Max() - energyInterp.Min()) / (nbins - 1);

      // get the oversampled array, small padding at ends to avoid range error
      double* energyResample = new double[nbins];
      energyResample[0] = energyInterp(energyInterp.Min() + step*1e-6);
      energyResample[nbins-1] = energyInterp(energyInterp.Max() - step*1e-6);
      for (int i=1; i<nbins-1; i++)
	      energyResample[i] = energyInterp(energyInterp.Min() + i * step);
      fMinEnergy = energyInterp.Min();
      fMaxEnergy = energyInterp.Max();

      if(fRndmEnergy)
	      delete fRndmEnergy;
      fRndmEnergy = new CLHEP::RandGeneral(energyResample, nbins);
    }
    else
      fEnergy = CLHEP::hbarc * CLHEP::twopi / (wavelength * CLHEP::nm);

    double exp = 0.0;
    is >> exp;
    if (exp < 0.0)
      Log::Die("VertexGen_PhotonBeam: Exponential time constant must be positive");

    fNumPhotons = num;
    fExpTime = exp;
    macroTheta = macroArgTheta/180*CLHEP::pi;
    macroPhi = macroArgPhi/180*CLHEP::pi;

  }

  G4String VertexGen_PhotonRay::GetState()
  {
    if(fRndmEnergy)
      return dformat("%d\t%s\t%f", fNumPhotons, fMaterial.c_str(), fExpTime);
    else
      return dformat("%d\t%f\t%f", fNumPhotons, fEnergy, fExpTime);
  }

} // namespace RAT
