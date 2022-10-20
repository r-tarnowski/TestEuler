#include <iostream>
#include <math.h>


const int SUCCESS_RETURN = 0;
const int ERROR_RETURN = -1;

const double pi = 3.141592653589793239;
const double TWOpi = 6.283185307179586478;
const double HALFpi = 1.570796326794896619;
const double NEGHALFpi = -HALFpi;
const double DtoR = 1.45329251994329577e-02;

using std::cout;
using std::endl;

// my additions starts here

struct AngularVelocity {
   double roll;
   double pitch;
   double yaw;
};

struct Orientation {
   double psi;
   double theta;
   double phi;
};

// my additions stops here


struct DblOrient {
   double psi;
   double theta;
   double phi;
};

struct DblRotRotate {
   double roll;
   double pitch;
   double yaw;
};

struct TrigValues {
   double sinPsi, sinTheta, sinPhi;
   double cosPsi, cosTheta, cosPhi;
   double sPsisPhi, sPsicPhi, cPsisPhi, cPsicPhi;
};

//faster than fabs:
#define ABS(x) ( (x) >= 0 ? (x) : ( -(x) ) )

//minimum significant rate = 1deg/5sec
const double  MIN_ROTATION_RATE = 0.2 * DtoR;

void printHeader() {
   cout << endl;
   cout << "===============================================" << endl;
   cout << "Tests of the Euler's angles computing" << endl;
   cout << "implemented by John Towers & Jack Hines" << endl;
   cout << "===============================================" << endl;
   cout << endl;
}

void printValues( const char * pTitle, const Orientation & orientation ) {
   cout << pTitle << " Values:" << endl;
   cout << "-------------------------" << endl;
   cout << "Orientation [deg] :" << endl;
   cout << "psi   = " << std::fixed << orientation.psi / DtoR << endl;
   cout << "theta = " << std::fixed << orientation.theta / DtoR << endl;
   cout << "phi   = " << std::fixed << orientation.phi / DtoR << endl;
}

void printDRParams( const AngularVelocity & angularVelocity ) {

   cout << "Dead Reckoning Parameters:"<< endl;
   cout << "--------------------------" << endl;
   cout << "Angular Velocity [deg/s]  :" << endl;
   cout << "roll  = " << std::fixed << angularVelocity.roll / DtoR << endl;
   cout << "pitch = " << std::fixed << angularVelocity.pitch / DtoR << endl;
   cout << "yaw   = " << std::fixed << angularVelocity.yaw / DtoR << endl;
}

int calcEulerAngles( const double deltaTime, const AngularVelocity & angularVelocity,
                     Orientation & orientation ) {


   if ( 0.0 == deltaTime ) {
      cout << ">>> calcEulerAngles: Invalid input parameter: deltaTime." << endl;
      return ERROR_RETURN;
   }
   //bool rotating = false;
   if ( ( ABS(angularVelocity.roll) >= MIN_ROTATION_RATE ) ||
        ( ABS(angularVelocity.pitch) >= MIN_ROTATION_RATE ) ||
        ( ABS(angularVelocity.yaw) >= MIN_ROTATION_RATE ) ) {
      cout << ">>> calcEulerAngles: object is rotating." << endl;
   }
   else {
      cout << ">>> calcEulerAngles: object is NOT rotating (angular velocities too small)." << endl;
      return SUCCESS_RETURN;
   }


   //Set up and compute intermediate variables
   cout << ">>> calcEulerAngles: Setting up and computing intermediate variables." << endl;
   // Convert all angular values to double precision
   cout << ">>> calcEulerAngles: Converting all angular values to double precision." << endl;
   DblOrient initOrient = { 0.0, 0.0, 0.0 };
   initOrient.psi   = orientation.psi;
   initOrient.theta = orientation.theta;
   initOrient.phi   = orientation.phi;

   //Compute all sine and cosine values for the Euler angles
   cout << ">>> calcEulerAngles: Computing all sine and cosine values for the Euler angles." << endl;
   TrigValues trigVals = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0};
   trigVals.sinPsi   = sin( initOrient.psi   );
   trigVals.sinTheta = sin( initOrient.theta );
   trigVals.sinPhi   = sin( initOrient.phi   );
   trigVals.cosPsi   = cos( initOrient.psi   );
   trigVals.cosTheta = cos( initOrient.theta );
   trigVals.cosPhi   = cos( initOrient.phi   );
   trigVals.sPsisPhi = trigVals.sinPsi * trigVals.sinPhi;
   trigVals.sPsicPhi = trigVals.sinPsi * trigVals.cosPhi;
   trigVals.cPsisPhi = trigVals.cosPsi * trigVals.sinPhi;
   trigVals.cPsicPhi = trigVals.cosPsi * trigVals.cosPhi;

   //Build rotation matrix for WCS to ECS
   cout << ">>> calcEulerAngles: Building rotation matrix for WCS to ECS." << endl;
   double WCStoECS1Mat[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
   WCStoECS1Mat[ 0 ][ 0 ] = trigVals.cosPsi * trigVals.cosTheta;
   WCStoECS1Mat[ 0 ][ 1 ] = trigVals.sinPsi * trigVals.cosTheta;
   WCStoECS1Mat[ 0 ][ 2 ] = - trigVals.sinTheta;
   WCStoECS1Mat[ 1 ][ 0 ] = ( trigVals.cPsisPhi * trigVals.sinTheta ) - trigVals.sPsicPhi;
   WCStoECS1Mat[ 1 ][ 1 ] = ( trigVals.sPsisPhi * trigVals.sinTheta ) + trigVals.cPsicPhi;
   WCStoECS1Mat[ 1 ][ 2 ] = trigVals.cosTheta * trigVals.sinPhi;
   WCStoECS1Mat[ 2 ][ 0 ] = ( trigVals.cPsicPhi * trigVals.sinTheta ) + trigVals.sPsisPhi;
   WCStoECS1Mat[ 2 ][ 1 ] = ( trigVals.sPsicPhi * trigVals.sinTheta ) - trigVals.cPsisPhi;
   WCStoECS1Mat[ 2 ][ 2 ] = trigVals.cosTheta * trigVals.cosPhi;

   //Convert rotation rates to double
   cout << ">>> calcEulerAngles: Converting rotation rates to double." << endl;
   DblRotRotate rotRate = { 0.0, 0.0, 0.0 };
   rotRate.roll = angularVelocity.roll;
   rotRate.pitch = angularVelocity.pitch;
   rotRate.yaw = angularVelocity.yaw;

   //Compute initial ECS to final ECS rotation matrix
   cout << ">>> calcEulerAngles: Computing initial ECS to final ECS rotation matrix." << endl;
   double rollSq = rotRate.roll * rotRate.roll;
   double rollPitch = rotRate.roll * rotRate.pitch;
   double rollYaw = rotRate.roll * rotRate.yaw;
   double pitchSq = rotRate.pitch * rotRate.pitch;
   double pitchYaw = rotRate.pitch * rotRate.yaw;
   double yawSq = rotRate.yaw * rotRate.yaw;
   double wMagSq = rollSq + pitchSq + yawSq;
   double wMag = sqrt( wMagSq );
   double wMagT = wMag * deltaTime;
   double cosWMagT = cos( wMagT );
   double term1 = ( 1.0 - cosWMagT ) / wMagSq;
   double term1RollPitch = term1 * rollPitch;
   double term1RollYaw = term1 * rollYaw;
   double term1PitchYaw = term1 * pitchYaw;
   double sinWMagT = sin( wMagT );
   double term3 = sinWMagT / wMag;
   double term3Roll = term3 * rotRate.roll;
   double term3Pitch = term3 * rotRate.pitch;
   double term3Yaw = term3 * rotRate.yaw;
   double ECS1toECS2Mat[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
   ECS1toECS2Mat[ 0 ][ 0 ] = ( term1 * rollSq ) + cosWMagT;
   ECS1toECS2Mat[ 0 ][ 1 ] = term1RollPitch + term3Yaw;
   ECS1toECS2Mat[ 0 ][ 2 ] = term1RollYaw - term3Pitch;
   ECS1toECS2Mat[ 1 ][ 0 ] = term1RollPitch - term3Yaw;
   ECS1toECS2Mat[ 1 ][ 1 ] = ( term1 * pitchSq ) + cosWMagT;
   ECS1toECS2Mat[ 1 ][ 2 ] = term1PitchYaw + term3Roll;
   ECS1toECS2Mat[ 2 ][ 0 ] = term1RollYaw + term3Pitch;
   ECS1toECS2Mat[ 2 ][ 1 ] = term1PitchYaw - term3Roll;
   ECS1toECS2Mat[ 2 ][ 2 ] = ( term1 * yawSq ) + cosWMagT;

   //Compute final Euler angles
   double WCStoECS2Mat[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
   DblOrient finalOrient = { 0.0, 0.0, 0.0 };
   cout << ">>> calcEulerAngles: Computing final Euler angles." << endl;

   //Compute theta
   cout << ">>> calcEulerAngles: Computing theta." << endl;
   WCStoECS2Mat[0][2] = ECS1toECS2Mat[0][0] * WCStoECS1Mat[0][2] +
                        ECS1toECS2Mat[0][1] * WCStoECS1Mat[1][2] +
                        ECS1toECS2Mat[0][2] * WCStoECS1Mat[2][2];
   //trap round-off domain errors
   if ( WCStoECS2Mat[0][2] >= 1.0 ) {
      finalOrient.theta = HALFpi;
   }
   else if ( WCStoECS2Mat[0][2] <= -1.0 ) {
      finalOrient.theta = NEGHALFpi;
   }
   else {
      finalOrient.theta = -asin( WCStoECS2Mat[0][2] );
   }

   double cosTheta = cos( finalOrient.theta );
   double acosArg = 0;
   if ( cosTheta != 0.0 ) {
      //normal case: theta is not + or - pi/2
      cout << ">>> calcEulerAngles: normal case: theta is not + or - pi/2." << endl;

      //Compute Psi
      cout << ">>> calcEulerAngles: Computing psi." << endl;
      WCStoECS2Mat[0][0] = ECS1toECS2Mat[0][0] * WCStoECS1Mat[0][0] +
                           ECS1toECS2Mat[0][1] * WCStoECS1Mat[1][0] +
                           ECS1toECS2Mat[0][2] * WCStoECS1Mat[2][0];
      acosArg =  WCStoECS2Mat[0][0] / cosTheta;
      //trap round-off domain errors
      if ( acosArg >= 1.0 ) {
         finalOrient.psi = 0.0;
      }
      else if ( acosArg <= -1.0 ) {
         finalOrient.psi = pi;
      }
      else {
         finalOrient.psi = acos( acosArg );
      }
      //Normalize Psi to range of + or - pi radians
      cout << ">>> calcEulerAngles: Normalizing Psi to range of + or - pi radians." << endl;
      WCStoECS2Mat[0][1] = ECS1toECS2Mat[0][0] * WCStoECS1Mat[0][1] +
                           ECS1toECS2Mat[0][1] * WCStoECS1Mat[1][1] +
                           ECS1toECS2Mat[0][2] * WCStoECS1Mat[2][1];
      if ( WCStoECS2Mat[0][1] < 0.0 ) {
         finalOrient.psi = -finalOrient.psi;
      }

      //Compute Phi
      cout << ">>> calcEulerAngles: Computing phi." << endl;
      WCStoECS2Mat[2][2] = ECS1toECS2Mat[2][0] * WCStoECS1Mat[0][2] +
                           ECS1toECS2Mat[2][1] * WCStoECS1Mat[1][2] +
                           ECS1toECS2Mat[2][2] * WCStoECS1Mat[2][2];
      acosArg = WCStoECS2Mat[2][2] / cosTheta;
      //trap round-off domain errors
      if ( acosArg >= 1.0 ) {
         finalOrient.phi = 0.0;
      }
      else if ( acosArg <= -1.0 ) {
         finalOrient.phi = pi;
      }
      else {
         finalOrient.phi = acos( acosArg );
      }
      //Normalize Phi to range of + or - pi radians
      cout << ">>> calcEulerAngles: Normalizing Phi to range of + or - pi radians." << endl;
      WCStoECS2Mat[1][2] = ECS1toECS2Mat[1][0] * WCStoECS1Mat[0][2] +
                           ECS1toECS2Mat[1][1] * WCStoECS1Mat[1][2] +
                           ECS1toECS2Mat[1][2] * WCStoECS1Mat[2][2];
      if ( WCStoECS2Mat[1][2] < 0.0 ) {
         finalOrient.phi = -finalOrient.phi;
      }
   }
   else {
      //special case: theta is + or - pi/2
      cout << ">>> calcEulerAngles: special case: theta is + or - pi/2." << endl;

      // Compute psi
      cout << ">>> calcEulerAngles: Computing psi." << endl;
      WCStoECS2Mat[1][1] = ECS1toECS2Mat[1][0] * WCStoECS1Mat[0][1] +
                           ECS1toECS2Mat[1][1] * WCStoECS1Mat[1][1] +
                           ECS1toECS2Mat[1][2] * WCStoECS1Mat[2][1];
      acosArg = WCStoECS2Mat[1][1];
      //trap round-off domain errors
      if ( acosArg >= 1.0 ) {
         finalOrient.psi = 0.0;
      }
      else if ( acosArg <= -1.0 ) {
         finalOrient.psi = pi;
      }
      else {
         finalOrient.psi = acos( acosArg );
      }
      //Normalize Psi to range of + or - pi radians
      cout << ">>> calcEulerAngles: Normalizing Psi to range of + or - pi radians." << endl;
      WCStoECS2Mat[1][0] = ECS1toECS2Mat[1][0] * WCStoECS1Mat[0][0] +
                           ECS1toECS2Mat[1][1] * WCStoECS1Mat[1][0] +
                           ECS1toECS2Mat[1][2] * WCStoECS1Mat[2][0];
      if ( WCStoECS2Mat[1][0] > 0.0 ) {
         finalOrient.psi = - finalOrient.psi;
      }

      //Set phi to zero - no further rotations required
      cout << ">>> calcEulerAngles: Setting phi to zero - no further rotations required." << endl;
      finalOrient.phi = 0.0;
   }

   //Update ESPDU Euler angle values
   orientation.psi   = finalOrient.psi;
   orientation.theta = finalOrient.theta;
   orientation.phi   = finalOrient.phi;

   return SUCCESS_RETURN;
}

int main( int argc, char *argv[] ) {
   printHeader();

   Orientation orientation = { 0.0, 0.0, 0.0 };
   AngularVelocity angularVelocity = { 0.0, 0.0, 0.0 };

   double deltaTime = 0.0;
   int retVal = SUCCESS_RETURN;

   cout << endl << "TEST 1:" << endl;
   orientation.psi = 15.0 * DtoR;
   orientation.theta = 20.0 * DtoR;
   orientation.phi = 25.0 * DtoR;

   angularVelocity.roll = 60.0 * DtoR;
   angularVelocity.pitch = 65.0 * DtoR;
   angularVelocity.yaw = 70.0 * DtoR;

   deltaTime = 0.5;

   cout << "Calling calcEulerAngles with parameters: deltaTime = " << deltaTime << endl;
   printValues( "Initial", orientation );
   printDRParams( angularVelocity );
   retVal = calcEulerAngles( deltaTime, angularVelocity, orientation );
   cout << "calcEulerAngles called with: parameter: deltaTime = " << deltaTime << " - returned: "
        << retVal << endl;
   printValues( "Final", orientation );

   cout << endl << "TEST 2:" << endl;
   orientation.psi = 15.0 * DtoR;
   orientation.theta = 20.0 * DtoR;
   orientation.phi = 25.0 * DtoR;

   angularVelocity.roll = -60.0 * DtoR;
   angularVelocity.pitch = -65.0 * DtoR;
   angularVelocity.yaw = -70.0 * DtoR;

   deltaTime = 0.5;

   cout << "Calling calcEulerAngles with parameters: deltaTime = " << deltaTime << endl;
   printValues( "Initial", orientation );
   printDRParams( angularVelocity );
   retVal = calcEulerAngles( deltaTime, angularVelocity, orientation );
   cout << "calcEulerAngles called with: parameter: deltaTime = " << deltaTime << " - returned: "
        << retVal << endl;
   printValues( "Final", orientation );

   cout << endl << "TEST 3:" << endl;
   orientation.psi = 15.0 * DtoR;
   orientation.theta = 0.0 * DtoR;
   orientation.phi = 0.0 * DtoR;

   angularVelocity.roll = 0.0 * DtoR;
   angularVelocity.pitch = 0.0 * DtoR;
   angularVelocity.yaw = 30.0 * DtoR;

   deltaTime = 1.0;

   cout << "Calling calcEulerAngles with parameters: deltaTime = " << deltaTime << endl;
   printValues( "Initial", orientation );
   printDRParams( angularVelocity );
   retVal = calcEulerAngles( deltaTime, angularVelocity, orientation );
   cout << "calcEulerAngles called with: parameter: deltaTime = " << deltaTime << " - returned: "
        << retVal << endl;
   printValues( "Final", orientation );

   cout << endl << "TEST 4:" << endl;
   orientation.psi = 0.0 * DtoR;
   orientation.theta = 15.0 * DtoR;
   orientation.phi = 0.0 * DtoR;

   angularVelocity.roll = 0.0 * DtoR;
   angularVelocity.pitch = 30.0 * DtoR;
   angularVelocity.yaw = 0.0 * DtoR;

   deltaTime = 1.0;

   cout << "Calling calcEulerAngles with parameters: deltaTime = " << deltaTime << endl;
   printValues( "Initial", orientation );
   printDRParams( angularVelocity );
   retVal = calcEulerAngles( deltaTime, angularVelocity, orientation );
   cout << "calcEulerAngles called with: parameter: deltaTime = " << deltaTime << " - returned: "
        << retVal << endl;
   printValues( "Final", orientation );

   cout << endl << "TEST 5:" << endl;
   orientation.psi = 0.0 * DtoR;
   orientation.theta = 0.0 * DtoR;
   orientation.phi = 15.0 * DtoR;

   angularVelocity.roll = 30.0 * DtoR;
   angularVelocity.pitch = 0.0 * DtoR;
   angularVelocity.yaw = 0.0 * DtoR;

   deltaTime = 1.0;

   cout << "Calling calcEulerAngles with parameters: deltaTime = " << deltaTime << endl;
   printValues( "Initial", orientation );
   printDRParams( angularVelocity );
   retVal = calcEulerAngles( deltaTime, angularVelocity, orientation );
   cout << "calcEulerAngles called with: parameter: deltaTime = " << deltaTime << " - returned: "
        << retVal << endl;
   printValues( "Final", orientation );

   cout << endl << "TEST 6:" << endl;
   orientation.psi = 0.0 * DtoR;
   orientation.theta = 0.0 * DtoR;
   orientation.phi = 0.0 * DtoR;

   angularVelocity.roll = 0.0 * DtoR;
   angularVelocity.pitch = 0.0 * DtoR;
   angularVelocity.yaw = 0.0 * DtoR;

   deltaTime = 1.0;

   cout << "Calling calcEulerAngles with parameters: deltaTime = " << deltaTime << endl;
   printValues( "Initial", orientation );
   printDRParams( angularVelocity );
   retVal = calcEulerAngles( deltaTime, angularVelocity, orientation );
   cout << "calcEulerAngles called with: parameter: deltaTime = " << deltaTime << " - returned: "
        << retVal << endl;
   printValues( "Final", orientation );

   return SUCCESS_RETURN;
}
