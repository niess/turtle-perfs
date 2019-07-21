/*
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <http://unlicense.org>
 */

#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserSteppingAction.hh"
#include "globals.hh"
/* C interface(s) */
extern "C" {
#include "tictac.h"
}


class G4ParticleGun;
class G4Track;


class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
        G4VPhysicalVolume * Construct();
};


class PhysicsList : public G4VModularPhysicsList {
  public:
        PhysicsList(G4double cut);
        ~PhysicsList();

        void ConstructParticle();
        void ConstructProcess();
        void SetCuts();

  private:
        G4VPhysicsConstructor * fDecayPhysics;
        G4VPhysicsConstructor * fEmPhysics;
        G4VPhysicsConstructor * fExtraPhysics;
};


class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
  public:
        PrimaryGenerator();
        ~PrimaryGenerator();
        void GeneratePrimaries(G4Event * event);

  protected:
        G4ParticleGun * fParticleGun;
};


class StackingAction : public G4UserStackingAction
{
  public:
        G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track * track);
};


class RunAction : public G4UserRunAction
{
  public:
          RunAction();
          void BeginOfRunAction(const G4Run *);
          void EndOfRunAction(const G4Run *);

  private:
          struct tictac fTictac;
};


class SteppingAction : public G4UserSteppingAction
{
  public:
          void UserSteppingAction(const G4Step * step);
};
