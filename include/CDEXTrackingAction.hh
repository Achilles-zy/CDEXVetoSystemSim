#ifndef CDEXTrackingAction_h
#define CDEXTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
class G4Track;
class CDEXDetectorConstruction;
class CDEXEventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CDEXTrackingAction : public G4UserTrackingAction
{
  public:
    CDEXTrackingAction(CDEXEventAction*, CDEXDetectorConstruction* cons);
   ~CDEXTrackingAction();

  public:
    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);

    void RecordTrackPos(G4ThreeVector trackpos);
    void AddEdepTrack(G4double edepstep);

    G4int GetParticleIntType(G4String name) {
		G4int particletype;
		if (name == "opticalphoton") {
			particletype = 0;
		}
		else if (name == "gamma") {
			particletype = 1;
		}
		else if (name == "e-") {
			particletype = 2;
		}
		else if (name == "e+") {
			particletype = 3;
		}
		else if (name == "alpha") {
			particletype = 4;
		}
		else {
			particletype = 5;
		}
		return particletype;
    }


	G4int GetCreatorProcessIntType(G4String name) {
		G4int creatorprocesstype;
		if (name == "RadioactiveDecay"|| name == "RadioactiveDecayBase") {
			creatorprocesstype = 0;
		}
		else if (name == "conv") {
			creatorprocesstype = 1;
		}
		else if (name == "phot") {
			creatorprocesstype = 2;
		}
		else if (name == "compt") {
			creatorprocesstype = 3;
		}
		else {
			creatorprocesstype = 4;
		}
		return creatorprocesstype;
	}
    
  private:
      G4double EdepTrack;
      G4ThreeVector TrackPos;
      CDEXEventAction* CDEXEvent;
      CDEXDetectorConstruction* CDEXCons;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
