#ifndef CDEXTrackingAction_h
#define CDEXTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class G4Track;
class CDEXDetectorConstruction;
class CDEXEventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CDEXTrackingAction : public G4UserTrackingAction
{
  public:
    CDEXTrackingAction(CDEXEventAction*);
   ~CDEXTrackingAction();

  public:
    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);

  private:

    CDEXEventAction*      CDEXEvent;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
