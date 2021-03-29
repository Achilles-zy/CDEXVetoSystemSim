// CDEXSteppingAction.hh

#ifndef CDEXSteppingAction_h
#define CDEXSteppingAction_h 1

#include "G4UserSteppingAction.hh"

#include "G4Types.hh"

class CDEXDetectorConstruction;
class CDEXEventAction;
class CDEXRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CDEXSteppingAction : public G4UserSteppingAction
{
  public:
      CDEXSteppingAction(CDEXEventAction*, CDEXRunAction*, CDEXDetectorConstruction*);
      ~CDEXSteppingAction() {};

    void UserSteppingAction(const G4Step*);
    inline G4double GetEfficiency(G4double wavelength);

  private:

    CDEXEventAction*      CDEXEvent;
    CDEXRunAction* CDEXRun;
    CDEXDetectorConstruction* CDEXCons;
    G4int SignalSiPMCount;
    G4int ContainerSignalSiPMCount;
    G4bool EnableAcc;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
