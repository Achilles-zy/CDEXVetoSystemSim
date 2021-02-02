//CDEXTrackingAction.cc

#include "CDEXTrackingAction.hh"
#include "CDEXDetectorConstruction.hh"
#include "CDEXEventAction.hh"

#include "G4Event.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
CDEXTrackingAction::CDEXTrackingAction(CDEXEventAction* evt)
    :CDEXEvent(evt)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
CDEXTrackingAction::~CDEXTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void CDEXTrackingAction::PreUserTrackingAction(const G4Track* track)
{

    
    G4double charge = track->GetDefinition()->GetPDGCharge();
    G4int ID = track->GetTrackID();
    G4int parentID = track->GetParentID();
    //G4cout << "photoncnt = " << CDEXEvent->GetPhotonCnt() << G4endl;
    G4double trackTime = track->GetGlobalTime();
    //if (CDEXEvent->GetPhotonCnt() > 3) {
        //track->GetStep()->GetTrack()->SetTrackStatus(fStopAndKill);
    //}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void CDEXTrackingAction::PostUserTrackingAction(const G4Track* trk)
{
    const CDEXDetectorConstruction* detectorConstruction
        = static_cast<const CDEXDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    G4String particlename = trk->GetParticleDefinition()->GetParticleName();
    if (trk->GetVolume() == detectorConstruction->GetEnv() && particlename == "opticalphoton") {
        //CDEXEvent->CountEscapedPhoton(1);
        //G4cout << trk->GetVolume()->GetName() << G4endl;
        //G4cout << "Next" << G4endl;
    }

    //G4cout << trk->GetVolume()->GetName() << G4endl;
    //G4cout << detectorConstruction->GetEnv()->GetName() << G4endl;
    //G4cout << trk->GetVolume()->GetName() << G4endl;
    //G4cout << trk->GetNextVolume()->GetName() << G4endl;
    //G4cout << "" <<G4endl;
    /*
    auto next_volume = trk->GetNextVolume();
    for(G4int i=0; i<16; i++)
    {
        auto sipm = CDEXDetCons->GetSiPM(i);
        if(sipm == next_volume){
            CDEXEvent->AddToSiPM(i);
        }
    }
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
