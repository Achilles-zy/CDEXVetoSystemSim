#ifndef CDEXPrimaryGeneratorAction_h
#define CDEXPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "CDEXDetectorConstruction.hh"
#include "CDEXPrimaryGeneratorMessenger.hh"

class G4GeneralParticleSource;
class G4Event;
class CDEXDetectorConstruction;
class CDEXPrimaryGeneratorMessenger;

class CDEXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

public:
    CDEXPrimaryGeneratorAction(CDEXDetectorConstruction*);
    ~CDEXPrimaryGeneratorAction();
    void GeneratePrimaries(G4Event* anEvent);
    G4double GetPrimaryE() {
        return PrimaryE;
    }

    G4String GetPrimaryName() {
        return PrimaryName;
    }

    void SetSrcType(G4String type) {
        SrcType = type;
    }

    void SetImprintID(G4int id) {
        ImprintID = id;
    }

    G4String GetSrcType() {
        return SrcType;
    }

private:
    G4double PrimaryE;
    G4double InitialE;
    G4String PrimaryName;
    G4String SrcType;
    G4String LEGENDSrcPos;
    G4int ImprintID;


    G4GeneralParticleSource* fCDEXGPS;
    CDEXDetectorConstruction* fDetCons;
    CDEXPrimaryGeneratorMessenger* fPrimaryMessenger;
};

#endif