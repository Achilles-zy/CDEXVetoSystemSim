#include "CDEXDetectorMessenger.hh"

#include "CDEXDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

#include "G4RunManager.hh"

class CDEXDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

CDEXDetectorMessenger::CDEXDetectorMessenger(CDEXDetectorConstruction* Det)
	:G4UImessenger(),
	fDetCons(Det),
	fCDEXDir(0),
	fDetDir(0),
	fMatDir(0),
	cmdSetWireType(0),
	cmdSetReflectorType(0),
	cmdSetConfine(0),
	cmdSetRunInfo(0),
	cmdSetLayerNb(0),
	cmdSetReadoutAngle(0),
	cmdSetPENPropertiesID(0),
	cmdSetOuterReflector(0),
	cmdSetInnerReflector(0)
{
	fDetDir = new G4UIdirectory("/CDEX/cons/set/");
	fDetDir->SetGuidance("Set construction parameters");

	fCDEXDir = new G4UIdirectory("/CDEX/sim/set/");
	fCDEXDir->SetGuidance("Set simulation parameters");

	fMatDir = new G4UIdirectory("/CDEX/mat/set/");
	fMatDir->SetGuidance("Set material parameters");

	cmdSetConfine = new G4UIcmdWithAString("/CDEX/sim/set/confine", this);
	cmdSetConfine->SetGuidance("Set confine name in file name.");
	cmdSetConfine->SetParameterName("ConfineInfo", false);
	cmdSetConfine->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetConfine->SetToBeBroadcasted(false);

	cmdSetRunInfo = new G4UIcmdWithAString("/CDEX/sim/set/runinfo", this);
	cmdSetRunInfo->SetGuidance("Add run info");
	cmdSetRunInfo->SetParameterName("RunInfo", false);
	cmdSetRunInfo->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetRunInfo->SetToBeBroadcasted(false);

	cmdSetMode = new G4UIcmdWithAString("/CDEX/sim/set/mode", this);
	cmdSetMode->SetGuidance("Select simulation mode.");
	cmdSetMode->SetParameterName("SimulationMode", false);
	cmdSetMode->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetMode->SetToBeBroadcasted(false);

	cmdSetWireType = new G4UIcmdWithAString("/CDEX/cons/set/wiretype", this);
	cmdSetWireType->SetGuidance("Select wire type.");
	cmdSetWireType->SetParameterName("WireType", false);
	cmdSetWireType->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetWireType->SetToBeBroadcasted(false);

	cmdSetReflectorType = new G4UIcmdWithAString("/CDEX/cons/set/reflectortype", this);
	cmdSetReflectorType->SetGuidance("Select reflector type.");
	cmdSetReflectorType->SetParameterName("ReflectorType", false);
	cmdSetReflectorType->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetReflectorType->SetToBeBroadcasted(false);

	cmdSetPENPropertiesID = new G4UIcmdWithAnInteger("/CDEX/mat/set/PENpropertiesID", this);
	cmdSetPENPropertiesID->SetGuidance("Set PEN properties ID");
	cmdSetPENPropertiesID->SetParameterName("ID", false);
	cmdSetPENPropertiesID->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetPENPropertiesID->SetToBeBroadcasted(false);

	cmdSetOuterReflector = new G4UIcmdWithABool("/CDEX/cons/set/outerreflector", this);
	cmdSetOuterReflector->SetGuidance("Enable/disable outer reflective surface on OuterPENShell.");
	cmdSetOuterReflector->SetParameterName("ifReflcetor", false);
	cmdSetOuterReflector->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetOuterReflector->SetToBeBroadcasted(false);

	cmdSetInnerReflector = new G4UIcmdWithABool("/CDEX/cons/set/innerreflector", this);
	cmdSetInnerReflector->SetGuidance("Enable/disable inner reflective surface on OuterPENShell.");
	cmdSetInnerReflector->SetParameterName("ifReflcetor", false);
	cmdSetInnerReflector->AvailableForStates(G4State_PreInit, G4State_Idle);
	cmdSetInnerReflector->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDEXDetectorMessenger::~CDEXDetectorMessenger()
{
	delete fDetDir;
	delete fCDEXDir;
	delete fMatDir;
	delete cmdSetWireType;
	delete cmdSetReflectorType;
	delete cmdSetConfine;
	delete cmdSetRunInfo;
	delete cmdSetPENPropertiesID;
	delete cmdSetOuterReflector;
	delete cmdSetInnerReflector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDEXDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if (command == cmdSetWireType) {
		fDetCons->SetWireType(newValue);
	}

	if (command == cmdSetReflectorType) {
		fDetCons->SetReflectorType(newValue);
	}

	if (command == cmdSetRunInfo) {
		fDetCons->SetRunInfo(newValue);
	}

	if (command == cmdSetConfine) {
		fDetCons->SetConfine(newValue);
	}

	if (command == cmdSetMode) {
		fDetCons->SetMode(newValue);
	}

	if (command == cmdSetPENPropertiesID) {
		fDetCons->SetPENPropertiesID(cmdSetPENPropertiesID->GetNewIntValue(newValue));
	}

	if (command == cmdSetOuterReflector) {
		fDetCons->SetOuterReflector(cmdSetOuterReflector->GetNewBoolValue(newValue));
	}

	if (command == cmdSetInnerReflector) {
		fDetCons->SetInnerReflector(cmdSetInnerReflector->GetNewBoolValue(newValue));
	}
}