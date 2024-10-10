////////////////////////////////////////////////////////////////////////
// Class:       LightCaloSBND
// Plugin Type: analyzer (art v3_04_00)
// File:        LightCaloSBND_module.cc
//
// Generated at Wed Aug 12 11:28:18 2020 by Francisco Nicolas-Arnaldos using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <stdlib.h>
#include <fstream>
#include <string.h>

//Geometry
#include "larcore/Geometry/Geometry.h"
// LArSoft includes
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/SimListUtils.h"
#include "lardataobj/Simulation/sim.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardataobj/RecoBase/OpFlash.h"

//BackTrackers
//#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "sbncode/OpDet/PDMapAlg.h"


#include <unordered_map>

//Random numbers
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "nurandom/RandomUtils/NuRandomService.h"

//Semi-Analytical model
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

namespace opana {
  class LightCaloSBND;
}


class opana::LightCaloSBND : public art::EDAnalyzer {
public:
  explicit LightCaloSBND(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LightCaloSBND(LightCaloSBND const&) = delete;
  LightCaloSBND(LightCaloSBND&&) = delete;
  LightCaloSBND& operator=(LightCaloSBND const&) = delete;
  LightCaloSBND& operator=(LightCaloSBND&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  void ResetVariables();
  double DriftPositionToDriftTime(double x_pos);
  //Services
  art::ServiceHandle<geo::Geometry> geo;

  //PDS map
  std::unique_ptr<opdet::PDMapAlg> fPDSMapPtr;

  double fDriftVelocity;
  double fWirePlanePosition;
  double fElectronLifetime;

  //Semi-Analytical model
  std::unique_ptr<phot::SemiAnalyticalModel> fSemiModel;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;

  //Conf fhicl variables
  bool fUseRecoPE;
  bool fUseRecoCharge;
  bool fUseAverageVisibility;
  std::vector<std::string> fSimPhotonsModuleLabel;
  std::string fSimEnergyDepositModuleLabel;
  std::string fSimEnergyDepositProvenanceLabel;
  std::vector<std::string> fOpFlashesModuleLabel;
  std::string fHitLabel;
  std::string fSpacePointLabel;
  double fDetectionEfficiency;
  std::vector<std::string> fPDTypes;

  // vectors for average visibilities
  std::vector<double> fAverageVisibilty;
  std::vector<double> fAverageVisibiltyReco;
  std::vector<double> fVisibilityMeanPoint;
  std::vector<double> fVisibilityMeanPointReco;

  // Declare member data here.
  TTree *fTree;
  unsigned int feventID, frunID, fsubrunID;

  // True mean position
  double fMeanX, fMeanY, fMeanZ;

  // Reco mean position
  double fMeanX_Reco, fMeanY_Reco, fMeanZ_Reco;

  // True and reconstructed photons
  std::vector <double> fTruePh_v;
  std::vector <double> fRecoPh_v;

  // Estimated photons from visibility
  std::vector <double> fHypoPhTrue_v;
  std::vector <double> fHypoPhRecoFlash_v;
  std::vector <double> fHypoPhRecoTrack_v;
  std::vector <double> fHypoPhReco_v;

  // Distance to PDs
  std::vector <double> fPDDistanceTrue_v;
  std::vector <double> fPDDistanceReco_v;

  // Resolutions
  std::vector <double> fHypoPhResolutionTrue_v;
  std::vector <double> fHypoPhResolutionRecoFlash_v;
  std::vector <double> fHypoPhResolutionRecoTrack_v;
  std::vector <double> fHypoPhResolutionReco_v;

  // produced photons and electrons for each energy deposition
  std::vector<double> fEnDep_v;
  std::vector<int> fEnDepPDG_v;
  std::vector<double> fEnDepStep_v;
  std::vector<double> fNScintPhotons_v;
  std::vector<double> fNDriftElectrons_v;

  double fDepEn;
  double fNScintPhotons;
  double fNDriftElectrons;
  double fRecoCharge;
  double fRecoChargeCorr;

  double fOpChAverageL;
};

opana::LightCaloSBND::LightCaloSBND(fhicl::ParameterSet const& p) :
  EDAnalyzer{p},
  fPDSMapPtr( art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDSMapTool")) ),
  fUseRecoPE ( p.get< bool >("UseRecoPE") ),
  fUseRecoCharge ( p.get< bool >("UseRecoCharge") ),
  fUseAverageVisibility ( p.get< bool >("UseAverageVisibility", true) ),
  fSimPhotonsModuleLabel ( p.get<std::vector<std::string>>("SimPhotonsModuleLabel",   {"pdfastsim", "pdfastsimout"}) ),
  fSimEnergyDepositModuleLabel ( p.get<std::string>("SimEnergyDepositModuleLabel",  "ionandscint") ),
  fSimEnergyDepositProvenanceLabel ( p.get<std::string>("SimEnergyDepositProvenanceLabel",  "priorSCE") ),
  fOpFlashesModuleLabel ( p.get<std::vector<std::string>>("OpFlashesModuleLabel",  {"opflashtpc0", "opflashtpc1"}) ),
  fHitLabel ( p.get<std::string>("HitLabel",  "gaushit") ),
  fSpacePointLabel ( p.get<std::string>("fSpacePointLabel",  "pandora") ),
  fDetectionEfficiency ( p.get< double >("DetectionEfficiency") ),
  fPDTypes ( p.get<std::vector<std::string>>("PDTypes",   {"pmt_coated", "pmt_uncoated"}) )
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  fWirePlanePosition = std::abs( geo->Plane(geo::PlaneID(0, 0, 0)).GetCenter().X() );
  fElectronLifetime = detProp.ElectronLifetime();

  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VISHits");
  fSemiModel = std::make_unique<phot::SemiAnalyticalModel>(_vuv_params, _vis_params, true, false);


  std::cout<<" LightCaloSBND module configured!\n";
  std::cout<<"  - Drift Velocity: "<<fDriftVelocity<<" cm/mus \n";
  std::cout<<"  - ElectronLifetime: "<<fElectronLifetime<<" mus \n";
}


void opana::LightCaloSBND::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("LightCaloSBND", "Optical Analyzer Tree");

  fTree->Branch("eventID", &feventID, "eventID/i");
  fTree->Branch("runID", &frunID, "runID/i");
  fTree->Branch("subrunID", &fsubrunID, "subrunID/i");

  fTree->Branch("DepEn", &fDepEn, "DepEn/d");
  fTree->Branch("NScintPhotons", &fNScintPhotons, "NScintPhotons/d");
  fTree->Branch("NDriftElectrons", &fNDriftElectrons, "NDriftElectrons/d");

  fTree->Branch("EnDep_v",&fEnDep_v);
  fTree->Branch("EnDepStep_v",&fEnDepStep_v);
  fTree->Branch("EnDepPDG_v",&fEnDepPDG_v);
  fTree->Branch("NScintPhotons_v",&fNScintPhotons_v);
  fTree->Branch("NDriftElectrons_v",&fNDriftElectrons_v);

  fTree->Branch("RecoCharge", &fRecoCharge, "RecoCharge/d");
  fTree->Branch("RecoChargeCorr", &fRecoChargeCorr, "RecoChargeCorr/d");

  fTree->Branch("MeanX", &fMeanX, "MeanX/d");
  fTree->Branch("MeanY", &fMeanY, "MeanY/d");
  fTree->Branch("MeanZ", &fMeanZ, "MeanZ/d");

  fTree->Branch("MeanX_Reco", &fMeanX_Reco, "MeanX_Reco/d");
  fTree->Branch("MeanY_Reco", &fMeanY_Reco, "MeanY_Reco/d");
  fTree->Branch("MeanZ_Reco", &fMeanZ_Reco, "MeanZ_Reco/d");

  fTree->Branch("TruePh_v",&fTruePh_v);
  fTree->Branch("RecoPh_v",&fRecoPh_v);

  fTree->Branch("HypoPhTrue_v",&fHypoPhTrue_v);
  fTree->Branch("HypoPhRecoFlash_v",&fHypoPhRecoFlash_v);
  fTree->Branch("HypoPhRecoTrack_v",&fHypoPhRecoTrack_v);
  fTree->Branch("HypoPhReco_v",&fHypoPhReco_v);

  fTree->Branch("PDDistanceTrue_v",&fPDDistanceTrue_v);
  fTree->Branch("PDDistanceReco_v",&fPDDistanceReco_v);

  fTree->Branch("HypoPhResolutionTrue_v",&fHypoPhResolutionTrue_v);
  fTree->Branch("HypoPhResolutionRecoFlash_v",&fHypoPhResolutionRecoFlash_v);
  fTree->Branch("HypoPhResolutionRecoTrack_v",&fHypoPhResolutionRecoTrack_v);
  fTree->Branch("HypoPhResolutionReco_v",&fHypoPhResolutionReco_v);

  fTree->Branch("OpChAverageL", &fOpChAverageL, "OpChAverageL/d");

}


void opana::LightCaloSBND::analyze(art::Event const& e)
{
  std::cout<<"Running LightCaloSBND---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<std::endl;

  //Event General Info
  feventID = e.id().event();
  frunID = e.id().run();
  fsubrunID = e.id().subRun();

  ResetVariables();

  // get the SimPhotons
  std::cout<<"Reading PE from SimPhotons..."<<std::endl;

  std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> fLitePhotonHandle_list;
  fLitePhotonHandle_list = e.getMany<std::vector<sim::SimPhotonsLite>>();


  // Get the true number of photons in each OpCh
  if( fLitePhotonHandle_list.size() == 0 )
    std::cout << "[Error] No available SimPhotonsLite" << std::endl;
  else{
    //loop over the Handle
    for ( const art::Handle<std::vector<sim::SimPhotonsLite>>& fLitePhotonHandle: (fLitePhotonHandle_list) ){
      std::string spLabel = fLitePhotonHandle.provenance()->moduleLabel();
      if(std::find(fSimPhotonsModuleLabel.begin(), fSimPhotonsModuleLabel.end(), spLabel)==fSimPhotonsModuleLabel.end()) continue;
      std::cout<<"*********Saving SimPhotonLite Handle: "<<fLitePhotonHandle.provenance()->moduleLabel()<<" "<<fLitePhotonHandle.provenance()->productInstanceName()<<std::endl;
      bool reflected = (fLitePhotonHandle.provenance()->productInstanceName()== "Reflected");
      for ( auto const& fLitePhotons : (*fLitePhotonHandle) ){
        std::string pd_type=fPDSMapPtr->pdType(fLitePhotons.OpChannel);
        if(reflected && pd_type=="xarapuca_vuv") continue;
        if(!reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;

        int opch=fLitePhotons.OpChannel;

        if(std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMapPtr->pdType(opch) ) != fPDTypes.end() ){
          std::map<int, int> fLitePhotons_map = fLitePhotons.DetectedPhotons;

          int nphotons=0;
          for(auto fphoton = fLitePhotons_map.begin(); fphoton!= fLitePhotons_map.end(); fphoton++){
            nphotons+=fphoton->second;
          }

          fTruePh_v.at(opch)+=nphotons;

        }

      }
    }
  }

  // Get the reconstructed number of photons in each OpCh
  if(fUseRecoPE){
    std::cout<<" Reading PE from OpFlashes...\n";
    art::Handle< std::vector<recob::OpFlash> > opflashListHandle;

    for (size_t s = 0; s < fOpFlashesModuleLabel.size(); s++) {
      std::cout<<"  --OpFlash Module Label:"<<fOpFlashesModuleLabel[s]<<std::endl;
      e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);

      for (unsigned int i = 0; i < opflashListHandle->size(); ++i) {
        // Get OpFlash
        art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
        recob::OpFlash Flash = *FlashPtr;

        std::vector<double> _flash_pe_v = Flash.PEs() ;

        for(size_t oc=0; oc<_flash_pe_v.size(); oc++){
          fRecoPh_v.at(oc)+=_flash_pe_v[oc];
        }

      }
    }
  }

  // Get the true track energy deposition points
  // and the average visibility for each PD
  art::Handle<std::vector<sim::SimEnergyDeposit> > SimEDHandle;
  std::vector<art::Ptr< sim::SimEnergyDeposit> > SimEDList;
  //if (e.getByLabel(fSimEnergyDepositLabel, fSimEnergyDepositInstanceLabel, SimEDHandle)) {
  if (e.getByLabel(fSimEnergyDepositModuleLabel, fSimEnergyDepositProvenanceLabel, SimEDHandle)) {
    art::fill_ptr_vector(SimEDList, SimEDHandle);
  }
  std::cout<<"  ---- Reading SimEnergyDeposition from handle: "<<SimEDHandle.provenance()->moduleLabel();
  std::cout<<":"<<SimEDHandle.provenance()->productInstanceName()<<" ----\n";

  for(auto &enedep:SimEDList){
    // Fill TTree vectors
    fEnDep_v.push_back(enedep->Energy());
    fEnDepStep_v.push_back(enedep->StepLength());
    fEnDepPDG_v.push_back(enedep->PdgCode());
    fNScintPhotons_v.push_back(enedep->NumPhotons());
    fNDriftElectrons_v.push_back(enedep->NumElectrons());

    fDepEn+=enedep->Energy();
    fNScintPhotons+=enedep->NumPhotons();
    fNDriftElectrons+=enedep->NumElectrons();
    fMeanX+=enedep->Energy()*enedep->MidPointX();
    fMeanY+=enedep->Energy()*enedep->MidPointY();
    fMeanZ+=enedep->Energy()*enedep->MidPointZ();

    geo::Point_t const p = {enedep->MidPointX(), enedep->MidPointY(), enedep->MidPointZ()};
    std::vector<double> point_visibility_vuv;
    std::vector<double> point_visibility_vis;
    fSemiModel->detectedDirectVisibilities(point_visibility_vuv, p);
    fSemiModel->detectedReflectedVisibilities(point_visibility_vis, p);

    for(size_t oc=0; oc<geo->NOpChannels(); oc++){
      std::string pd_type=fPDSMapPtr->pdType(oc);
      if(pd_type=="xarapuca_vuv")
        fAverageVisibilty[oc] = fAverageVisibilty[oc] + enedep->Energy()*point_visibility_vuv[oc];
      else if(pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")
        fAverageVisibilty[oc] = fAverageVisibilty[oc] + enedep->Energy()*point_visibility_vis[oc];
      else
        fAverageVisibilty[oc] = fAverageVisibilty[oc] + enedep->Energy()*(point_visibility_vuv[oc]+point_visibility_vis[oc]);
    }
  }
  fMeanX/=fDepEn;
  fMeanY/=fDepEn;
  fMeanZ/=fDepEn;

  for(size_t oc=0; oc<geo->NOpChannels(); oc++){
    fAverageVisibilty[oc] = fAverageVisibilty[oc]/fDepEn;
    std::cout<<fAverageVisibilty[oc]<<" ";
  }

  std::cout<<" Total Deposited Energy="<<fDepEn<<"  Total NScintPhotons="<<fNScintPhotons<<std::endl;
  std::cout<<" True XYZ: "<<fMeanX<<" "<<fMeanY<<" "<<fMeanZ<<std::endl;

  if(fUseRecoCharge){
    art::Handle<std::vector<recob::Hit>> hitsHandle;
    std::vector<art::Ptr<recob::Hit>> hitsVect;
    std::cout<<" --- Saving recob::Hit\n";
    e.getByLabel(fHitLabel, hitsHandle);
    art::fill_ptr_vector(hitsVect, hitsHandle);

    art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
    std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;
    std::cout<<" --- Saving recob::SpacePoints\n";

    e.getByLabel(fSpacePointLabel, eventSpacePoints);
    art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);

    art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, fSpacePointLabel);

    for (const art::Ptr<recob::SpacePoint> &SP: eventSpacePointsVect){

      std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());

      if (SPHit.at(0)->WireID().Plane==2){
        fRecoCharge+=SPHit.at(0)->Integral();

        double drift_time = DriftPositionToDriftTime( SP->position().X() ) ;
        double attenuation_corr = std::exp( drift_time/fElectronLifetime );

        std::cout<<" SPx="<< SP->position().X() << " DriftTime=" << drift_time <<" AttFactor="
        <<attenuation_corr<<std::endl;
        fRecoChargeCorr+=attenuation_corr * SPHit.at(0)->Integral();


        fMeanX_Reco+=SPHit.at(0)->Integral() * SP->position().X();
        fMeanY_Reco+=SPHit.at(0)->Integral() * SP->position().Y();
        fMeanZ_Reco+=SPHit.at(0)->Integral() * SP->position().Z();

        geo::Point_t const p = { SP->position().X(), SP->position().Y(),  SP->position().Z()};
        std::vector<double> point_visibility_vuv;
        std::vector<double> point_visibility_vis;
        fSemiModel->detectedDirectVisibilities(point_visibility_vuv, p);
        fSemiModel->detectedReflectedVisibilities(point_visibility_vis, p);

        for(size_t oc=0; oc<geo->NOpChannels(); oc++){
          std::string pd_type=fPDSMapPtr->pdType(oc);
          if(pd_type=="xarapuca_vuv")
            fAverageVisibiltyReco[oc] = fAverageVisibiltyReco[oc] + SPHit.at(0)->Integral()*point_visibility_vuv[oc];
          else if(pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")
            fAverageVisibiltyReco[oc] = fAverageVisibiltyReco[oc] + SPHit.at(0)->Integral()*point_visibility_vis[oc];
          else
            fAverageVisibiltyReco[oc] = fAverageVisibiltyReco[oc] + SPHit.at(0)->Integral()*(point_visibility_vuv[oc]+point_visibility_vis[oc]);
        }
      }


    }
    fMeanX_Reco/=fRecoCharge;
    fMeanY_Reco/=fRecoCharge;
    fMeanZ_Reco/=fRecoCharge;
    std::cout<<" Reco XYZ: "<<fMeanX_Reco<<" "<<fMeanY_Reco<<" "<<fMeanZ_Reco<<std::endl;

    for(size_t oc=0; oc<geo->NOpChannels(); oc++){
      fAverageVisibiltyReco[oc] = fAverageVisibiltyReco[oc]/fRecoCharge;
      std::cout<<fAverageVisibiltyReco[oc]<<" ";
    }
  }

  geo::Point_t const MeanXYZ = {fMeanX, fMeanY, fMeanZ};
  geo::Point_t const MeanXYZ_Reco = {fMeanX_Reco, fMeanY_Reco, fMeanZ_Reco};

  fSemiModel->detectedDirectVisibilities(fVisibilityMeanPoint, MeanXYZ);
  if(fUseRecoCharge){
    fSemiModel->detectedDirectVisibilities(fVisibilityMeanPointReco, MeanXYZ_Reco);
  }


  int opch_cont=0;
  for(size_t oc=0; oc<geo->NOpChannels(); oc++){
    if(fTruePh_v[oc]<=0) continue;

    geo::Point_t xyz_opch = geo->OpDetGeoFromOpChannel(oc).GetCenter();

    double d = std::sqrt( std::pow(xyz_opch.X()-fMeanX, 2) + std::pow(xyz_opch.Y()-fMeanY, 2) + std::pow(xyz_opch.Z()-fMeanZ, 2) );
    double dReco = std::sqrt( std::pow(xyz_opch.X()-fMeanX_Reco, 2) + std::pow(xyz_opch.Y()-fMeanY_Reco, 2) + std::pow(xyz_opch.Z()-fMeanZ_Reco, 2) );

    double hypope_true=0, hypope_reco=0, hypope_recoflash=0, hypope_recotrack=0;

    if(fUseAverageVisibility && fAverageVisibilty[oc]!=0){
      std::cout<<"Using Average Visibilities\n";
      hypope_true= fTruePh_v[oc]/(fAverageVisibilty[oc]*fDetectionEfficiency);

      if(fUseRecoCharge){
        hypope_recotrack = fTruePh_v[oc]/(fAverageVisibiltyReco[oc]*fDetectionEfficiency);
      }
      if(fUseRecoPE){
        hypope_recoflash = fRecoPh_v[oc]/(fAverageVisibilty[oc]*fDetectionEfficiency);
      }
      if(fUseRecoPE && fUseRecoCharge){
        hypope_reco = fRecoPh_v[oc]/(fAverageVisibiltyReco[oc]*fDetectionEfficiency);
      }
    }
    else if(fAverageVisibiltyReco[oc]!=0){
      std::cout<<"Using Visibilities from mean point\n";
      hypope_true= fTruePh_v[oc]/(fVisibilityMeanPoint[oc]*fDetectionEfficiency);

      if(fUseRecoCharge){
        hypope_recotrack = fTruePh_v[oc]/(fVisibilityMeanPointReco[oc]*fDetectionEfficiency);
      }
      if(fUseRecoPE){
        hypope_recoflash = fRecoPh_v[oc]/(fVisibilityMeanPoint[oc]*fDetectionEfficiency);
      }
      if(fUseRecoPE && fUseRecoCharge){
        hypope_reco = fRecoPh_v[oc]/(fVisibilityMeanPointReco[oc]*fDetectionEfficiency);
      }
    }

    std::cout <<" OpCh="<<oc<<" TruePh="<<fTruePh_v[oc]<<"   Hypo PE="<<hypope_reco<<std::endl;

    if(hypope_reco!=0){
      opch_cont++;
      fOpChAverageL+=hypope_reco;
    }

    fHypoPhTrue_v[oc] = hypope_true;
    fHypoPhResolutionTrue_v[oc] = (hypope_true-fNScintPhotons)/fNScintPhotons;

    fHypoPhRecoFlash_v[oc] = hypope_recoflash;
    fHypoPhResolutionRecoFlash_v[oc] = (hypope_recoflash-fNScintPhotons)/fNScintPhotons;

    fHypoPhRecoTrack_v[oc] = hypope_recotrack;
    fHypoPhResolutionRecoTrack_v[oc] = (hypope_recotrack-fNScintPhotons)/fNScintPhotons;

    fHypoPhReco_v[oc] = hypope_reco;
    fHypoPhResolutionReco_v[oc] = (hypope_reco-fNScintPhotons)/fNScintPhotons;

    fPDDistanceTrue_v[oc] = d;
    fPDDistanceReco_v[oc] = dReco;

  }
  fOpChAverageL/=opch_cont;


  //Fill Tree
  fTree->Fill();
}

void opana::LightCaloSBND::ResetVariables(){

  fMeanX=0;
  fMeanY=0;
  fMeanZ=0;

  fMeanX_Reco=0;
  fMeanY_Reco=0;
  fMeanZ_Reco=0;

  fTruePh_v.clear(); fTruePh_v.resize(geo->NOpChannels(), 0);
  fRecoPh_v.clear(); fRecoPh_v.resize(geo->NOpChannels(), 0);

  // Estimated photons from visibility
  fHypoPhTrue_v.clear(); fHypoPhTrue_v.resize(geo->NOpChannels(), 0);
  fHypoPhRecoFlash_v.clear(); fHypoPhRecoFlash_v.resize(geo->NOpChannels(), 0);
  fHypoPhRecoTrack_v.clear(); fHypoPhRecoTrack_v.resize(geo->NOpChannels(), 0);
  fHypoPhReco_v.clear(); fHypoPhReco_v.resize(geo->NOpChannels(), 0);

  // Distance to PDs
  fPDDistanceTrue_v.clear(); fPDDistanceTrue_v.resize(geo->NOpChannels(), 0);
  fPDDistanceReco_v.clear(); fPDDistanceReco_v.resize(geo->NOpChannels(), 0);

  // Resolutions
  fHypoPhResolutionTrue_v.clear(); fHypoPhResolutionTrue_v.resize(geo->NOpChannels(), 0);
  fHypoPhResolutionRecoFlash_v.clear(); fHypoPhResolutionRecoFlash_v.resize(geo->NOpChannels(), 0);
  fHypoPhResolutionRecoTrack_v.clear(); fHypoPhResolutionRecoTrack_v.resize(geo->NOpChannels(), 0);
  fHypoPhResolutionReco_v.clear(); fHypoPhResolutionReco_v.resize(geo->NOpChannels(), 0);

  fDepEn = 0.;
  fNScintPhotons = 0.;
  fNDriftElectrons = 0.;

  fRecoCharge = 0.;
  fRecoChargeCorr = 0.;

  fOpChAverageL = 0.;

  fAverageVisibilty.clear(); fAverageVisibilty.resize(geo->NOpChannels(), 0);
  fAverageVisibiltyReco.clear(); fAverageVisibiltyReco.resize(geo->NOpChannels(), 0);
  fVisibilityMeanPoint.clear(); fVisibilityMeanPointReco.resize(geo->NOpChannels(), 0);

  fEnDep_v.clear();
  fEnDepStep_v.clear();
  fEnDepPDG_v.clear();
  fNScintPhotons_v.clear();
  fNDriftElectrons_v.clear();

}


double opana::LightCaloSBND::DriftPositionToDriftTime(double x_pos){
  return ( fWirePlanePosition - std::abs(x_pos) )/fDriftVelocity;
}



DEFINE_ART_MODULE(opana::LightCaloSBND)
