set(MODULE_LIBRARIES
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities

    ROOT::Tree
    ROOT::Core
    
    larsim::Utils
    larsim::PhotonPropagation
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::MCBase
    lardataobj::Simulation
    nusimdata::SimulationBase
    nug4::ParticleNavigation

    sbnobj::Common_Reco
    sbndcode_OpDetSim

)

cet_build_plugin(LightCaloSBND art::Module SOURCE LightCaloSBND_module.cc LIBRARIES ${MODULE_LIBRARIES} )


install_fhicl()
install_source()