#include "StFwdMip/StFwdMipAnalysisMaker.h"
#include "StFwdTrackMaker/Common.h"

#include "TMath.h"

#include <limits>
#include <map>
#include <string>
#include <string>
#include <vector>

#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StPrimaryVertex.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StTriggerData.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StChain/StChainOpt.h"

#include "StEventUtilities/StEventHelper.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"


#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"


#include "TROOT.h"
#include "TLorentzVector.h"
#include "StEvent/StFwdTrack.h"
#include "StFcsDbMaker/StFcsDb.h"

//________________________________________________________________________
StFwdMipAnalysisMaker::StFwdMipAnalysisMaker() : StMaker("fwdAna"){};
int StFwdMipAnalysisMaker::Finish() { 
    
    auto prevDir = gDirectory;
        
    // output file name
    string name = "FwdMip.root";
    TFile *fOutput = new TFile(name.c_str(), "RECREATE");
    fOutput->cd();
    for (auto nh : mHists) {
        nh.second->SetDirectory(gDirectory);
        nh.second->Write();
    }

    // restore previous directory
    gDirectory = prevDir;

    LOG_INFO << "Writing FwdAna output" << endm;

    m_tf_output_ttree -> cd() ;
    m_tt_output_ttree -> Write() ;
    m_tf_output_ttree -> Write() ;

    return kStOk; 
}
//________________________________________________________________________
int StFwdMipAnalysisMaker::Init() { 
    LOG_DEBUG << "StFwdMipAnalysisMaker::Init" << endm; 

    mHists["fwdMultFailed"] = new TH1F("fwdMultFailed", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultAll"] = new TH1F("fwdMultAll", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultGood"] = new TH1F("fwdMultGood", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultFST"] = new TH1F("fwdMultFST", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["nHitsFit"] = new TH1F("nHitsFit", ";nHitsFit; counts", 10, 0, 10);
    mHists["fwdMultEcalMatch"] = new TH1F("fwdMultEcalMatch", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultHcalMatch"] = new TH1F("fwdMultHcalMatch", ";N_{ch}^{FWD}; counts", 100, 0, 100);

    mHists["fwdMultEcalClusters"] = new TH1F("fwdMultEcalClusters", ";N_{Clu}^{ECAL}; counts", 100, 0, 100);
    mHists["fwdMultHcalClusters"] = new TH1F("fwdMultHcalClusters", ";N_{Clu}^{HCAL}; counts", 100, 0, 100);

    mHists["eta"] = new TH1F("eta", ";#eta; counts", 100, 0, 5);
    mHists["phi"] = new TH1F("phi", ";#phi; counts", 100, -3.1415926, 3.1415926);
    mHists["pt"] = new TH1F("pt", "; pT; counts", 500, 0, 10);

    mHists["ecalMatchPerTrack"] = new TH1F("ecalMatchPerTrack", ";N_{match} / track; counts", 5, 0, 5);
    mHists["hcalMatchPerTrack"] = new TH1F("hcalMatchPerTrack", ";N_{match} / track; counts", 5, 0, 5);

    mHists["matchedEcalEnergy"] = new TH1F("matchedEcalEnergy", ";Energy; counts", 100, 0, 15);
    mHists["matchedHcalEnergy"] = new TH1F("matchedHcalEnergy", ";Energy; counts", 100, 0, 15);

    mHists["ecalEnergy"] = new TH1F("ecalEnergy", ";Energy; counts", 100, 0, 15);
    mHists["hcalEnergy"] = new TH1F("hcalEnergy", ";Energy; counts", 100, 0, 15);

    mHists["ecalXY"] = new TH2F( "ecalXY", ";ecalX;ecalY", 200, -200, 200, 200, -200, 200 );
    mHists["hcalXY"] = new TH2F( "hcalXY", ";hcalX;hcalY", 200, 0, 50, 200, 0, 50 );

    mHists[ "ecaldX" ] = new TH1F( "ecaldX", ";dx (trk - ecal); counts", 400, -200, 200 );
    mHists[ "matchedEcaldX" ] = new TH1F( "matchedEcaldX", ";dx (trk - ecal); counts", 400, -200, 200 );
    mHists[ "ecaldY" ] = new TH1F( "ecaldY", ";dy (trk - ecal); counts", 400, -200, 200 );
    mHists[ "matchedEcaldY" ] = new TH1F( "matchedEcaldY", ";dy (trk - ecal); counts", 400, -200, 200 );
    mHists[ "ecaldR" ] = new TH1F( "ecaldR", ";dr (trk - ecal); counts", 400, 0, 400 );
    mHists[ "ecalMindR" ] = new TH1F( "ecalMindR", ";dr (trk - ecal); counts", 400, 0, 400 );
    mHists[ "matchedEcaldR" ] = new TH1F( "matchedEcaldR", ";dr (trk - ecal); counts", 400, 0, 400 );

    mHists[ "hcaldX" ] = new TH1F( "hcaldX", ";dx (trk - hcal); counts", 400, -200, 200 );
    mHists[ "matchedHcaldX" ] = new TH1F( "matchedHcaldX", ";dx (trk - hcal); counts", 400, -200, 200 );
    mHists[ "hcaldY" ] = new TH1F( "hcaldY", ";dy (trk - hcal); counts", 400, -200, 200 );
    mHists[ "matchedHcaldY" ] = new TH1F( "matchedHcaldY", ";dy (trk - hcal); counts", 400, -200, 200 );
    mHists[ "hcaldR" ] = new TH1F( "hcaldR", ";dr (trk - hcal); counts", 400, 0, 400 );
    mHists[ "hcalMindR" ] = new TH1F( "hcalMindR", ";dr (trk - hcal); counts", 400, 0, 400 );
    mHists[ "matchedHcaldR" ] = new TH1F( "matchedHcaldR", ";dr (trk - hcal); counts", 400, 0, 400 );

    mHists[ "trkEcalX" ] = new TH2F( "trkEcalX", ";trkX;ecalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkEcalY" ] = new TH2F( "trkEcalY", ";trkY;ecalY", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkEcalMinX" ] = new TH2F( "trkEcalMinX", ";trkX;ecalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkEcalMinY" ] = new TH2F( "trkEcalMinY", ";trkY;ecalY", 300, -150, 150, 300, -150, 150 );

    mHists[ "trkHcalX" ] = new TH2F( "trkHcalX", ";trkX;hcalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkHcalY" ] = new TH2F( "trkHcalY", ";trkY;hcalY", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkHcalMinX" ] = new TH2F( "trkHcalMinX", ";trkX;hcalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkHcalMinY" ] = new TH2F( "trkHcalMinY", ";trkY;hcalY", 300, -150, 150, 300, -150, 150 );



    //--owen: add a TTree
    m_tf_output_ttree = new TFile( "mip-analysis-ttree.root", "recreate" ) ;
    m_tt_output_ttree = new TTree( "mipanalysis", "FWD MIP analysis TTree" ) ;

    gROOT -> ProcessLine("#include <vector>") ; // owen: why is this necessary!?

//  m_tt_output_ttree -> Branch( "rcN",   &m_ttree_data. rcN, "rcN/I" ) ;
//  m_tt_output_ttree -> Branch( "rcPt",  &m_ttree_data. rcPt ) ;
//  m_tt_output_ttree -> Branch( "rcEta", &m_ttree_data. rcEta ) ;
//  m_tt_output_ttree -> Branch( "rcPhi", &m_ttree_data. rcPhi ) ;

//  m_tt_output_ttree -> Branch( "rcNumFST", &m_ttree_data. rcNumFST ) ;
//  m_tt_output_ttree -> Branch( "rcNumFTT", &m_ttree_data. rcNumFTT ) ;

//  m_tt_output_ttree -> Branch( "rcProjEcalx", &m_ttree_data. rcProjEcalx ) ;
//  m_tt_output_ttree -> Branch( "rcProjEcaly", &m_ttree_data. rcProjEcaly ) ;
//  m_tt_output_ttree -> Branch( "rcProjEcalz", &m_ttree_data. rcProjEcalz ) ;
//  m_tt_output_ttree -> Branch( "rcProjHcalx", &m_ttree_data. rcProjHcalx ) ;
//  m_tt_output_ttree -> Branch( "rcProjHcaly", &m_ttree_data. rcProjHcaly ) ;
//  m_tt_output_ttree -> Branch( "rcProjHcalz", &m_ttree_data. rcProjHcalz ) ;
//  m_tt_output_ttree -> Branch( "rcProjEcalPx", &m_ttree_data. rcProjEcalPx ) ;
//  m_tt_output_ttree -> Branch( "rcProjEcalPy", &m_ttree_data. rcProjEcalPy ) ;
//  m_tt_output_ttree -> Branch( "rcProjEcalPz", &m_ttree_data. rcProjEcalPz ) ;

//  m_tt_output_ttree -> Branch( "rcEcalClIndex", &m_ttree_data. rcEcalClIndex ) ;
//  m_tt_output_ttree -> Branch( "rcHcalClIndex", &m_ttree_data. rcHcalClIndex ) ;

//  m_tt_output_ttree -> Branch( "rcEcalClNhit", &m_ttree_data. rcEcalClNhit ) ;
//  m_tt_output_ttree -> Branch( "rcHcalClNhit", &m_ttree_data. rcHcalClNhit ) ;

//  m_tt_output_ttree -> Branch( "rcEcalClE", &m_ttree_data. rcEcalClE ) ;
//  m_tt_output_ttree -> Branch( "rcHcalClE", &m_ttree_data. rcHcalClE ) ;
//  m_tt_output_ttree -> Branch( "rcEcalClDx", &m_ttree_data. rcEcalClDx ) ;
//  m_tt_output_ttree -> Branch( "rcEcalClDy", &m_ttree_data. rcEcalClDy ) ;
//  m_tt_output_ttree -> Branch( "rcEcalClDr", &m_ttree_data. rcEcalClDr ) ;
//  m_tt_output_ttree -> Branch( "rcHcalClDx", &m_ttree_data. rcHcalClDx ) ;
//  m_tt_output_ttree -> Branch( "rcHcalClDy", &m_ttree_data. rcHcalClDy ) ;
//  m_tt_output_ttree -> Branch( "rcHcalClDr", &m_ttree_data. rcHcalClDr ) ;

//  m_tt_output_ttree -> Branch( "fcs_cl_ecalN", &m_ttree_data. fcs_cl_ecalN, "fcs_cl_ecalN/I" ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_ecalE", &m_ttree_data. fcs_cl_ecalE ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_ecalX", &m_ttree_data. fcs_cl_ecalX ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_ecalY", &m_ttree_data. fcs_cl_ecalY ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_ecalZ", &m_ttree_data. fcs_cl_ecalZ ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_ecalNhit", &m_ttree_data. fcs_cl_ecalNhit ) ;

//  m_tt_output_ttree -> Branch( "fcs_cl_hcalN", &m_ttree_data. fcs_cl_hcalN, "fcs_cl_hcalN/I" ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_hcalE", &m_ttree_data. fcs_cl_hcalE ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_hcalX", &m_ttree_data. fcs_cl_hcalX ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_hcalY", &m_ttree_data. fcs_cl_hcalY ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_hcalZ", &m_ttree_data. fcs_cl_hcalZ ) ;
//  m_tt_output_ttree -> Branch( "fcs_cl_hcalNhit", &m_ttree_data. fcs_cl_hcalNhit ) ;





//  m_tt_output_ttree ->Branch("fcs_rec_ecalN", &m_ttree_data.fcs_rec_ecalN, "fcs_rec_ecalN/I" ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalX", &m_ttree_data.fcs_rec_ecalX ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalY", &m_ttree_data.fcs_rec_ecalY ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalZ", &m_ttree_data.fcs_rec_ecalZ ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalLX", &m_ttree_data.fcs_rec_ecalLX ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalLY", &m_ttree_data.fcs_rec_ecalLY ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalE", &m_ttree_data.fcs_rec_ecalE ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalId", &m_ttree_data.fcs_rec_ecalId ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalDet", &m_ttree_data.fcs_rec_ecalDet ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalClIndex", &m_ttree_data.fcs_rec_ecalClIndex ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_ecalTrkIndex", &m_ttree_data.fcs_rec_ecalTrkIndex ) ;

//  m_tt_output_ttree ->Branch("fcs_rec_hcalN", &m_ttree_data.fcs_rec_hcalN, "fcs_rec_hcalN/I" ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalX", &m_ttree_data.fcs_rec_hcalX ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalY", &m_ttree_data.fcs_rec_hcalY ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalZ", &m_ttree_data.fcs_rec_hcalZ ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalLX", &m_ttree_data.fcs_rec_hcalLX ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalLY", &m_ttree_data.fcs_rec_hcalLY ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalE", &m_ttree_data.fcs_rec_hcalE ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalId", &m_ttree_data.fcs_rec_hcalId ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalDet", &m_ttree_data.fcs_rec_hcalDet ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalClIndex", &m_ttree_data.fcs_rec_hcalClIndex ) ;
//  m_tt_output_ttree ->Branch("fcs_rec_hcalTrkIndex", &m_ttree_data.fcs_rec_hcalTrkIndex ) ;


    m_tt_output_ttree ->Branch("n_mip", &m_ttree_data.n_mip, "n_mip/I" ) ;

    m_tt_output_ttree ->Branch("mip_ecal_E", &m_ttree_data.mip_ecal_E ) ;
    m_tt_output_ttree ->Branch("mip_ecal_nhit", &m_ttree_data.mip_ecal_nhit ) ;
    m_tt_output_ttree ->Branch("mip_ecal_x", &m_ttree_data.mip_ecal_x ) ;
    m_tt_output_ttree ->Branch("mip_ecal_y", &m_ttree_data.mip_ecal_y ) ;
    m_tt_output_ttree ->Branch("mip_ecal_iso20", &m_ttree_data.mip_ecal_iso20 ) ;
    m_tt_output_ttree ->Branch("mip_ecal_iso30", &m_ttree_data.mip_ecal_iso30 ) ;
    m_tt_output_ttree ->Branch("mip_trk_nsys", &m_ttree_data.mip_trk_nsys ) ;
    m_tt_output_ttree ->Branch("mip_trk_x", &m_ttree_data.mip_trk_x ) ;
    m_tt_output_ttree ->Branch("mip_trk_y", &m_ttree_data.mip_trk_y ) ;
    m_tt_output_ttree ->Branch("mip_hcal_E", &m_ttree_data.mip_hcal_E ) ;
    m_tt_output_ttree ->Branch("mip_hcal_nhit", &m_ttree_data.mip_hcal_nhit ) ;
    m_tt_output_ttree ->Branch("mip_hcal_x", &m_ttree_data.mip_hcal_x ) ;
    m_tt_output_ttree ->Branch("mip_hcal_y", &m_ttree_data.mip_hcal_y ) ;
    m_tt_output_ttree ->Branch("mip_hcal_iso20", &m_ttree_data.mip_hcal_iso20 ) ;
    m_tt_output_ttree ->Branch("mip_hcal_iso30", &m_ttree_data.mip_hcal_iso30 ) ;
    m_tt_output_ttree ->Branch("mip_hcal_ecal_cos_angle", &m_ttree_data.mip_hcal_ecal_cos_angle ) ;
    m_tt_output_ttree ->Branch("mip_ecal_hit0ind", &m_ttree_data.mip_ecal_hit0ind ) ;
    m_tt_output_ttree ->Branch("mip_ecal_hit1ind", &m_ttree_data.mip_ecal_hit1ind ) ;
    m_tt_output_ttree ->Branch("mip_hcal_hit0ind", &m_ttree_data.mip_hcal_hit0ind ) ;
    m_tt_output_ttree ->Branch("mip_hcal_hit1ind", &m_ttree_data.mip_hcal_hit1ind ) ;
    m_tt_output_ttree ->Branch("mip_hcal_hit2ind", &m_ttree_data.mip_hcal_hit2ind ) ;
    m_tt_output_ttree ->Branch("mip_hcal_hit3ind", &m_ttree_data.mip_hcal_hit3ind ) ;
    m_tt_output_ttree ->Branch("mip_hcal_hit4ind", &m_ttree_data.mip_hcal_hit4ind ) ;

    m_tt_output_ttree ->Branch("mip_trk_ind", &m_ttree_data.mip_trk_ind ) ;



    m_tt_output_ttree -> Print("") ;

    m_tt_output_ttree -> SetAutoFlush(0) ;

    //quiet = false ;
    quiet = true ;

    return kStOK;
}
//________________________________________________________________________
int StFwdMipAnalysisMaker::Make() {
    LOG_DEBUG << "StFwdMipAnalysisMaker::Make" << endm;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (event){
        StFttCollection *fttCol = event->fttCollection();
        if (fttCol){
            if ( ! quiet ) LOG_INFO << "The Ftt Collection has " << fttCol->numberOfPoints() << " points" << endm;
        }
    }
    long long itStart = FwdTrackerUtils::nowNanoSecond();
    ProcessFwdTracks();
    // ProcessFwdMuTracks();
    LOG_DEBUG << "Processing Fwd Tracks took: " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e6 << " ms" << endm;
    return kStOK;
} // Make
//________________________________________________________________________
void StFwdMipAnalysisMaker::Clear(const Option_t *opts) {

   LOG_DEBUG << "StFwdMipAnalysisMaker::CLEAR" << endm;

   m_ttree_data. rcPt.clear() ;
   m_ttree_data. rcEta.clear() ;
   m_ttree_data. rcPhi.clear() ;
   m_ttree_data. rcNumFST.clear() ;
   m_ttree_data. rcNumFTT.clear() ;

   m_ttree_data. rcProjEcalx.clear() ;
   m_ttree_data. rcProjEcaly.clear() ;
   m_ttree_data. rcProjEcalz.clear() ;
   m_ttree_data. rcProjHcalx.clear() ;
   m_ttree_data. rcProjHcaly.clear() ;
   m_ttree_data. rcProjHcalz.clear() ;
   m_ttree_data. rcProjEcalPx.clear() ;
   m_ttree_data. rcProjEcalPy.clear() ;
   m_ttree_data. rcProjEcalPz.clear() ;

   m_ttree_data. rcEcalClIndex.clear() ;
   m_ttree_data. rcHcalClIndex.clear() ;
   m_ttree_data. rcEcalClNhit.clear() ;
   m_ttree_data. rcHcalClNhit.clear() ;
   m_ttree_data. rcEcalClE.clear() ;
   m_ttree_data. rcHcalClE.clear() ;
   m_ttree_data. rcEcalClDx.clear() ;
   m_ttree_data. rcEcalClDy.clear() ;
   m_ttree_data. rcEcalClDr.clear() ;
   m_ttree_data. rcHcalClDx.clear() ;
   m_ttree_data. rcHcalClDy.clear() ;
   m_ttree_data. rcHcalClDr.clear() ;

   m_ttree_data. fcs_cl_ecalE.clear() ;
   m_ttree_data. fcs_cl_ecalX.clear() ;
   m_ttree_data. fcs_cl_ecalY.clear() ;
   m_ttree_data. fcs_cl_ecalZ.clear() ;
   m_ttree_data. fcs_cl_ecalNhit.clear() ;

   m_ttree_data. fcs_cl_hcalE.clear() ;
   m_ttree_data. fcs_cl_hcalX.clear() ;
   m_ttree_data. fcs_cl_hcalY.clear() ;
   m_ttree_data. fcs_cl_hcalZ.clear() ;
   m_ttree_data. fcs_cl_hcalNhit.clear() ;


   m_ttree_data. fcs_rec_ecalX.clear() ;
   m_ttree_data. fcs_rec_ecalY.clear() ;
   m_ttree_data. fcs_rec_ecalZ.clear() ;
   m_ttree_data. fcs_rec_ecalLX.clear() ;
   m_ttree_data. fcs_rec_ecalLY.clear() ;
   m_ttree_data. fcs_rec_ecalE.clear() ;
   m_ttree_data. fcs_rec_ecalId.clear() ;
   m_ttree_data. fcs_rec_ecalDet.clear() ;
   m_ttree_data. fcs_rec_ecalClIndex.clear() ;
   m_ttree_data. fcs_rec_ecalTrkIndex.clear() ;

   m_ttree_data. fcs_rec_hcalX.clear() ;
   m_ttree_data. fcs_rec_hcalY.clear() ;
   m_ttree_data. fcs_rec_hcalZ.clear() ;
   m_ttree_data. fcs_rec_hcalLX.clear() ;
   m_ttree_data. fcs_rec_hcalLY.clear() ;
   m_ttree_data. fcs_rec_hcalE.clear() ;
   m_ttree_data. fcs_rec_hcalId.clear() ;
   m_ttree_data. fcs_rec_hcalDet.clear() ;
   m_ttree_data. fcs_rec_hcalClIndex.clear() ;
   m_ttree_data. fcs_rec_hcalTrkIndex.clear() ;

   m_ttree_data. rcN = 0 ;
   m_ttree_data. fcs_rec_ecalN = 0 ;
   m_ttree_data. fcs_rec_hcalN = 0 ;
   m_ttree_data. fcs_cl_ecalN = 0 ;
   m_ttree_data. fcs_cl_hcalN = 0 ;

   m_ttree_data. n_mip = 0 ;

   m_ttree_data. mip_ecal_E.clear() ;
   m_ttree_data. mip_ecal_nhit.clear() ;
   m_ttree_data. mip_ecal_x.clear() ;
   m_ttree_data. mip_ecal_y.clear() ;
   m_ttree_data. mip_ecal_iso20.clear() ;
   m_ttree_data. mip_ecal_iso30.clear() ;
   m_ttree_data. mip_trk_nsys.clear() ;
   m_ttree_data. mip_trk_x.clear() ;
   m_ttree_data. mip_trk_y.clear() ;
   m_ttree_data. mip_hcal_E.clear() ;
   m_ttree_data. mip_hcal_nhit.clear() ;
   m_ttree_data. mip_hcal_x.clear() ;
   m_ttree_data. mip_hcal_y.clear() ;
   m_ttree_data. mip_hcal_iso20.clear() ;
   m_ttree_data. mip_hcal_iso30.clear() ;
   m_ttree_data. mip_hcal_ecal_cos_angle.clear() ;

   m_ttree_data. mip_ecal_hit0ind.clear() ;
   m_ttree_data. mip_ecal_hit1ind.clear() ;

   m_ttree_data. mip_hcal_hit0ind.clear() ;
   m_ttree_data. mip_hcal_hit1ind.clear() ;
   m_ttree_data. mip_hcal_hit2ind.clear() ;
   m_ttree_data. mip_hcal_hit3ind.clear() ;
   m_ttree_data. mip_hcal_hit4ind.clear() ;

   m_ttree_data. mip_trk_ind.clear() ;

}
//________________________________________________________________________
void StFwdMipAnalysisMaker::ProcessFwdTracks(  ){

    gMessMgr -> SwitchOff("I") ;
    gMessMgr -> SwitchOff("Q") ;
    gMessMgr -> SwitchOff("W") ;

    gMessMgr -> SetLimit("I",0) ;
    ///int info_limit = gMessMgr -> GetLimit("I") ;
    ///printf("  info limit is: %d\n", info_limit ) ;

    bool check_info = gMessMgr -> isInfoEnabled() ;
    if ( check_info ) {
       ///printf("  Info is enabled.\n") ;
    } else {
       ///printf("  Info is disabled.\n") ;
    }



    // This is an example of how to process fwd track collection
    if ( ! quiet ) LOG_INFO << "\n\n\nStFwdMipAnalysisMaker::ProcessFwdTracks" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if ( ! quiet ) LOG_INFO << Form("StFwdMipAnalysisMaker::ProcessFwdTracks : stEvent pointer %p ", stEvent) << endm ;
    if ( stEvent == 0x0 ) {
       if ( ! quiet ) LOG_INFO << "StFwdMipAnalysisMaker::ProcessFwdTracks : *** stEvent null pointer! " << endm ;
    }
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if ( ! quiet ) LOG_INFO << Form("StFwdMipAnalysisMaker::ProcessFwdTracks : stEvent->fwdTrackCollection() pointer %p ", ftc) << endm ;
    if ( ftc == 0x0 ) {
       if ( ! quiet ) LOG_INFO << "StFwdMipAnalysisMaker::ProcessFwdTracks : *** stEvent->fwdTrackCollection() returned null pointer! " << endm ;
    }
    if (!ftc)
        return;

    if ( ! quiet ) LOG_INFO << "Checking FcsCollection" << endm;
    StFcsCollection *fcs = stEvent->fcsCollection();
    if (!fcs) return;

    StFcsDb *mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));

    // if (ftc->tracks().size() > 4) return;

    size_t fwdMultEcalMatch = 0;
    size_t fwdMultHcalMatch = 0;
    size_t fwdMultFST = 0;

    if ( ! quiet ) LOG_INFO << "FwdTrackCollection has: " << ftc->tracks().size() << " tracks" << endm;


    mHists[ "fwdMultAll" ]->Fill( ftc->tracks().size() );

    // Cluster info (independen t of tracks)
    size_t fwdMultEcalClusters = 0;
    size_t fwdMultHcalClusters = 0;
    vector< StFcsCluster* > ecal_cl_pointers ;
    vector< StFcsCluster* > hcal_cl_pointers ;
    for ( int iDet = 0; iDet < 4; iDet++ ){
        for( size_t i = 0; i < fcs->clusters(iDet).size(); i++){
            StFcsCluster * clu = fcs->clusters(iDet)[i];

            StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());

            if ( iDet < 2 ) {
               m_ttree_data. fcs_cl_ecalE.push_back( clu->energy() ) ;
               m_ttree_data. fcs_cl_ecalX.push_back( xyz.x() ) ;
               m_ttree_data. fcs_cl_ecalY.push_back( xyz.y() ) ;
               m_ttree_data. fcs_cl_ecalZ.push_back( xyz.z() ) ;
               m_ttree_data. fcs_cl_ecalNhit.push_back( clu->hits().size() ) ;
               ecal_cl_pointers.push_back( clu ) ;
            } else {
               m_ttree_data. fcs_cl_hcalE.push_back( clu->energy() ) ;
               m_ttree_data. fcs_cl_hcalX.push_back( xyz.x() ) ;
               m_ttree_data. fcs_cl_hcalY.push_back( xyz.y() ) ;
               m_ttree_data. fcs_cl_hcalZ.push_back( xyz.z() ) ;
               m_ttree_data. fcs_cl_hcalNhit.push_back( clu->hits().size() ) ;
               hcal_cl_pointers.push_back( clu ) ;
            }
            if ( ! quiet ) LOG_INFO << Form( "  iDet = %d, cluster %2lu, cluster id = %2d,  E = %7.3f, X = %7.3f, Y = %7.3f, Z = %7.3f, Nhit = %lu",
                   iDet, i, clu->id(), clu->energy(), xyz.x(), xyz.y(), xyz.z(), clu->hits().size() ) << endm ;


            if ( iDet < 2 ){
                fwdMultEcalClusters++;
                mHists[ "ecalEnergy" ]->Fill( clu->energy() );
            } else if ( iDet < 4 ){
                fwdMultHcalClusters++;
                mHists[ "hcalEnergy" ]->Fill( clu->energy() );
            }
        }
    }
    m_ttree_data. fcs_cl_ecalN = fwdMultEcalClusters  ;
    m_ttree_data. fcs_cl_hcalN = fwdMultHcalClusters  ;

    mHists[ "fwdMultEcalClusters" ]->Fill( fwdMultEcalClusters );
    mHists[ "fwdMultHcalClusters" ]->Fill( fwdMultHcalClusters );


    size_t nGood = 0;
    size_t nFailed = 0;
    for ( auto fwdTrack : ftc->tracks() ){
        if ( !fwdTrack->didFitConvergeFully() ) {
            nFailed++;
            continue;
        }
        nGood++;


        if ( ! quiet ) LOG_INFO << TString::Format("StFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f ]", fwdTrack->mProjections.size(), fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size(), fwdTrack->momentum().perp()) << endm;
        
        if ( ! quiet ) LOG_INFO << "StFwdTrack has " << fwdTrack->ecalClusters().size() << " ecal matches" << endm;
        if ( ! quiet ) LOG_INFO << "StFwdTrack has " << fwdTrack->hcalClusters().size() << " hcal matches" << endm;

        mHists["ecalMatchPerTrack"]->Fill( fwdTrack->ecalClusters().size() );
        mHists["hcalMatchPerTrack"]->Fill( fwdTrack->hcalClusters().size() );
        
        mHists[ "nHitsFit" ]->Fill( fwdTrack->numberOfFitPoints() );

        if (fwdTrack->mFSTPoints.size() > 0){
            fwdMultFST ++;
        }


        //-- owen:
        for ( size_t ci=0; ci<fwdTrack->ecalClusters().size(); ci++ ) {
           StFcsCluster *clu = fwdTrack->ecalClusters()[ci];
           if ( ! quiet ) LOG_INFO << Form("    ECAL cluster %lu, energy = %7.3f,  ID = %d", ci, clu->energy(), clu->id() ) << endm ;
           if ( ! quiet ) LOG_INFO << Form("       cluster %lu has %lu ECAL hits", ci, clu->hits().size() ) << endm ;
           for ( size_t hi=0; hi<clu->hits().size(); hi++ ) {
              StFcsHit* hit = clu->hits()[hi] ;
              if ( ! quiet ) LOG_INFO << Form("         ECAL cluster %lu, id = %2d, hit %3lu :  detectorId = %d, id = %4d,  E = %7.3f", ci, clu->id(), hi, hit->detectorId(), hit->id(), hit->energy() ) << endm ;
           } // hi
        } // ci


        for ( size_t ci=0; ci<fwdTrack->hcalClusters().size(); ci++ ) {
           StFcsCluster *clu = fwdTrack->hcalClusters()[ci];
           if ( ! quiet ) LOG_INFO << Form("    HCAL cluster %lu, energy = %7.3f,  ID = %d", ci, clu->energy(), clu->id() ) << endm ;
           if ( ! quiet ) LOG_INFO << Form("       cluster %lu has %lu HCAL hits", ci, clu->hits().size() ) << endm ;
           for ( size_t hi=0; hi<clu->hits().size(); hi++ ) {
              StFcsHit* hit = clu->hits()[hi] ;
              if ( ! quiet ) LOG_INFO << Form("         HCAL cluster %lu, id = %2d, hit %3lu :  detectorId = %d, id = %4d,  E = %7.3f", ci, clu->id(), hi, hit->detectorId(), hit->id(), hit->energy() ) << endm ;
           } // hi
        } // ci


        mHists["eta"]->Fill( fwdTrack->momentum().pseudoRapidity() );
        mHists["phi"]->Fill( fwdTrack->momentum().phi() );
        mHists["pt"]->Fill( fwdTrack->momentum().perp() );
    

        if ( ! quiet ) LOG_INFO << Form("  Owen:  track %2lu :  Pt = %7.2f  Eta = %7.3f,  Phi = %7.3f\n", (nGood-1), fwdTrack->momentum().perp(), fwdTrack->momentum().pseudoRapidity(), fwdTrack->momentum().phi() ) << endm ;

        // ecal proj
        //////////float c[9];
        int detId = kFcsWcalId;
        TVector3 ecalXYZ;
        TVector3 ecapP;

        StFwdTrackProjection ecalProj = fwdTrack->getProjectionFor( detId, 0 );
        StFwdTrackProjection hcalProj = fwdTrack->getProjectionFor( kFcsHcalId, 0 );

        if ( ! quiet ) LOG_INFO << Form(" Owen:  ECAL detId = %d,   proj (x,y,z) = (%8.3f, %8.3f, %8.3f)", kFcsWcalId, ecalProj.mXYZ.x(), ecalProj.mXYZ.y(), ecalProj.mXYZ.z()  ) << endm ;
        if ( ! quiet ) LOG_INFO << Form(" Owen:  HCAL detId = %d,   proj (x,y,z) = (%8.3f, %8.3f, %8.3f)", kFcsHcalId, hcalProj.mXYZ.x(), hcalProj.mXYZ.y(), hcalProj.mXYZ.z()  ) << endm ;

        if ( ! quiet ) LOG_INFO << "EcalProj z= " << ecalProj.mXYZ.z() << endm;
        if ( ! quiet ) LOG_INFO << "HcalProj z= " << hcalProj.mXYZ.z() << endm;

        for ( size_t iEcal = 0; iEcal < fwdTrack->ecalClusters().size(); iEcal++ ){
            StFcsCluster *clu = fwdTrack->ecalClusters()[iEcal];
            if ( ! quiet ) LOG_INFO << "Ecal clu detId = " << clu->detectorId() << endm;
            mHists["matchedEcalEnergy"]->Fill( clu->energy() );

            StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());
            float dx = ecalProj.mXYZ.x() - xyz.x();
            float dy = ecalProj.mXYZ.y() - xyz.y();
            float dr = sqrt(dx*dx + dy*dy);
            mHists["matchedEcaldX"]->Fill( dx );
            mHists["matchedEcaldY"]->Fill( dy );
            mHists["matchedEcaldR"]->Fill( dr );
        }

        float ecal_cl_e = 0. ;
        float ecal_cl_dx = -999. ;
        float ecal_cl_dy = -999. ;
        float ecal_cl_dr = -9. ;
        int   ecal_cl_index = -1 ;
        int   ecal_cl_nhit = 0 ;

        if (ecalProj.mXYZ.z() > 500){
            double mindR = 999;
            StFcsCluster * cclu = nullptr; // closet cluster
            for ( int iDet = 0; iDet < 2; iDet++ ){
                for( size_t i = 0; i < fcs->clusters(iDet).size(); i++){
                    StFcsCluster * clu = fcs->clusters(iDet)[i];

                    StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());
                    mHists["ecalXY"]->Fill( xyz.x(), xyz.y() );

                    float dx = ecalProj.mXYZ.x() - xyz.x();
                    float dy = ecalProj.mXYZ.y() - xyz.y();
                    float dr = sqrt(dx*dx + dy*dy);

                    if ( fabs(dy) < 25 )
                        mHists[ "ecaldX" ]->Fill( dx );
                    if ( fabs(dx) < 25 )
                        mHists[ "ecaldY" ]->Fill( dy );
                    mHists[ "ecaldR" ]->Fill( dr );
                    if ( dr < mindR ){
                        mindR = dr;
                        cclu = clu;
                    }

                    mHists[ "trkEcalX" ] -> Fill( ecalProj.mXYZ.x(), xyz.x() );
                    mHists[ "trkEcalY" ] -> Fill( ecalProj.mXYZ.y(), xyz.y() );

                }
            }
            if ( cclu != nullptr ) {
               StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(cclu->detectorId(), cclu->x(), cclu->y());
               float dx = ecalProj.mXYZ.x() - xyz.x();
               float dy = ecalProj.mXYZ.y() - xyz.y();
               float dr = sqrt(dx*dx + dy*dy);
               ecal_cl_e = cclu->energy() ;
               ecal_cl_dx = dx ;
               ecal_cl_dy = dy ;
               ecal_cl_dr = dr ;
               for ( size_t pi=0; pi<ecal_cl_pointers.size(); pi++ ) {
                  if ( ecal_cl_pointers[pi] == cclu ) {
                     ecal_cl_index = pi ;
                     ecal_cl_nhit = cclu->hits().size() ;
                     if ( ! quiet ) LOG_INFO << Form("  ECAL cluster pointer %p is index %lu", cclu, pi ) << endm ;
                     break ;
                  }
               } // pi
            }



            mHists[ "ecalMindR" ]->Fill( mindR );
            if (cclu){
                StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(cclu->detectorId(), cclu->x(), cclu->y());
                mHists[ "trkEcalMinX" ] -> Fill( ecalProj.mXYZ.x(), xyz.x() );
                mHists[ "trkEcalMinY" ] -> Fill( ecalProj.mXYZ.y(), xyz.y() );
            }
        }

        float hcal_cl_e = 0. ;
        float hcal_cl_dx = -999. ;
        float hcal_cl_dy = -999. ;
        float hcal_cl_dr = -9. ;
        int   hcal_cl_index = -1 ;
        int   hcal_cl_nhit = 0 ;

        if (hcalProj.mXYZ.z() > 500){
            
            double mindR = 999;
            StFcsCluster * cclu = nullptr;
            for ( int iDet = 2; iDet < 4; iDet++ ){
                for( size_t i = 0; i < fcs->clusters(iDet).size(); i++){
                    StFcsCluster * clu = fcs->clusters(iDet)[i];
                    if (!clu) continue;
                    StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());
                    mHists["hcalXY"]->Fill( xyz.x(), xyz.y() );

                    float dx = hcalProj.mXYZ.x() - xyz.x();
                    float dy = hcalProj.mXYZ.y() - xyz.y();
                    float dr = sqrt(dx*dx + dy*dy);

                    if ( fabs(dy) < 25 )
                        mHists[ "hcaldX" ]->Fill( dx );
                    if ( fabs(dx) < 25 )
                        mHists[ "hcaldY" ]->Fill( dy );
                    mHists[ "hcaldR" ]->Fill( dr );

                    if ( dr < mindR ){
                        mindR = dr;
                        cclu = clu;
                    }

                    mHists[ "trkHcalX" ] -> Fill( hcalProj.mXYZ.x(), xyz.x() );
                    mHists[ "trkHcalY" ] -> Fill( hcalProj.mXYZ.y(), xyz.y() );
                }
            }

            if ( cclu != nullptr ) {
               StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(cclu->detectorId(), cclu->x(), cclu->y());
               float dx = hcalProj.mXYZ.x() - xyz.x();
               float dy = hcalProj.mXYZ.y() - xyz.y();
               float dr = sqrt(dx*dx + dy*dy);
               hcal_cl_e = cclu->energy() ;
               hcal_cl_dx = dx ;
               hcal_cl_dy = dy ;
               hcal_cl_dr = dr ;
               for ( size_t pi=0; pi<hcal_cl_pointers.size(); pi++ ) {
                  if ( hcal_cl_pointers[pi] == cclu ) {
                     hcal_cl_index = pi ;
                     hcal_cl_nhit = cclu->hits().size() ;
                     if ( ! quiet ) LOG_INFO << Form("  HCAL cluster pointer %p is index %lu", cclu, pi ) << endm ;
                     break ;
                  }
               } // pi
            }

            mHists[ "hcalMindR" ]->Fill( mindR );
            if (cclu){
                StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(cclu->detectorId(), cclu->x(), cclu->y());
                mHists[ "trkHcalMinX" ] -> Fill( hcalProj.mXYZ.x(), xyz.x() );
                mHists[ "trkHcalMinY" ] -> Fill( hcalProj.mXYZ.y(), xyz.y() );
            }
        }

        if (fwdTrack->ecalClusters().size() > 0)
            fwdMultEcalMatch++;
        if (fwdTrack->hcalClusters().size() > 0)
            fwdMultHcalMatch++;



       //-- TTree vector entries

        m_ttree_data.rcPt.push_back( fwdTrack->momentum().perp() ) ;
        m_ttree_data.rcEta.push_back( fwdTrack->momentum().pseudoRapidity() ) ;
        m_ttree_data.rcPhi.push_back( fwdTrack->momentum().phi() ) ;

        m_ttree_data.rcNumFST.push_back( fwdTrack->mFSTPoints.size() ) ;
        m_ttree_data.rcNumFTT.push_back( fwdTrack->mFTTPoints.size() ) ;

        m_ttree_data.rcProjEcalx.push_back( ecalProj.mXYZ.x() ) ;
        m_ttree_data.rcProjEcaly.push_back( ecalProj.mXYZ.y() ) ;
        m_ttree_data.rcProjEcalz.push_back( ecalProj.mXYZ.z() ) ;

        m_ttree_data.rcProjHcalx.push_back( hcalProj.mXYZ.x() ) ;
        m_ttree_data.rcProjHcaly.push_back( hcalProj.mXYZ.y() ) ;
        m_ttree_data.rcProjHcalz.push_back( hcalProj.mXYZ.z() ) ;

        m_ttree_data.rcProjEcalPx.push_back( ecalProj.mMom.x() ) ;
        m_ttree_data.rcProjEcalPy.push_back( ecalProj.mMom.y() ) ;
        m_ttree_data.rcProjEcalPz.push_back( ecalProj.mMom.z() ) ;

        m_ttree_data. rcEcalClE.push_back( ecal_cl_e ) ;
        m_ttree_data. rcEcalClDx.push_back( ecal_cl_dx ) ;
        m_ttree_data. rcEcalClDy.push_back( ecal_cl_dy ) ;
        m_ttree_data. rcEcalClDr.push_back( ecal_cl_dr ) ;
        m_ttree_data. rcEcalClIndex.push_back( ecal_cl_index ) ;
        m_ttree_data. rcEcalClNhit.push_back( ecal_cl_nhit ) ;

        m_ttree_data. rcHcalClE.push_back( hcal_cl_e ) ;
        m_ttree_data. rcHcalClDx.push_back( hcal_cl_dx ) ;
        m_ttree_data. rcHcalClDy.push_back( hcal_cl_dy ) ;
        m_ttree_data. rcHcalClDr.push_back( hcal_cl_dr ) ;
        m_ttree_data. rcHcalClIndex.push_back( hcal_cl_index ) ;
        m_ttree_data. rcHcalClNhit.push_back( hcal_cl_nhit ) ;


    } // Loop ftc->tracks()

    m_ttree_data.rcN = nGood ;

    mHists[ "fwdMultGood" ]->Fill( nGood );
    mHists[ "fwdMultFailed" ]->Fill( nFailed );
    mHists["fwdMultFST"]->Fill( fwdMultFST );
    mHists["fwdMultHcalMatch"]->Fill( fwdMultHcalMatch );
    mHists["fwdMultEcalMatch"]->Fill( fwdMultEcalMatch );

    if ( ! quiet ) LOG_INFO << "Found " << nFailed << " failed track fits out of " << ftc->tracks().size()  << endm;


   //--- owen: add FCS ecal and hcal hits.

    for ( int det = 0; det <4; det++ ) {
       StSPtrVecFcsHit&  hits = fcs->hits(det);
       if ( ! quiet ) LOG_INFO << Form("  Number of FCS hits:  %lu", hits.size()) << endm ;
       for ( unsigned int hi=0; hi<hits.size(); hi++ ) {
          StFcsHit* hit = hits[hi];
          float lx, ly ;
          float sx, sy, sz ;
          mFcsDb -> getLocalXYinCell( hit, lx, ly ) ;
          ///////////////StThreeVectorD s3v = mFcsDb -> getStarXYZ( hit->detectorId(), lx, ly ) ;
          StThreeVectorD s3v = mFcsDb -> getStarXYZfromColumnRow( hit->detectorId(), lx, ly ) ;
          sx = s3v.x() ;
          sy = s3v.y() ;
          sz = s3v.z() ;
          ////////  printf("  owen:  det=%d, id=%4d, ns=%d, ehp=%d, local x,y  (%7.4f, %7.4f), star x,y,z = (%7.4f, %7.4f, %7.4f)  E = %7.3f\n",
          ////////     hit->detectorId(), hit->id(), hit->ns(), hit->ehp(), lx, ly, sx, sy, sz, hit->energy()    ) ;
          if ( det < 2 ) {
             m_ttree_data.fcs_rec_ecalX.push_back( sx ) ;
             m_ttree_data.fcs_rec_ecalY.push_back( sy ) ;
             m_ttree_data.fcs_rec_ecalZ.push_back( sz ) ;
             m_ttree_data.fcs_rec_ecalLX.push_back( lx ) ;
             m_ttree_data.fcs_rec_ecalLY.push_back( ly ) ;
             m_ttree_data.fcs_rec_ecalE.push_back( hit->energy() ) ;
             m_ttree_data.fcs_rec_ecalId.push_back( hit->id() ) ;
             m_ttree_data.fcs_rec_ecalDet.push_back( hit->detectorId() ) ;
             m_ttree_data.fcs_rec_ecalN++ ;
             int cl_index = -1 ;
             //-- assuming a hit can be on only one cluster.
             for ( size_t pi=0; pi<ecal_cl_pointers.size(); pi++ ) {
                StFcsCluster* clu = ecal_cl_pointers[pi] ;
                for ( size_t hi=0; hi<clu->hits().size(); hi++ ) {
                   StFcsHit* thishit = clu->hits()[hi] ;
                   if ( thishit == hit ) {
                      cl_index = pi ;
                      break ;
                   }
                } // hi
                if ( cl_index >= 0 ) break ;
             } // pi
             m_ttree_data.fcs_rec_ecalClIndex.push_back( cl_index ) ;
             int trk_index = -1 ;
             if ( cl_index >= 0 ) {
                for ( size_t ti=0; ti < m_ttree_data. rcEcalClIndex.size(); ti ++ ) {
                   if ( m_ttree_data. rcEcalClIndex[ti] == cl_index ) {
                      trk_index = ti ;
                      break ;
                   }
                }
             }
             m_ttree_data.fcs_rec_ecalTrkIndex.push_back( trk_index ) ;
          } else {
             m_ttree_data.fcs_rec_hcalX.push_back( sx ) ;
             m_ttree_data.fcs_rec_hcalY.push_back( sy ) ;
             m_ttree_data.fcs_rec_hcalZ.push_back( sz ) ;
             m_ttree_data.fcs_rec_hcalLX.push_back( lx ) ;
             m_ttree_data.fcs_rec_hcalLY.push_back( ly ) ;
             m_ttree_data.fcs_rec_hcalE.push_back( hit->energy() ) ;
             m_ttree_data.fcs_rec_hcalId.push_back( hit->id() ) ;
             m_ttree_data.fcs_rec_hcalDet.push_back( hit->detectorId() ) ;
             m_ttree_data.fcs_rec_hcalN++ ;
             int cl_index = -1 ;
             //-- assuming a hit can be on only one cluster.
             for ( size_t pi=0; pi<hcal_cl_pointers.size(); pi++ ) {
                StFcsCluster* clu = hcal_cl_pointers[pi] ;
                for ( size_t hi=0; hi<clu->hits().size(); hi++ ) {
                   StFcsHit* thishit = clu->hits()[hi] ;
                   if ( thishit == hit ) {
                      cl_index = pi ;
                      break ;
                   }
                } // hi
                if ( cl_index >= 0 ) break ;
             } // pi
             m_ttree_data.fcs_rec_hcalClIndex.push_back( cl_index ) ;
             int trk_index = -1 ;
             if ( cl_index >= 0 ) {
                for ( size_t ti=0; ti < m_ttree_data. rcHcalClIndex.size(); ti ++ ) {
                   if ( m_ttree_data. rcHcalClIndex[ti] == cl_index ) {
                      trk_index = ti ;
                      break ;
                   }
                }
             }
             m_ttree_data.fcs_rec_hcalTrkIndex.push_back( trk_index ) ;
          }
       } // hi
    } // det








   //--- MIP analysis below here



    if ( ! quiet ) LOG_INFO << Form("\n\n Begin MIP analysis for this event.") << endm ;

    vector<int> all_used_ecal_hits ;

    vector< vector<int> > ecal_mip_clusters ;


    for ( unsigned int hi=0; hi<m_ttree_data.fcs_rec_ecalE.size(); hi++ ) {

       if ( std::find( all_used_ecal_hits.begin(), all_used_ecal_hits.end(), hi ) != all_used_ecal_hits.end() ) continue ; // already in list.  Skip

       all_used_ecal_hits.push_back( hi ) ; // add to all used hits list.

       float energy = m_ttree_data.fcs_rec_ecalE.at(hi) ;

       if ( energy < 0.001 ) continue ;

       vector<int> all_neighbors ;
       ecal_check_neighbors( hi, all_neighbors ) ;

       for ( unsigned int nhi=0; nhi<all_neighbors.size(); nhi++ ) { all_used_ecal_hits.push_back( all_neighbors.at( nhi ) ) ; }

       ecal_mip_clusters.push_back( all_neighbors ) ;

    } // hi
    if ( ! quiet ) LOG_INFO << Form(" Number of ECAL MIP cluster candidates: %lu", ecal_mip_clusters.size() ) << endm ;


    m_ttree_data. n_mip = 0 ;

    for ( unsigned int ci=0; ci< ecal_mip_clusters.size(); ci++ ) {

       vector<int>& cluster = ecal_mip_clusters.at(ci) ;

       float sumE = 0. ;
       float sumEx = 0. ;
       float sumEy = 0. ;

       for ( unsigned int hii=0; hii<cluster.size(); hii++ ) {
          int hi = cluster.at(hii) ;
          float x = m_ttree_data.fcs_rec_ecalX.at(hi) ;
          float y = m_ttree_data.fcs_rec_ecalY.at(hi) ;
          float energy = m_ttree_data.fcs_rec_ecalE.at(hi) ;
          //int id = m_ttree_data.fcs_rec_ecalId.at(hi) ;
          //int row = id/22 + 1 ; // correlated with Y
          //int col = id%22 + 1 ; // correlated with X
          //int det = m_ttree_data.fcs_rec_ecalDet.at(hi) ;
          //if (verbose) printf("  ci %3d :  hii %3d  hi %3d   (%3d, %3d) det %d   (%6.1f, %6.1f)   E = %6.3f\n", ci, hii, hi, row, col, det, x, y, energy  ) ;
          sumE += energy ;
          sumEx += energy * x ;
          sumEy += energy * y ;
       } // hii
       float ecal_clusterx = sumEx / sumE ;
       float ecal_clustery = sumEy / sumE ;
       int ecal_nhits = cluster.size() ;




      //--- calculate ECAL isolation

       float ecal_iso20 = 0. ;
       float ecal_iso30 = 0. ;

       for ( unsigned int hi=0; hi<m_ttree_data.fcs_rec_ecalE.size(); hi++ ) {

          if ( std::find( cluster.begin(), cluster.end(), hi ) != cluster.end() ) continue ; // skip hits that are in the cluster.

          float energy = m_ttree_data.fcs_rec_ecalE.at(hi) ;

          if ( energy < 0.001 ) continue ;

          float hitx = m_ttree_data.fcs_rec_ecalX.at(hi) ;
          float hity = m_ttree_data.fcs_rec_ecalY.at(hi) ;

          float dx = hitx - ecal_clusterx ;
          float dy = hity - ecal_clustery ;

          float dr = sqrt( dx*dx + dy*dy ) ;
          if ( dr < 20. ) ecal_iso20 += energy ;
          if ( dr < 30. ) ecal_iso30 += energy ;

       } // hi



      //--- find all tracks within 20 cm of ECAL cluster.

       int trk_ind(-1) ;
       int nmatch(0) ;

       for ( int ti=0; ti<m_ttree_data.rcN; ti++ ) {

          float tpx = m_ttree_data.rcProjEcalx.at(ti) ;
          float tpy = m_ttree_data.rcProjEcaly.at(ti) ;
          float dx = tpx - ecal_clusterx  ;
          float dy = tpy - ecal_clustery  ;
          float dr = sqrt( dx*dx + dy*dy ) ;

          if ( dr < 20. ) {
             nmatch ++ ;
             trk_ind = ti ;
          }

       } // ti

       if ( nmatch > 1 ) trk_ind = -2 ;



         int trk_nsys(0) ;
         float trk_x(-999.) ;
         float trk_y(-999.) ;
         if ( trk_ind >= 0 ) {
            if ( m_ttree_data.rcNumFST.at(trk_ind) > 0 ) trk_nsys++ ;
            if ( m_ttree_data.rcNumFTT.at(trk_ind) > 0 ) trk_nsys++ ;
            trk_x = m_ttree_data.rcProjEcalx.at(trk_ind) ;
            trk_y = m_ttree_data.rcProjEcaly.at(trk_ind) ;
         }




      //--- find closest HCAL hit to ECAL cluster.


       float min_dr = 999. ;
       int   closest_hcal_hit_index(-1) ;
       for ( unsigned int hi=0; hi<m_ttree_data.fcs_rec_hcalE.size(); hi++ ) {
          float energy_hcal = m_ttree_data.fcs_rec_hcalE.at(hi) ;
          if ( energy_hcal < 0.1 ) continue ;
          float x = m_ttree_data.fcs_rec_hcalX.at(hi) ;
          float y = m_ttree_data.fcs_rec_hcalY.at(hi) ;
          float rx = x - ecal_clusterx  ;
          float ry = y - ecal_clustery  ;
          float dr = sqrt( rx*rx + ry*ry ) ;
          if ( dr > 20 ) continue ;
          if ( dr < min_dr ) {
             min_dr = dr ;
             closest_hcal_hit_index = hi ;
          }
       } // hi




      //--- make the HCAL cluster associated with this ECAL cluster

       vector<int> all_neighbors_hcal ;
       hcal_check_neighbors( closest_hcal_hit_index, all_neighbors_hcal ) ;




      //--- calculate HCAL MIP cluster positon and energy

       int nhits_hcal(0) ;
       float sumE_hcal(0.00000000001) ;
       float sumEx_hcal(0.) ;
       float sumEy_hcal(0.) ;
       for ( unsigned int ni=0; ni<all_neighbors_hcal.size(); ni++ ) {
          int hi = all_neighbors_hcal.at(ni) ;
          float energy_hcal = m_ttree_data.fcs_rec_hcalE.at(hi) ;
          if ( energy_hcal < 0.01 ) continue ;
          float x = m_ttree_data.fcs_rec_hcalX.at(hi) ;
          float y = m_ttree_data.fcs_rec_hcalY.at(hi) ;
          nhits_hcal ++ ;
          sumE_hcal += energy_hcal ;
          sumEx_hcal += energy_hcal * x ;
          sumEy_hcal += energy_hcal * y ;
       } // hi
       float clusterx_hcal = sumEx_hcal / sumE_hcal ;
       float clustery_hcal = sumEy_hcal / sumE_hcal ;




      //--- calculate HCAL isolation

       float hcal_iso20(0.) ;
       float hcal_iso30(0.) ;
       for ( unsigned int hi=0; hi<m_ttree_data.fcs_rec_hcalE.size(); hi++ ) {
          if ( std::find( all_neighbors_hcal.begin(), all_neighbors_hcal.end(), hi ) != all_neighbors_hcal.end() ) continue ; // don't include hcal cluster h
          float energy_hcal = m_ttree_data.fcs_rec_hcalE.at(hi) ;
          if ( energy_hcal < 0.01 ) continue ;
          float x = m_ttree_data.fcs_rec_hcalX.at(hi) ;
          float y = m_ttree_data.fcs_rec_hcalY.at(hi) ;
          float dx = ecal_clusterx - x ;
          float dy = ecal_clustery - y ;
          float dr = sqrt( dx*dx + dy*dy ) ;
          if ( dr > 30 ) continue ;
          hcal_iso30 += energy_hcal ;
          if ( dr > 20 ) continue ;
          hcal_iso20 += energy_hcal ;
       } // hi





      //--- Calculate the directional consistency between the ECAL cluster direction
      //      and the HCAL cluster displacement relative to the ECAL cluster.

       float hcal_ecal_cos_angle(-2) ;
       if ( nhits_hcal > 0 && ecal_nhits > 0 ) {
          float cluster_dir_dx = clusterx_hcal - ecal_clusterx ;
          float cluster_dir_dy = clustery_hcal - ecal_clustery ;
          float cluster_dir_Lxy = sqrt( cluster_dir_dx*cluster_dir_dx + cluster_dir_dy*cluster_dir_dy  ) ;
          float cluster_ecal_Lxy = sqrt( ecal_clusterx*ecal_clusterx + ecal_clustery*ecal_clustery ) ;
          hcal_ecal_cos_angle = ( ecal_clusterx * cluster_dir_dx + ecal_clustery * cluster_dir_dy ) / ( cluster_ecal_Lxy * cluster_dir_Lxy ) ;
       }


       if ( ! quiet ) LOG_INFO << Form(" ECAL MIP cand %2u :  Nhit = %2d, E = %7.3f, trk_ind %3d  iso20 %7.3f, iso30 %7.3f, HCAL E = %7.3f, Nhit = %2d",
           ci, ecal_nhits, sumE, trk_ind, ecal_iso20, ecal_iso30, sumE_hcal, nhits_hcal ) << endm ;




      //--- save the hit indices.

       int tmp_ecal_hit0ind(-1) ;
       int tmp_ecal_hit1ind(-1) ;

       int tmp_hcal_hit0ind(-1) ;
       int tmp_hcal_hit1ind(-1) ;
       int tmp_hcal_hit2ind(-1) ;
       int tmp_hcal_hit3ind(-1) ;
       int tmp_hcal_hit4ind(-1) ;

       if ( cluster.size() > 0 ) tmp_ecal_hit0ind = cluster.at(0) ;
       if ( cluster.size() > 1 ) tmp_ecal_hit1ind = cluster.at(1) ;

       if ( all_neighbors_hcal.size() > 0 ) tmp_hcal_hit0ind = all_neighbors_hcal.at(0) ;
       if ( all_neighbors_hcal.size() > 1 ) tmp_hcal_hit1ind = all_neighbors_hcal.at(1) ;
       if ( all_neighbors_hcal.size() > 2 ) tmp_hcal_hit2ind = all_neighbors_hcal.at(2) ;
       if ( all_neighbors_hcal.size() > 3 ) tmp_hcal_hit3ind = all_neighbors_hcal.at(3) ;
       if ( all_neighbors_hcal.size() > 4 ) tmp_hcal_hit4ind = all_neighbors_hcal.at(4) ;



      //--- make entries for this ECAL cluster in the output TTree vectors.

      //--- add some selection requirements to avoid bloated output files.

       if (    sumE < 2.0
            && sumE > 0.05
            && ecal_nhits < 5   ) {

          m_ttree_data. mip_ecal_E.push_back( sumE ) ;
          m_ttree_data. mip_ecal_nhit.push_back( ecal_nhits ) ;
          m_ttree_data. mip_ecal_x.push_back( ecal_clusterx ) ;
          m_ttree_data. mip_ecal_y.push_back( ecal_clustery ) ;
          m_ttree_data. mip_ecal_iso20.push_back( ecal_iso20 ) ;
          m_ttree_data. mip_ecal_iso30.push_back( ecal_iso30 ) ;
          m_ttree_data. mip_trk_nsys.push_back( trk_nsys ) ;
          m_ttree_data. mip_trk_x.push_back( trk_x ) ;
          m_ttree_data. mip_trk_y.push_back( trk_y ) ;
          m_ttree_data. mip_hcal_E.push_back( sumE_hcal ) ;
          m_ttree_data. mip_hcal_nhit.push_back( nhits_hcal ) ;
          m_ttree_data. mip_hcal_x.push_back( clusterx_hcal ) ;
          m_ttree_data. mip_hcal_y.push_back( clustery_hcal ) ;
          m_ttree_data. mip_hcal_iso20.push_back( hcal_iso20 ) ;
          m_ttree_data. mip_hcal_iso30.push_back( hcal_iso30 ) ;

          m_ttree_data. mip_hcal_ecal_cos_angle.push_back( hcal_ecal_cos_angle ) ;

          m_ttree_data. mip_ecal_hit0ind.push_back( tmp_ecal_hit0ind ) ;
          m_ttree_data. mip_ecal_hit1ind.push_back( tmp_ecal_hit1ind ) ;

          m_ttree_data. mip_hcal_hit0ind.push_back( tmp_hcal_hit0ind ) ;
          m_ttree_data. mip_hcal_hit1ind.push_back( tmp_hcal_hit1ind ) ;
          m_ttree_data. mip_hcal_hit2ind.push_back( tmp_hcal_hit2ind ) ;
          m_ttree_data. mip_hcal_hit3ind.push_back( tmp_hcal_hit3ind ) ;
          m_ttree_data. mip_hcal_hit4ind.push_back( tmp_hcal_hit4ind ) ;

          m_ttree_data. mip_trk_ind.push_back( trk_ind ) ;

          m_ttree_data. n_mip ++ ;

       }


    } // ci ( loop over ECAL MIP clusters )





   //--- Look at MC info, if present.

    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if ( g2t_track != 0x0 ) {
       for (int irow = 0; irow < g2t_track->GetNRows(); irow++) {
           g2t_track_st *track = (g2t_track_st *)g2t_track->At(irow);

           if (0 == track)
               continue;

           int track_id = track->id;
           float pt2 = track->p[0] * track->p[0] + track->p[1] * track->p[1];
           float pt = std::sqrt(pt2);
           float eta = track->eta;
           float phi = std::atan2(track->p[1], track->p[0]); //track->phi;
           int q = track->charge;
           int pid = track->ge_pid ;

           LOG_INFO << Form( "  MC track %3d :  PID %7d  pt = %7.3f, eta = %7.3f, phi = %7.3f\n", irow, pid, pt, eta, phi ) << endm ;

        } // irow
    }






    m_tt_output_ttree -> Fill() ;

} // ProcessFwdTracks

//________________________________________________________________________
void StFwdMipAnalysisMaker::ProcessFwdMuTracks(  ){
    // This is an example of how to process fwd track collection
    if ( ! quiet ) LOG_INFO << "StFwdMipAnalysisMaker::ProcessFwdMuTracks" << endm;
    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker) {
        LOG_WARN << " No MuDstMaker ... bye-bye" << endm;
        return;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        LOG_WARN << " No MuDst ... bye-bye" << endm;
        return;
    }
    StMuFwdTrackCollection * ftc = mMuDst->muFwdTrackCollection();
    if (!ftc) return;

    StMuFcsCollection *fcs = mMuDst->muFcsCollection();
    if (!fcs) return;

    cout << "Number of StMuFwdTracks: " << ftc->numberOfFwdTracks() << endl;

    if (ftc->numberOfFwdTracks() > 4 ) return;

    StFcsDb *mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));
    

    size_t fwdMultFST = 0;
    size_t fwdMultEcalMatch = 0;
    size_t fwdMultHcalMatch = 0;

    for ( size_t iTrack = 0; iTrack < ftc->numberOfFwdTracks(); iTrack++ ){
        StMuFwdTrack * muFwdTrack = ftc->getFwdTrack( iTrack );
        // LOG_DEBUG << TString::Format("StMuFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f ]", muFwdTrack->mProjections.size(), muFwdTrack->mFTTPoints.size(), muFwdTrack->mFSTPoints.size(), muFwdTrack->momentum().Pt()) << endm;

        if ( ! quiet ) LOG_INFO << "StMuFwdTrack has " << muFwdTrack->mEcalClusters.GetEntries() << " Ecal matched" << endm;
        if ( ! quiet ) LOG_INFO << "StMuFwdTrack has " << muFwdTrack->mHcalClusters.GetEntries() << " Hcal matched" << endm;

        mHists["eta"]->Fill( muFwdTrack->momentum().Eta() );
        mHists["phi"]->Fill( muFwdTrack->momentum().Phi() );

        if (muFwdTrack->mFSTPoints.size() > 0){
            fwdMultFST ++;
        }

        if (muFwdTrack->mEcalClusters.GetEntries() > 0)
            fwdMultEcalMatch++;
        if (muFwdTrack->mHcalClusters.GetEntries() > 0)
            fwdMultHcalMatch++;

        
        // ecal proj
        ////////////float c[9];
        int detId = kFcsWcalId;
        TVector3 ecalXYZ;
        TVector3 ecapP;

        StMuFwdTrackProjection ecalProj;
        bool foundEcalProj = muFwdTrack->getProjectionFor( detId, ecalProj, 0 );

        if (foundEcalProj){

            // LOG_INFO << "EcalProj z= " << ecalProj.mXYZ.Z() << endm;
            for( size_t i = 0; i < fcs->numberOfClusters(); i++){
                StMuFcsCluster * clu = fcs->getCluster(i);

                if ( clu->detectorId() > 1 ) continue;

                if ( clu->energy() < 1 ) continue;
                StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());

                float dx = ecalProj.mXYZ.X() - xyz.x();
                float dy = ecalProj.mXYZ.Y() - xyz.y();
                float dr = sqrt(dx*dx + dy*dy);

                mHists[ "ecaldX" ]->Fill( dx );
                mHists[ "ecaldY" ]->Fill( dy );
                mHists[ "ecaldR" ]->Fill( dr );

                mHists[ "trkEcalX" ] -> Fill( ecalProj.mXYZ.X(), xyz.x() );

            }

            
        }

        
        for ( int i = 0; i < muFwdTrack->mEcalClusters.GetEntries(); i++ ){
            auto c = (StMuFcsCluster*) muFwdTrack->mEcalClusters.At(i);
            if (!c) continue;
            mHists["ecalEnergy"]->Fill( c->energy() );
            
            if ( ! quiet ) LOG_INFO << "eCal Cluster detId = " << c->detectorId() << endm;
            StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(c->detectorId(), c->x(), c->y());
            mHists["ecalXY"]->Fill( xyz.x(), xyz.y() );

            if (foundEcalProj)
                mHists["matchedEcaldX"]->Fill( ecalProj.mXYZ.X() - xyz.x() );

            
        }

        mHists["ecalMatchPerTrack"]->Fill( muFwdTrack->mEcalClusters.GetEntries() );
        mHists["hcalMatchPerTrack"]->Fill( muFwdTrack->mHcalClusters.GetEntries() );

        for ( int i = 0; i < muFwdTrack->mHcalClusters.GetEntries(); i++ ){
            auto c = (StMuFcsCluster*) muFwdTrack->mHcalClusters.At(i);
            if (!c) continue;
            mHists["hcalEnergy"]->Fill( c->energy() );

            mHists["hcalXY"]->Fill( c->x(), c->y() );
        }

    }
    mHists["fwdMult"]->Fill( ftc->numberOfFwdTracks() );
    mHists["fwdMultFST"]->Fill( fwdMultFST );
    mHists["fwdMultHcalMatch"]->Fill( fwdMultHcalMatch );
    mHists["fwdMultEcalMatch"]->Fill( fwdMultEcalMatch );

}

//---------------------------------------------

void StFwdMipAnalysisMaker::ecal_check_neighbors( int hit_index, vector<int>& all_neighbors ) {

   if ( all_neighbors.size() == 0 ) all_neighbors.push_back( hit_index ) ; // add the hit itself if necessary.

   int hit_id = m_ttree_data.fcs_rec_ecalId.at(hit_index) ;
   int hit_row = hit_id/22 + 1 ; // correlated with Y
   int hit_col = hit_id%22 + 1 ; // correlated with X
   int hit_det = m_ttree_data.fcs_rec_ecalDet.at(hit_index) ;

   for ( unsigned int hi=0; hi<m_ttree_data.fcs_rec_ecalE.size(); hi++ ) {

      unsigned int ui_hit_index = hit_index ; // avoid compiler warning
      if ( hi == ui_hit_index ) continue ; // skip hit itself

      if ( std::find( all_neighbors.begin(), all_neighbors.end(), hi ) != all_neighbors.end() ) continue ; // already in neighbor list.

      int det = m_ttree_data.fcs_rec_ecalDet.at(hi) ;
      if ( det != hit_det ) continue ; // not in same side of FCS.



      float energy = m_ttree_data.fcs_rec_ecalE.at(hi) ;

      if ( energy < 0.001 ) continue ;

      int id = m_ttree_data.fcs_rec_ecalId.at(hi) ;
      int row = id/22 + 1 ; // correlated with Y
      int col = id%22 + 1 ; // correlated with X

      if ( abs( hit_row - row ) <=1   && abs( hit_col - col ) <= 1 ) {
         ///if ( verbose ) printf("  hit %3d (%3d,%3d) neighbors hit %3d (%3d,%3d)\n", hit_index, hit_row, hit_col, hi, row, col ) ;
         all_neighbors.push_back( hi ) ;
         ecal_check_neighbors( hi, all_neighbors ) ;
      }


   } // hi

} // ecal_check_neighbors

//---------------------------------------------


void StFwdMipAnalysisMaker::hcal_check_neighbors( int hit_index, vector<int>& all_neighbors ) {

   if ( hit_index < 0 ) return ;

   if ( all_neighbors.size() == 0 ) all_neighbors.push_back( hit_index ) ; // add the hit itself if necessary.

   if ( all_neighbors.size() >= 15 ) return ; // cut off the size of the cluster.

   int hit_id = m_ttree_data.fcs_rec_hcalId.at(hit_index) ;
   int hit_row = hit_id/13 + 1 ; // correlated with Y
   int hit_col = hit_id%13 + 1 ; // correlated with X
   int hit_det = m_ttree_data.fcs_rec_hcalDet.at(hit_index) ;

   for ( unsigned int hi=0; hi<m_ttree_data.fcs_rec_hcalE.size(); hi++ ) {

      unsigned int ui_hit_index = hit_index ; // avoid compiler warning
      if ( hi == ui_hit_index ) continue ; // skip hit itself

      if ( std::find( all_neighbors.begin(), all_neighbors.end(), hi ) != all_neighbors.end() ) continue ; // already in neighbor list.

      int det = m_ttree_data.fcs_rec_hcalDet.at(hi) ;
      if ( det != hit_det ) continue ; // not in same side of FCS.



      float energy = m_ttree_data.fcs_rec_hcalE.at(hi) ;

      if ( energy < 0.01 ) continue ;

      int id = m_ttree_data.fcs_rec_hcalId.at(hi) ;
      int row = id/13 + 1 ; // correlated with Y
      int col = id%13 + 1 ; // correlated with X

    //--- see if they share a side (no corner connections).
      if (  ( abs( hit_row - row ) == 1 && (hit_col == col) ) ||
            ( abs( hit_col - col ) == 1 && (hit_row == row) )      )  {
         //////if ( verbose ) printf("  HCAL hit %3d (%3d,%3d) neighbors hit %3d (%3d,%3d)\n", hit_index, hit_row, hit_col, hi, row, col ) ;
         if ( all_neighbors.size() >= 15 ) return ; // cut off the size of the cluster.
         all_neighbors.push_back( hi ) ;
         hcal_check_neighbors( hi, all_neighbors ) ;
      }

   } // hi

} // hcal_check_neighbors

//-----------------------------------------------------------------------------------------------------------




