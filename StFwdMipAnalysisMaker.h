#ifndef ST_FWD_MIP_ANALYSIS_MAKER_H
#define ST_FWD_MIP_ANALYSIS_MAKER_H

#include "StChain/StMaker.h"
#include "TVector3.h"
// ROOT includes
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
// STL includes
#include <vector>
#include <memory>
#include <map>

class StFwdTrack;

struct fwd_mip_ttree_data {

   int rcN ;

   std::vector<float> rcPt, rcEta, rcPhi ;
   std::vector<int>   rcNumFST, rcNumFTT ;
   vector<float> rcProjEcalx, rcProjEcaly, rcProjEcalz ;
   vector<float> rcProjHcalx, rcProjHcaly, rcProjHcalz ;
   vector<float> rcProjEcalPx, rcProjEcalPy, rcProjEcalPz ;
   vector<int> rcEcalClIndex, rcHcalClIndex ;
   vector<int> rcEcalClNhit, rcHcalClNhit ;
   vector<float> rcEcalClE, rcHcalClE ;
   vector<float> rcEcalClDx, rcEcalClDy, rcEcalClDr ;
   vector<float> rcHcalClDx, rcHcalClDy, rcHcalClDr ;


   int fcs_cl_ecalN ;
   vector<float> fcs_cl_ecalE, fcs_cl_ecalX, fcs_cl_ecalY, fcs_cl_ecalZ ;
   vector<int> fcs_cl_ecalNhit ;

   int fcs_cl_hcalN ;
   vector<float> fcs_cl_hcalE, fcs_cl_hcalX, fcs_cl_hcalY, fcs_cl_hcalZ ;
   vector<int> fcs_cl_hcalNhit ;



   int fcs_rec_ecalN ;
   vector<float> fcs_rec_ecalX, fcs_rec_ecalY, fcs_rec_ecalZ, fcs_rec_ecalE ;
   vector<float> fcs_rec_ecalLX, fcs_rec_ecalLY ;
   vector<int> fcs_rec_ecalId ;
   vector<int> fcs_rec_ecalDet ;
   vector<int> fcs_rec_ecalClIndex ;
   vector<int> fcs_rec_ecalTrkIndex ;

   int fcs_rec_hcalN ;
   vector<float> fcs_rec_hcalX, fcs_rec_hcalY, fcs_rec_hcalZ, fcs_rec_hcalE ;
   vector<float> fcs_rec_hcalLX, fcs_rec_hcalLY ;
   vector<int> fcs_rec_hcalId ;
   vector<int> fcs_rec_hcalDet ;
   vector<int> fcs_rec_hcalClIndex ;
   vector<int> fcs_rec_hcalTrkIndex ;


  //-- new analysis stuff
   int n_mip ;
   vector<float> mip_ecal_E ;
   vector<int>   mip_ecal_nhit ;
   vector<float> mip_ecal_x ;
   vector<float> mip_ecal_y ;
   vector<float> mip_ecal_iso20 ;
   vector<float> mip_ecal_iso30 ;
   vector<int>   mip_trk_nsys ;
   vector<float> mip_trk_x ;
   vector<float> mip_trk_y ;
   vector<float> mip_hcal_E ;
   vector<int>   mip_hcal_nhit ;
   vector<float> mip_hcal_x ;
   vector<float> mip_hcal_y ;
   vector<float> mip_hcal_iso20 ;
   vector<float> mip_hcal_iso30 ;
   vector<float> mip_hcal_ecal_cos_angle ;

   vector<int>   mip_ecal_hit0ind ;
   vector<int>   mip_ecal_hit1ind ;

   vector<int>   mip_hcal_hit0ind ;
   vector<int>   mip_hcal_hit1ind ;
   vector<int>   mip_hcal_hit2ind ;
   vector<int>   mip_hcal_hit3ind ;
   vector<int>   mip_hcal_hit4ind ;

   vector<int>   mip_trk_ind ;

} ;


class StFwdMipAnalysisMaker : public StMaker {

    ClassDef(StFwdMipAnalysisMaker, 0);

  public:
    StFwdMipAnalysisMaker();
    ~StFwdMipAnalysisMaker(){/* nada */};

    int Init();
    int Finish();
    int Make();
    void Clear(const Option_t *opts = "");
    void ProcessFwdTracks();
    void ProcessFwdMuTracks();

    void setQuiet( bool quietVal = true ) { quiet = quietVal ;  printf("\n\n Owen : Set quiet to %d\n\n", quietVal ) ; }

  private:
  protected:

    bool quiet ;

    void ecal_check_neighbors( int hit_index, vector<int>& all_neighbors ) ;
    void hcal_check_neighbors( int hit_index, vector<int>& all_neighbors ) ;


    TFile* m_tf_output_ttree = nullptr ;
    TTree* m_tt_output_ttree = nullptr ;

    fwd_mip_ttree_data    m_ttree_data   ;


  std::map<TString, TH1*> mHists;
};

#endif
