//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

void myanalysis_event(
                    const char *inFile = "/gpfs01/star/pwg_tasks/FwdCalib/PROD/ZeroField-07-21-2023/st_physics_23072003_raw_6000001.event.root",
                    int n = 100000000,
                    std::string configFile = "StRoot/StFwdMip/macro/event_track.xml",
                    const char *geom = "y2023") {
    TString _chain;
    gSystem->Load( "libStarRoot.so" );

    cout << "LL0" << endl;
	gSystem->Load("libStarClassLibrary.so");
	gSystem->Load("libStarRoot.so");
	cout << "LL1" << endl;
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	cout << "LL2" << endl;

    bool rerun_tracking = false;

    string fwdTrackOpt = "fwdTrack";
    if (rerun_tracking == false)
        fwdTrackOpt = "";
        
    _chain = Form("in, %s, useXgeom, AgML, db, StEvent, %s", geom, fwdTrackOpt.c_str());
    // "in, y2023, useXgeom, AgML, db, StEvent, MakeEvent"
    
    // needed in this wonky spack environment 
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);
    // bfc( 50, "in, y2023, useXgeom, AgML, db, StEvent, MakeEvent", "/gpfs01/star/pwg_tasks/FwdCalib/PROD/zfa/out/st_physics_23072003_raw_1500007.event.root" );
    // return;
    //  Extra configuration  for the Forward Tracking
/// if (rerun_tracking) {
///     StFwdTrackMaker *fwdTrack = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
///     if ( fwdTrack ){ //if it is in the chain
///         fwdTrack->SetConfigFile( configFile );
///         // write debug histograms and ttree?
///         fwdTrack->SetGenerateTree( false );
///         fwdTrack->SetGenerateHistograms( false );
///         // write out wavefront OBJ files
///         fwdTrack->SetVisualize( false );
///     }
/// }

    
    

    gSystem->Load( "libStFcsDbMaker.so" ) ;
    StFcsDbMaker * fcsDb = new StFcsDbMaker();
    chain->AddMaker(fcsDb);
	fcsDb->SetDebug();

//  if (rerun_tracking){
//      gSystem->Load("StFcsTrackMatchMaker");
//      StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
//      match->setMaxDistance(6,10);
//      match->setFileName("fcstrk.root");
//      match->SetDebug();
//      chain->AddAfter( "fwdTrack", match);
//  }


    /////////gSystem->Load("StFwdUtils.so");
    /////////StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
    /////////chain->AddMaker(fwdAna);

    gSystem->Load("StFwdMip.so");
    StFwdMipAnalysisMaker * fwdMipAna = new StFwdMipAnalysisMaker();
    chain->AddMaker(fwdMipAna);

    // Initialize the chain
    chain->Init();

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {
        chain->Clear();
        if (kStOK != chain->Make())
            break;
    }
}
