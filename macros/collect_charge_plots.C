/// collect_charge_plots
/// Blair Jamieson (Feb. 2020)
///
/// Make uniform set of plots from each PTF run
/// No inputs.  All files hard coded.
//#include "find_circle.hpp"

// --- PTF style ---
TStyle* SetPTFStyle(Int_t WhichStyle = 1, TString styleName = "PTF") {
  TStyle *ptfStyle= new TStyle(styleName, "PTF approved plots style");
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t FontStyle = 22;
  Float_t FontSizeLabel = 0.035;
  Float_t FontSizeTitle = 0.05;
  Float_t YOffsetTitle = 1.3;
  switch(WhichStyle) {
    case 1:
      FontStyle = 42;
      FontSizeLabel = 0.05;
      FontSizeTitle = 0.065;
      YOffsetTitle = 1.19;
      break;
    case 2:
      FontStyle = 42;
      FontSizeLabel = 0.035;
      FontSizeTitle = 0.05;
      YOffsetTitle = 1.6;
      break;
    case 3:
      FontStyle = 132;
      FontSizeLabel = 0.03;
      FontSizeTitle = 0.035;
      YOffsetTitle = 1.3;
      break;
  }
  // use plain black on white colors
  ptfStyle->SetFrameBorderMode(0);
  ptfStyle->SetCanvasBorderMode(0);
  ptfStyle->SetCanvasBorderSize(0);
  ptfStyle->SetPadBorderMode(0);
  ptfStyle->SetPadColor(0);
  ptfStyle->SetCanvasColor(0);
  ptfStyle->SetStatColor(0);
  ptfStyle->SetFillColor(0);
  ptfStyle->SetEndErrorSize(4);
  ptfStyle->SetStripDecimals(kFALSE);
  ptfStyle->SetLegendBorderSize(0);
  ptfStyle->SetLegendFont(FontStyle);
  // set the paper & margin sizes
  ptfStyle->SetPaperSize(20, 26);
  ptfStyle->SetPadTopMargin(0.1);
  ptfStyle->SetPadBottomMargin(0.15);
  ptfStyle->SetPadRightMargin(0.13);
  // 0.075 -> 0.13 for colz option
  ptfStyle->SetPadLeftMargin(0.16);
  //to include both large/small font options
  // Fonts, sizes, offsets
  ptfStyle->SetTextFont(FontStyle);
  ptfStyle->SetTextSize(0.08);
  ptfStyle->SetLabelFont(FontStyle, "x");
  ptfStyle->SetLabelFont(FontStyle, "y");
  ptfStyle->SetLabelFont(FontStyle, "z");
  ptfStyle->SetLabelFont(FontStyle, "t");
  ptfStyle->SetLabelSize(FontSizeLabel, "x");
  ptfStyle->SetLabelSize(FontSizeLabel, "y");
  ptfStyle->SetLabelSize(FontSizeLabel, "z");
  ptfStyle->SetLabelOffset(0.015, "x");
  ptfStyle->SetLabelOffset(0.015, "y");
  ptfStyle->SetLabelOffset(0.015, "z");
  ptfStyle->SetTitleFont(FontStyle, "x");
  ptfStyle->SetTitleFont(FontStyle, "y");
  ptfStyle->SetTitleFont(FontStyle, "z");
  ptfStyle->SetTitleFont(FontStyle, "t");
  ptfStyle->SetTitleSize(FontSizeTitle, "y");
  ptfStyle->SetTitleSize(FontSizeTitle, "x");
  ptfStyle->SetTitleSize(FontSizeTitle, "z");
  ptfStyle->SetTitleSize(0.05, "t");
  ptfStyle->SetTitleOffset(1.14, "x");
  ptfStyle->SetTitleOffset(YOffsetTitle, "y");
  ptfStyle->SetTitleOffset(1.2, "z");
  ptfStyle->SetTitleStyle(0);
  //ptfStyle->SetTitleFontSize(0.05);
  //0.08
  ptfStyle->SetTitleFont(FontStyle, "pad");
  ptfStyle->SetTitleBorderSize(0);
  ptfStyle->SetTitleX(0.1f);
  ptfStyle->SetTitleW(0.8f);
  ptfStyle->SetTitleY(0.93f);
  // use bold lines and markers
  ptfStyle->SetMarkerStyle(20);
  ptfStyle->SetHistLineWidth( Width_t(2.5) );
  ptfStyle->SetLineStyleString(2, "[12 12]");
  // postscript dashes
  // get rid of X error bars and y error bar caps
  ptfStyle->SetErrorX(0.001);
  // do not display any of the standard histogram decorations
  //ptfStyle->SetOptTitle(0);
  ptfStyle->SetOptStat(0);
  ptfStyle->SetOptFit(0);
  // put tick marks on top and RHS of plots
  ptfStyle->SetPadTickX(1);
  ptfStyle->SetPadTickY(1);
  // -- color --
  ptfStyle->SetFillColor(1);
  // make color fillings (not white)
  // - color setup for 2D -
  // - "cold"/ blue-ish -
  Double_t red[] = { 0.00, 0.00, 0.00 };
  Double_t green[] = { 1.00, 0.00, 0.00 };
  Double_t blue[] = { 1.00, 1.00, 0.25 };
  // - "warm" red-ish colors -
  // Double_t red[] = {1.00, 1.00, 0.25 };
  // Double_t green[] = {1.00, 0.00, 0.00 };
  // Double_t blue[] = {0.00, 0.00, 0.00 };
  Double_t stops[] = { 0.25, 0.75, 1.00 };
  const Int_t NRGBs = 3;
  const Int_t NCont = 500;
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  ptfStyle->SetNumberContours(NCont);
  // - inbuilt color schemes -
  // ptfStyle->SetPalette(1); // use the rainbow color set
  //ptfStyle->SetPalette(kViridis); // use the viridis color set
  ptfStyle->SetPalette(kSunset); // use the sunset color set
  //ptfStyle->SetPalette(kColorPrintableOnGrey); // use the colorPrintableOnGrey color set
  //ptfStyle->SetPalette(kCubehelix); // use the cubehelix color set
  //ptfStyle->SetPalette(kLightTemperature);
  TColor::InvertPalette(); // invert color palette

  return ptfStyle;
}



/// make_1d_from_2d
/// Take z-axis values of TH2D and use them to fill a newly constructed TH1D
/// Last four arguments of function are used to build the TH1D
/// Ignores values that are less than 0.001
/// Returns a pointer to the newly created TH1D
TH1D* make_1d_from_2d( TH2D* hin, std::string title, int nbins, double xmin, double xmax ){
  static int count=0;
  ++count;
  std::ostringstream os;
  os<<"h1d_from_2d_"<<count;
  TH1D* h = new TH1D(os.str().c_str(), title.c_str(), nbins, xmin, xmax );

  for (unsigned ix = 1; ix <= hin->GetNbinsX(); ++ix ){
    for (unsigned iy = 1; iy <= hin->GetNbinsY(); ++iy ){
      if ( fabs( hin->GetBinContent( ix, iy) ) > 0.001 ){
	h->Fill( hin->GetBinContent( ix, iy ) );
      }
    }
  } 
  return h;
}


void collect_charge_plots(){
  //std::vector< int > colors = { kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring-1, kTeal, kAzure, kViolet, kPink, kRed+2 };
  //gStyle->SetPalette( kInvertedDarkBodyRadiator );
  
  TString localStyleName = "PTF";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t localWhichStyle = 3;
  TStyle* ptfstyle = SetPTFStyle(localWhichStyle, localStyleName);
  gROOT->SetStyle(ptfstyle->GetName());

  // -- margin --
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.15); //more space if no colz option
  gStyle->SetPadLeftMargin(0.15);
  // -- canvas size --
  gStyle->SetCanvasDefW(500);
  gStyle->SetCanvasDefH(500);
  
  std::vector< std::string > bscan2_tags = {
    // magnetic field scans with cover on
    //"4397 Cover Air 0mG           ",              
    //"4398 Cover Air X=+100mG      ",		  
    //"4399 Cover Air X=+50mG       ",		  
    //"4400 Cover Air X=-100mG      ",		  
    //"4401 Cover Air X=-50mG       ",		  
    //"4402 Cover Air Y=+100mG      ",		  
    //"4403 Cover Air Y=+50mG       ",		  
    //"4404 Cover Air Y=-100mG      ",		  
    //"4405 Cover Air Y=-50mG       ",		  
    //"4406 Cover Air Z=+100mG      ",		  
    //"4407 Cover Air Z=+50mG       ",		  
    //"4408 Cover Air Z=-100mG      ",		  
    //"4409 Cover Air Z=-50mG       ",		  
         					  
    // scans with cover off		
    "4427 Bare Air 0mG            ",              
    "4431 Bare Air 0mG            ",		  
    "4433 Bare Air X=+100mG       ",		  
    "4434 Bare Air X=+50mG        ",		  
    "4442 Bare Air Z=+100mG       ",		  
    "4443 Bare Air Z=+50mG        ",		  
    "4444 Bare Air Z=-100mG       ",		  
    "4445 Bare Air Z=-50mG        ",		  
    "4446 Bare Air X=-100mG       ",		  
    "4447 Bare Air X=-50mG        ",		  
    "4448 Bare Air Y=+100mG       ",		  
    "4458 Bare Air Y=-100mG       ",		  
    "4459 Bare Air Y=-50mG        ",		  
    "4493 Bare Air 0mG            ",              
    "4507 Bare Air -70deg X=+100mG",		  
    "4510 Bare Air -70deg, 0mG    ",		  
    "4511 Bare Air X=+100mG       ",		  
     	 				
    // scans in water with cover off	
    "4525 Bare Water 0mG          ",		  
    "4526 Bare Water 0mG          ",
    "4530 Bare Water X=+100mG *   ",
    "4531 Bare Water X=+50mG  *   ",
    "4533 Bare Water X=-100mG     ",
    "4534 Bare Water X=-50mG  *   ",
    "4536 Bare Water Y=+100mG     ",		  
    "4537 Bare Water Y=+50mG      ",		  
    "4538 Bare Water Y=-100mG     ",		  
    "4540 Bare Water Y=-50mG      ",		  
    "4562 Bare Water coils off    ",              
         					  
    // scans in water with cover on	
    "4551 Cover Water coils off   ",
    "4581 Cover Water coils off   " };            

   
  std::vector< TFile* > bscan2_files = {
    //new TFile( "ptf_charge_analysis_run04397.root" ),
    //new TFile( "ptf_charge_analysis_run04398.root" ),
    //new TFile( "ptf_charge_analysis_run04399.root" ),
    //new TFile( "ptf_charge_analysis_run04400.root" ),
    //new TFile( "ptf_charge_analysis_run04401.root" ),
    //new TFile( "ptf_charge_analysis_run04402.root" ),
    //new TFile( "ptf_charge_analysis_run04403.root" ),
    //new TFile( "ptf_charge_analysis_run04404.root" ),
    //new TFile( "ptf_charge_analysis_run04405.root" ),
    //new TFile( "ptf_charge_analysis_run04406.root" ),
    //new TFile( "ptf_charge_analysis_run04407.root" ),
    //new TFile( "ptf_charge_analysis_run04408.root" ),
    //new TFile( "ptf_charge_analysis_run04409.root" ),
    				          	   
                                                     
    new TFile( "ptf_charge_analysis_run04427.root" ),
    new TFile( "ptf_charge_analysis_run04431.root" ),
    new TFile( "ptf_charge_analysis_run04433.root" ),
    new TFile( "ptf_charge_analysis_run04434.root" ),
    new TFile( "ptf_charge_analysis_run04442.root" ),
    new TFile( "ptf_charge_analysis_run04443.root" ),
    new TFile( "ptf_charge_analysis_run04444.root" ),
    new TFile( "ptf_charge_analysis_run04445.root" ),
    new TFile( "ptf_charge_analysis_run04446.root" ),
    new TFile( "ptf_charge_analysis_run04447.root" ),
    new TFile( "ptf_charge_analysis_run04448.root" ),
    new TFile( "ptf_charge_analysis_run04458.root" ),
    new TFile( "ptf_charge_analysis_run04459.root" ),
    new TFile( "ptf_charge_analysis_run04493.root" ),
    new TFile( "ptf_charge_analysis_run04507.root" ),
    new TFile( "ptf_charge_analysis_run04510.root" ),
    new TFile( "ptf_charge_analysis_run04511.root" ),
                                                     
                                                     
    new TFile( "ptf_charge_analysis_run04525.root" ),
    new TFile( "ptf_charge_analysis_run04526.root" ),
    new TFile( "ptf_charge_analysis_run04530.root" ),
    new TFile( "ptf_charge_analysis_run04531.root" ),
    new TFile( "ptf_charge_analysis_run04533.root" ),
    new TFile( "ptf_charge_analysis_run04534.root" ),
    new TFile( "ptf_charge_analysis_run04536.root" ),
    new TFile( "ptf_charge_analysis_run04537.root" ),
    new TFile( "ptf_charge_analysis_run04538.root" ),
    new TFile( "ptf_charge_analysis_run04540.root" ),
    new TFile( "ptf_charge_analysis_run04562.root" ),
    				          	   
                                                     
    new TFile( "ptf_charge_analysis_run04551.root" ), 
    new TFile( "ptf_charge_analysis_run04581.root" )  };
  

  // Make a set of 1D histograms to see how fitted parameters look for each run
  TH1D* h1drun_N      = new TH1D("hrun_N"      ," ; ; N"         , bscan2_tags.size(), 0, bscan2_tags.size() );
  TH1D* h1drun_Q1     = new TH1D("hrun_Q1"     ," ; ; Q_{1}"     , bscan2_tags.size(), 0, bscan2_tags.size() );
  TH1D* h1drun_sigma1 = new TH1D("hrun_sigma1" ," ; ; #sigma_{1}", bscan2_tags.size(), 0, bscan2_tags.size() );
  TH1D* h1drun_mu     = new TH1D("hrun_mu"     ," ; ; #mu"       , bscan2_tags.size(), 0, bscan2_tags.size() );
  TH1D* h1drun_w      = new TH1D("hrun_w"      ," ; ; w"         , bscan2_tags.size(), 0, bscan2_tags.size() );
  TH1D* h1drun_alpha  = new TH1D("hrun_alpha"  ," ; ; #alpha"    , bscan2_tags.size(), 0, bscan2_tags.size() );
  TH1D* hpars[6] = { h1drun_N, h1drun_Q1, h1drun_sigma1, h1drun_mu, h1drun_w, h1drun_alpha };
  for ( unsigned ipar=0; ipar<6; ++ipar ) {
    for (unsigned ibin = 1; ibin <= bscan2_tags.size(); ++ibin ){
      hpars[ipar]->GetXaxis()->SetBinLabel( ibin, bscan2_tags[ ibin-1 ].c_str() );
    }
  }
  
  
  //KEY: TH1D  hqall;  1Charge deposit for all events
  //KEY: TH1D  hqsum;  1Charge deposit for good fit events
  //KEY: TH1D  hqped;  1Charge deposit for pedestal events
  //KEY: TH1D  hphall;  1Pulse height for all events
  //KEY: TH2D  hscanpt;  1Scan point
  //KEY: TH2D  hqavg;  1Average charge
  //KEY: TH2D  hq1pe;  1Q_{1} pe charge
  //KEY: TH2D  hq2pe;  1#sigma_{1} charge width
  //KEY: TH2D  hq21pe;  1#mu average npe
  TLegend * tleg = new TLegend(0.6,0.35,0.95,0.9);
  TCanvas * tc_b2_hqsum = new TCanvas();

  cout<< "\\begin{center}" << endl;
  cout<< "\\begin{tabular}{ l | c c c c c c }" << endl;
  cout << "Run and Description\t&  $q_{1}$\t& $\\sigma_{1}$\t& $\\mu$ \t&  $w$ \t& $\\alpha$ \t& $\\chi/d.o.f.$ \\\\ \\hline" <<endl;

  tc_b2_hqsum->cd();
  tleg->SetFillColor(kWhite);
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){
    TH1D* htmp = (TH1D*)bscan2_files[i]->Get("hqall");
    /*htmp->SetLineWidth(3);
    htmp->SetBit(TF1::kNotDraw);
    htmp->SetLineColor( colors[i] );
    htmp->SetMarkerColor( colors[i] );
    double scale = htmp->GetMaximum();
    htmp->Scale( 1./scale );
    tleg->AddEntry( htmp, bscan2_tags[i].c_str(), "lp" );
    if ( i==0 ) htmp->Draw("hist");
    else htmp->Draw("hist same");
    */
    // print table of fit results from the overall TH1D (Fit was done in previous layer of analysis)
    TF1* ftmp = (TF1*)htmp->GetFunction("pmt_response2");

    cout<< bscan2_tags[i] <<"\t& ";
    //cout<< scientific  << setprecision(4) << setw(12) << ftmp->GetParameter(0) << " $\\pm$ " << ftmp->GetParError(0) <<"\t& ";
    cout<< fixed << setprecision(1) << setw(12) << ftmp->GetParameter(1) << " $\\pm$ " << ftmp->GetParError(1) <<"\t& ";
    cout<< fixed << setprecision(2) << setw(12) << ftmp->GetParameter(2) << " $\\pm$ " << ftmp->GetParError(2) <<"\t& ";
    cout<< fixed << setprecision(4) << setw(12) << ftmp->GetParameter(3) << " $\\pm$ " << ftmp->GetParError(3) <<"\t& ";
    cout<< fixed << setprecision(4) << setw(12) << ftmp->GetParameter(4) <<"\t& ";
    cout<< fixed << setprecision(5) << setw(12) << ftmp->GetParameter(5) <<"\t& ";
    cout<< fixed << setprecision(1) << setw(12) << ftmp->GetChisquare()  << " / " << ftmp->GetNDF();
    cout<< "\\\\ "<<endl;

    h1drun_N      ->SetBinContent( i+1, ftmp->GetParameter(0) );
    h1drun_N      ->SetBinError  ( i+1, ftmp->GetParError (0) );
    h1drun_Q1     ->SetBinContent( i+1, ftmp->GetParameter(1) );
    h1drun_Q1     ->SetBinError  ( i+1, ftmp->GetParError (1) );
    h1drun_sigma1 ->SetBinContent( i+1, ftmp->GetParameter(2) );
    h1drun_sigma1 ->SetBinError  ( i+1, ftmp->GetParError (2) );
    h1drun_mu     ->SetBinContent( i+1, ftmp->GetParameter(3) );
    h1drun_mu     ->SetBinError  ( i+1, ftmp->GetParError (3) );
    h1drun_w      ->SetBinContent( i+1, ftmp->GetParameter(4) );
    h1drun_w      ->SetBinError  ( i+1, ftmp->GetParError (4) );
    h1drun_alpha  ->SetBinContent( i+1, ftmp->GetParameter(5) );
    h1drun_alpha  ->SetBinError  ( i+1, ftmp->GetParError (5) );
    
  }
  tleg->Draw();


  //KEY: TH2D  hq1pe;  1Q_{1} pe charge
  /*
  TCanvas * tc_b2_hq1pe = new TCanvas();
  tc_b2_hq1pe->cd();
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){
    TH2D* htmp2 = (TH2D*)bscan2_files[i]->Get("hq1pe");
    TH1D* htmp = make_1d_from_2d( htmp2, " ; Q_{1} ; scan points / bin", 100, 0., 850. );
    htmp->SetLineWidth(3);
    htmp->SetLineColor( colors[i] );
    htmp->SetMarkerColor( colors[i] );
    double scale = htmp->GetMaximum();
    htmp->Scale( 1./scale );
    if ( i==0 ) htmp->Draw("hist");
    else htmp->Draw("hist same");
    std::cout<<"i="<<i<<" entries="<<htmp->GetEntries()<<std::endl;
  }
  tleg->Draw();
  */


  //KEY: TH2D  hq1pe;  1Q_{1} pe charge
  //KEY: TH2D  hq2pe;  1#sigma_{1} charge width
  //KEY: TH2D  hq21pe;  1#mu average npe
  std::vector< TH1D* > vec_hq1;
  std::vector< TH1D* > vec_hsig;
  std::vector< TH1D* > vec_hmu;
  TH2D* gradhist;
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){
    TH2D* hq1_2d = (TH2D*)bscan2_files[i]->Get("hq1pe");
    TH2D* hsig_2d = (TH2D*)bscan2_files[i]->Get("hq1sig");
    TH2D* hmu_2d = (TH2D*)bscan2_files[i]->Get("hqmupe");

    // Below commented out ... already handled in original code
    //Circle_st pmtedge = find_circle_max_grad( hq1_2d, gradhist );
    //zero_outside_circle( hq1_2d, pmtedge );
    //zero_outside_circle( hsig_2d, pmtedge );
    //zero_outside_circle( hmu_2d, pmtedge );
    
    TH1D* hq1  = make_1d_from_2d( hq1_2d, " ; Q_{1} ; scan points / bin", 100, 0., 850. );
    double scale = hq1->GetMaximum();
    hq1->Scale( 1.0/scale );
    TH1D* hsig = make_1d_from_2d( hsig_2d, " ; #sigma_{1} ; scan points / bin", 100, 0., 500. );
    scale = hsig->GetMaximum();
    hsig->Scale( 1.0/scale );
    TH1D* hmu  = make_1d_from_2d( hmu_2d, " ; #mu ; scan points /bin", 100, 0., 0.5 );
    scale = hmu->GetMaximum();
    hmu->Scale( 1.0/scale );

    vec_hq1.push_back( hq1 );
    vec_hsig.push_back( hsig );
    vec_hmu.push_back( hmu );

    ostringstream os;
    os << "tc_" << i;
    TCanvas * tc = new TCanvas( os.str().c_str(), bscan2_tags[i].c_str(), 1200.0, 800.0 );
    tc->Divide( 3, 2 );
    tc->cd(1);
    hq1_2d->SetMinimum(100.0);
    hq1_2d->SetMaximum(700.0);
    hq1_2d->Draw("colz");
    tc->cd(2);
    hsig_2d->SetMinimum(100.0);
    hsig_2d->SetMaximum(500.0);
    hsig_2d->Draw("colz");
    tc->cd(3);
    hmu_2d->SetMinimum(0.05);
    hmu_2d->SetMaximum(0.5);
    hmu_2d->Draw("colz");
    tc->cd(4);
    hq1->Draw("e1");
    tc->cd(5);
    hsig->Draw("e1");
    tc->cd(6);
    hmu->Draw("e1");

    ostringstream os2;
    os2 << bscan2_files[i]->GetName() << ".pdf";
    tc->Print( os2.str().c_str() );
  }


  // print table of results from the TH1D's
  cout<< "\\begin{center}" << endl;
  cout<< "\\begin{tabular}{ l | c c c }" << endl;
  cout << "Run and Description\t&  $q_{1}$\t& $\sigma_{1}$\t& $\mu$ \\\\ \\hline" <<endl;
  for ( unsigned i=0 ; i < bscan2_files.size() ; ++i ){

    cout<< bscan2_tags[i] <<"\t& ";
    cout<< fixed << setprecision(1) << setw(12 )<< vec_hq1[i]->GetMean() << " $\\pm$ " << vec_hq1[i]->GetMeanError() <<"\t& ";
    cout<< fixed << setprecision(2) << setw(12 )<< vec_hsig[i]->GetMean() << " $\\pm$ " << vec_hsig[i]->GetMeanError()<<"\t& ";
    cout<< fixed << setprecision(4) << setw(12 )<< vec_hmu[i]->GetMean() << " $\\pm$ " << vec_hmu[i]->GetMeanError() ;
    cout<< "\\\\ "<<endl;
  }
  cout<< "\\end{tabular}" <<endl;
  cout<< "\\end{center}" <<endl;



  for (int ipar=0; ipar<6; ++ipar){
    ostringstream os;
    os<<"par_"<<ipar;
    TCanvas* cpar = new TCanvas( os.str().c_str(), os.str().c_str(), 1600, 400 );
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1000000001);
    hpars[ipar]->SetLineWidth(2);
    hpars[ipar]->SetMarkerStyle(24);
    hpars[ipar]->Fit("pol0");
    hpars[ipar]->Draw("e1");
  }
}
