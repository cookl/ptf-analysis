
/// PlotElements
/// Simple structure to hold the element that make up a plot
/// In particular the TGraphErrors, TCanvas, and TF1
struct PlotElements {
  TGraphErrors* tg = nullptr;
  TCanvas* tc = nullptr;
  TF1* f1 = nullptr;
};


PlotElements* make_plot( string name, string title, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xe, const std::vector<double>& ye, int color=kBlack  ){
  PlotElements* pe = new PlotElements();
  pe->tg = new TGraphErrors( x.size(), &x[0], &y[0], &xe[0], &ye[0] );
  pe->tg->SetTitle( title.c_str() );

  pe->tc = new TCanvas(name.c_str(),name.c_str(), 500, 500 );
  pe->tc->cd();
  pe->tg->SetMarkerStyle(33);
  pe->tg->SetLineWidth(2);
  pe->tg->SetMarkerColor( color );
  pe->tg->SetLineColor( color );
  string fname = string("f1_")+name;				    
  pe->f1 = new TF1(fname.c_str(),"pol1");
  if (name=="pe_s1_bwy"){
    pe->f1->SetParameters(250.0, -0.09 );
  }
  
  if (name=="pe_q1_bax"){
    std::cout<<"pe_q1_bax"<<std::endl;
    pe->f1->SetParameters(380, -0.17 );
  }
  pe->f1->SetLineColor(color+2);
  pe->f1->SetLineWidth(3);
  pe->tg->Fit( pe->f1 );
  pe->tg->Draw("AP");
  return pe;
}

void waveform_fit_results2(){
    
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Plots of response when Bare in Air as Magnetic field in X is changed
 vector< double > x_ba   = { -100., -50., 0., 0., 0., 50., 100., 100. };
  vector< double > xe_ba  = { 10. ,  10., 10.  , 10., 10., 10., 10., 10. };
  vector< double > q1_ba  = {399.9, 392.2, 388.3, 387.8, 383.0, 379.2, 364.0, 364.8};
  vector< double > q1e_ba = { 0.3, 0.3, 0.3, 0.2, 0.3, 0.3, 0.4, 0.2 };
  vector< double > s1_ba  = { 254.30, 246.32, 247.47, 245.58, 243.28, 245.10, 250.54, 251.64};
  vector< double > s1e_ba = { 0.36, 0.35, 0.36, 0.18, 0.34, 0.36, 0.39, 0.19 };
  vector< double > chi2_ba = {   4336.9 / 47, 3963.5 / 47, 3609.9 / 47, 14023.8 / 47, 3820.6 / 47, 3769.3 / 47, 4002.7 / 47, 20069.8 / 47 };

  // adjust errors to account for rather poor model fit chi2/ndof?
  bool inflate_errors = true;
  if (inflate_errors) {
    for ( unsigned i=0; i< q1e_ba.size(); ++i ){
      double inflation = sqrt( chi2_ba[i] );
      q1e_ba[i] *= inflation;
      s1e_ba[i] *= inflation;
    }
  }

  cout<<"Check sizes : "<<x_ba.size()<<", "<<q1_ba.size()<<", "<<xe_ba.size()<<", "<<q1e_ba.size()<<std::endl;
  
  PlotElements * pe_q1_bax = make_plot( "pe_q1_bax", "SK PMT in Bare Air; Bx (Gauss); Q_{1}"     , x_ba, q1_ba, xe_ba, q1e_ba, kBlue );
  PlotElements * pe_s1_bax = make_plot( "pe_s1_bax", "SK PMT in Bare Air; Bx (Gauss); #sigma_{1}", x_ba, s1_ba, xe_ba, s1e_ba, kBlue );



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Plots of response when Bare in Air as Magnetic field in Z is changed
  vector< double > z_ba   = { -100, -50 , 0   , 0   , 0   , +50 , +100 };
  vector< double > ze_ba  = {   40.,   40., 40.  , 40.  , 40.  ,  40. ,    40.};
  vector< double > q1_baz  = { 382.1, 381.7, 388.3, 387.8, 383.0, 383.3, 386.3 };
  vector< double > q1e_baz = { 0.4, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3 };
  vector< double > s1_baz  = { 249.67, 246.57, 247.47, 245.58, 243.28, 246.17, 248.16 };
  vector< double > s1e_baz = { 0.36, 0.35, 0.36, 0.18, 0.34, 0.35, 0.35 };
  vector< double > chi2_baz = {   4402.0 / 47, 3925.0 / 47, 3609.9 / 47, 14023.8 / 47, 3820.6 / 47, 3843.7 / 47, 3845.3 / 47 };
    
  // adjust errors to account for rather poor model fit chi2/ndof?
  if (inflate_errors) {
    for ( unsigned i=0; i< q1e_baz.size(); ++i ){
      double inflation = sqrt( chi2_baz[i] );
      q1e_baz[i] *= inflation;
      s1e_baz[i] *= inflation;
    }
  }

  PlotElements * pe_q1_baz = make_plot( "pe_q1_baz", "SK PMT in Bare Air; Bz (Gauss); Q_{1}"     , z_ba, q1_baz, ze_ba, q1e_baz, kBlue );
  PlotElements * pe_s1_baz = make_plot( "pe_s1_baz", "SK PMT in Bare Air; Bz (Gauss); #sigma_{1}", z_ba, s1_baz, ze_ba, s1e_baz, kBlue );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Plots of response when Bare in Air as Magnetic field in Y is changed
  vector< double > y_ba   = { -100, -50   , 0   , 0   , 0 , +100 };
  vector< double > ye_ba  = {   10., 10.  , 10.  , 10.  ,  10. ,    10.};
  vector< double > q1_bay = { 388.5, 384.9, 388.3, 387.8, 383.0, 381.5 };
  vector< double > q1e_bay = {  0.4, 0.3, 0.3, 0.2, 0.3, 0.3 };
  vector< double > s1_bay  = { 261.04, 251.98, 247.47, 245.58, 243.28, 241.60 };
  vector< double > s1e_bay = { 0.38, 0.36, 0.36, 0.18, 0.34, 0.34 };
  vector< double > chi2_bay = {   4500.7 / 47, 4372.9 / 47, 3609.9 / 47, 14023.8 / 47, 3820.6 / 47, 3240.9 / 47 };

  // adjust errors to account for rather poor model fit chi2/ndof?
  if (inflate_errors) {
    for ( unsigned i=0; i< q1e_bay.size(); ++i ){
      double inflation = sqrt( chi2_bay[i] );
      q1e_bay[i] *= inflation;
      s1e_bay[i] *= inflation;
    }
  }

  cout<<"Check sizes : "<<y_ba.size()<<", "<<q1_bay.size()<<", "<<ye_ba.size()<<", "<<q1e_bay.size()<<std::endl;

  PlotElements * pe_q1_bay = make_plot( "pe_q1_bay", "SK PMT in Bare Air; By (Gauss); Q_{1}"     , y_ba, q1_bay, ye_ba, q1e_bay, kBlue );
  PlotElements * pe_s1_bay = make_plot( "pe_s1_bay", "SK PMT in Bare Air; By (Gauss); #sigma_{1}", y_ba, s1_bay, ye_ba, s1e_bay, kBlue );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Plots of response when Bare in Water as Magnetic field in Y is changed
  vector< double > y_bw   = { -100, -50., 0.  , 0   , 50 , 100 };
  vector< double > ye_bw  = {   10., 10.  , 10.  , 10.  ,  10. ,    10.};
  vector< double > q1_bwy  = { 392.4, 384.2, 382.9, 378.1, 383.2, 381.5 };
  vector< double > q1e_bwy = { 0.2, 0.1, 0.1, 0.1, 0.1, 0.1 };
  vector< double > s1_bwy  = { 252.44, 241.49, 236.78, 235.12, 233.73, 232.49 };
  vector< double > s1e_bwy = { 0.17, 0.16, 0.15, 0.16, 0.15, 0.15 };
  vector< double > chi2_bwy = { 19942.1 / 47, 16469.8 / 47, 14871.5 / 47, 14314.8 / 47, 12696.0 / 47, 12316.6 / 47 };
    
  // adjust errors to account for rather poor model fit chi2/ndof?
  if (inflate_errors) {
    for ( unsigned i=0; i< q1e_bwy.size(); ++i ){
      double inflation = sqrt( chi2_bwy[i] );
      q1e_bwy[i] *= inflation;
      s1e_bwy[i] *= inflation;
    }
  }

  PlotElements * pe_q1_bwy = make_plot( "pe_q1_bwy", "SK PMT Bare in Water; By (Gauss); Q_{1}"     , y_bw, q1_bwy, ye_bw, q1e_bwy, kRed );

  cout<<"Check sizes : "<<y_bw.size()<<", "<<s1_bwy.size()<<", "<<ye_bw.size()<<", "<<s1e_bwy.size()<<std::endl;
  
  PlotElements * pe_s1_bwy = make_plot( "pe_s1_bwy", "SK PMT Bare in Water; By (Gauss); #sigma_{1}", y_bw, s1_bwy, ye_bw, s1e_bwy, kRed );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Plots of response when Bare in Water as Magnetic field in X is changed
  vector< double > x_bw   = {-100., -50.,  0.,  0., 50., 100. }; //{  -100., 0, 0 };
  vector< double > xe_bw  = {  10.,  10., 10., 10., 10.,  10. };
  vector< double > q1_bwx  = { 395.9, 392.6, 382.9, 378.1, 368.7, 361.8 };
  vector< double > q1e_bwx = {   0.2,   0.1, 0.1,   0.1,   0.2,   0.2 };
  vector< double > s1_bwx  = { 249.55, 239.00, 236.78, 235.12, 242.02, 245.35 };    
  vector< double > s1e_bwx = { 0.17, 0.16, 0.15, 0.16, 0.17, 0.18 };
  vector< double > chi2_bwx = { 17281.5 / 47,
				15187.4 / 47,
				14871.5 / 47,
				14314.8 / 47,
				21116.1 / 47,
    				31036.2 / 47 };
    
    
  // adjust errors to account for rather poor model fit chi2/ndof?
  if (inflate_errors) {
    for ( unsigned i=0; i< q1e_bwx.size(); ++i ){
      double inflation = sqrt( chi2_bwx[i] );
      q1e_bwx[i] *= inflation;
      s1e_bwx[i] *= inflation;
    }
  }

  PlotElements * pe_q1_bwx = make_plot( "pe_q1_bwx", "SK PMT Bare in Water; Bx (Gauss); Q_{1}"     , x_bw, q1_bwx, xe_bw, q1e_bwx, kRed );
  PlotElements * pe_s1_bwx = make_plot( "pe_s1_bwx", "SK PMT Bare in Water; Bx (Gauss); #sigma_{1}", x_bw, s1_bwx, xe_bw, s1e_bwx, kRed );



    // Overlay "Bare in Water" onto "Bare in Air"
  TLegend * tleg = new TLegend(0.3,0.6,0.7,0.9);
  tleg->SetFillColor(kWhite);
  tleg->AddEntry( pe_q1_bay->tg, "Air", "lp");
  tleg->AddEntry( pe_q1_bwy->tg, "Water", "lp");

  pe_q1_bay->tc->cd();
  pe_q1_bwy->tg->Draw("P");
  tleg->Draw();

  pe_s1_bay->tc->cd();
  pe_s1_bwy->tg->Draw("P");
  tleg->Draw();


  /// and in X  
  pe_q1_bax->tc->cd();
  pe_q1_bwx->tg->Draw("P");
  tleg->Draw();

  pe_s1_bax->tc->cd();
  pe_s1_bwx->tg->Draw("P");
  tleg->Draw();




}
