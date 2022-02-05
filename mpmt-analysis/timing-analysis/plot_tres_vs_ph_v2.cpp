{

  {


    TGraph *gr = new TGraph();

    gr->SetPoint(0,1,2.24);
    gr->SetPoint(1,2,1.68);
    gr->SetPoint(2,3,1.4);
    gr->SetPoint(3,4,1.3);
    gr->SetPoint(4,6,1.12);
    gr->SetPoint(5,8,1.08);
    gr->SetPoint(6,10,1.02);
    gr->SetPoint(7,15,0.82);
    gr->SetPoint(8,20,0.8);
    gr->SetPoint(9,25,0.7);
    gr->SetPoint(10,30,0.64);

    TCanvas *c = new TCanvas("C","C",0,0,1000,600);
    gr->Draw("AP*");
    gr->SetMarkerStyle(22);
    gr->GetYaxis()->SetRangeUser(0,2.5);
    gr->GetYaxis()->SetTitle("Measured FWHM (ns)");
    gr->GetXaxis()->SetTitle("Pulse Height (PE)");
    gr->SetTitle("PMT Timing resolution vs Pulse Height");
    TF1 *f = new TF1("f1","2.24/sqrt(x)",0.5,35);
    f->Draw("SAME");

    TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
    leg->AddEntry(gr,"data","P");
    leg->AddEntry(f,"2.24/sqrt(PH)");
    leg->Draw("SAME");
    c->SaveAs("tres_vs_ph_injected.png");


    TGraph *gr2 = new TGraph();

    gr2->SetPoint(0,2, 322.8);
    gr2->SetPoint(1,4, 322.8);
    gr2->SetPoint(2,6, 322.9);
    gr2->SetPoint(3,8,322.8);
    gr2->SetPoint(4,10,322.8);
    gr2->SetPoint(5,15,323.2);
    gr2->SetPoint(6,20,323.4);
    gr2->SetPoint(7,25,323.7);
    gr2->SetPoint(8,30,324);
    gr2->SetPoint(9,1,322.6);

    TCanvas *c2 = new TCanvas("C2","C2",0,0,1000,600);
    //gr->Draw("AP*");
    gr2->Draw("AP*");
    gr2->GetYaxis()->SetTitle("Timing distro mean (ns)");
    gr2->GetXaxis()->SetTitle("Pulse Height (PE)");
    gr2->SetTitle("Timing distribution mean vs Pulse Height");
    // gr->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    //gr2->SetMarkerColor(2);

    c2->SaveAs("tres_mean_vs_ph.png");


  }

}
