{

  {


    TGraph *gr = new TGraph();

    gr->SetPoint(0,1,1.53/sqrt(2));
    gr->SetPoint(1,2,1.25/sqrt(2));
    gr->SetPoint(2,3,1.18/sqrt(2));
    gr->SetPoint(3,4,1.04/sqrt(2));
    gr->SetPoint(4,8,0.93/sqrt(2));
    gr->SetPoint(5,12,0.84/sqrt(2));
    gr->SetPoint(6,20,0.64/sqrt(2));

    TCanvas *c = new TCanvas("C","C",0,0,1000,600);
    gr->Draw("AP*");
    gr->SetMarkerStyle(22);
    gr->GetYaxis()->SetRangeUser(0,1.2);
    gr->GetYaxis()->SetTitle("Measured Sigma (ns)");
    gr->GetXaxis()->SetTitle("Pulse Height (PE)");
    gr->SetTitle("PMT Timing resolution vs Pulse Height");
    TF1 *f = new TF1("f1","1.08/sqrt(x)",0.5,24);
    f->Draw("SAME");

    TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
    leg->AddEntry(gr,"data","P");
    leg->AddEntry(f,"1.08/sqrt(PH)");
    leg->Draw("SAME");
    c->SaveAs("tres_vs_ph.png");


    TGraph *gr2 = new TGraph();

    gr2->SetPoint(0,1, 1.44/sqrt(2));
    gr2->SetPoint(1,2, 1.19/sqrt(2));
    gr2->SetPoint(2,3, 1.04/sqrt(2));
    gr2->SetPoint(3,4, 0.94/sqrt(2));
    gr2->SetPoint(4,8, 0.78/sqrt(2));
    gr2->SetPoint(5,12,0.74/sqrt(2));
    gr2->SetPoint(6,20,0.62/sqrt(2));

    TCanvas *c2 = new TCanvas("C2","C2",0,0,1000,600);
    gr->Draw("AP*");
    gr2->Draw("*");
    gr->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(2);

    TLegend *leg2 = new TLegend(0.6,0.7,0.89,0.89);
    leg2->AddEntry(gr,"1E7 gain","P");
    leg2->AddEntry(gr2,"2E7 gain","p");
    leg2->Draw("SAME");
    c2->SaveAs("tres_vs_ph_hv.png");


  }

}
