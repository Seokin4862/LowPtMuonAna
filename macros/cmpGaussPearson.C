void cmpGaussPearson() {

TFile * f = new TFile( "/home/sy34/workspace/LowPtMuonAna/bin/out/fits_gauss.root" );
TFile * f2 = new TFile( "/home/sy34/workspace/LowPtMuonAna/bin/out/fits_pearson.root" );


TCanvas * c = new TCanvas("c","c");

TH1 *h = (TH1*)f->Get( "mu_lambda_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "mu_lambda_vs_p" );

h->Draw(  );
h2->SetLineColor( kRed );
h2->Draw( "same" );
c->Print( "mu_lambda_vs_p_superimposed.png" );

}
