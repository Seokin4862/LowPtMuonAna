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


TH1 *h = (TH1*)f->Get( "mu_sigma_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "mu_sigma_vs_p" );

h->Draw(  );
h2->SetLineColor( kRed );
h2->Draw( "same" );
c->Print( "mu_sigma_vs_p_superimposed.png" );


TH1 *h = (TH1*)f->Get( "mu_yield_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "mu_yield_vs_p" );

h->Draw(  );
h2->SetLineColor( kRed );
h2->Draw( "same" );
c->Print( "mu_yield_vs_p_superimposed.png" );

//pions

TH1 *h = (TH1*)f->Get( "pi_lambda_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "pi_lambda_vs_p" );

h->Draw(  );
h2->SetLineColor( kRed );
h2->Draw( "same" );
c->Print( "pi_lambda_vs_p_superimposed.png" );


TH1 *h = (TH1*)f->Get( "pi_sigma_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "pi_sigma_vs_p" );

h->Draw(  );
h2->SetLineColor( kRed );
h2->Draw( "same" );
c->Print( "pi_lambda_vs_p_superimposed.png" );


TH1 *h = (TH1*)f->Get( "pi_yield_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "pi_yield_vs_p" );

h->Draw(  );
h2->SetLineColor( kRed );
h2->Draw( "same" );
c->Print( "pi_lambda_vs_p_superimposed.png" );

}
