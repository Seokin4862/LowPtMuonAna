void cmpGaussPearson() {

TFile * f = new TFile( "" );
TFile * f2 = new TFile( "" );


TCanvas * c = new TCanvas("c","c");

TH1 *h = (TH1*)f->Get( "mu_lambda_vs_p" );
TH1 *h2 = (TH1*)f2->Get( "mu_lambda_vs_p" );

h->Draw(  );
h2->Draw( "same" );
c->Print( "mu_lambda_vs_p_superimposed" );

}
