// testing nick's new correlation implementations

#include <cmath>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom3.h"

int main() {
  TH1D* h = new TH1D("test", "geometric correction only;#eta", 100,-1, 1);
  TRandom3 rand;
  for ( int i = 0; i < 10000; ++i ) {
    double triggereta = rand.Rndm()*0.6;
    for ( int j = 0; j < 10000; ++j ) {
      double assoceta = rand.Rndm();
      h->Fill( triggereta - assoceta );
    }
  }
  
  TCanvas c1;
  h->Scale(1/h->GetMaximum() );
  h->Draw();
  c1.SaveAs("testgeom.pdf");
	
	return 0;
}
