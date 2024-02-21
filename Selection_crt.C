#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

{
  gROOT->Reset();
  const TString SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";

  //  gStyle->SetOptStat(0); //Removing the stats box
  gStyle->SetPalette(kCandy);
  TColor::InvertPalette(); 
  TFile  myFile("/sbnd/data/users/arellano/tpc_comm/production_231127_brazilCM/sbnd_crt_lifetime/trees/muon_hitdumper_NS_elifetime15ms_sce.root");
  TTree* myTree = (TTree*) myFile.Get("hitdumper/hitdumpertree");
  vector<double>*  ct_x1;
  vector<double>*  ct_x2;
  vector<double>*  ct_y1;
  vector<double>*  ct_y2;
  vector<double>*  ct_z1;
  vector<double>*  ct_z2;
  vector<double>*  chit_x;
  vector<double>*  chit_y;
  vector<double>*  chit_z;
  vector<double>* muontrk_type;
  vector<double>*muontrk_theta_xz;
  vector<double>*muontrk_theta_yz;
  vector<double>*muontrk_tpc;
  double dx;
  double dy;
  double dz;
  double CRT_theta_xz;
  double CRT_theta_yz;
  double modified_CRT_theta_yz;
  myTree->SetBranchAddress("ct_x1", &ct_x1);
  myTree->SetBranchAddress("ct_x2", &ct_x2);
  myTree->SetBranchAddress("ct_y1", &ct_y1);
  myTree->SetBranchAddress("ct_y2", &ct_y2);
  myTree->SetBranchAddress("ct_z1", &ct_z1);
  myTree->SetBranchAddress("ct_z2", &ct_z2); 
  myTree->SetBranchAddress("chit_x", &chit_x);
  myTree->SetBranchAddress("chit_y", &chit_y);
  myTree->SetBranchAddress("chit_z", &chit_z); 
  myTree->SetBranchAddress("muontrk_type", &muontrk_type);
  myTree->SetBranchAddress("muontrk_theta_xz", &muontrk_theta_xz);
  myTree->SetBranchAddress("muontrk_theta_yz", &muontrk_theta_yz);
  myTree->SetBranchAddress("muontrk_tpc", &muontrk_tpc);
  gROOT->cd(0);
  TH2F *x_y_Hit_1_South = new TH2F("xy1","xy1",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_1_North = new TH2F("xy2","xy2",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_2_South_TPC1 = new TH2F("xy3","xy3",70, -360, 0, 70, -360, 0);
  TH2F *x_y_Hit_2_North_TPC1 = new TH2F("xy4","xy4",70, -360, 0, 70, -360, 0);
  TH2F *x_y_Hit_2_South_TPC2 = new TH2F("xy7","xy7",70, 0, 360, 70, 0, 360);
  TH2F *x_y_Hit_2_North_TPC2 = new TH2F("xy8","xy8",70, 0, 360, 70, 0, 360);
  TH2F *x_z_Hit_1_Bottom = new TH2F("xz1","xz1",70, -420, 420, 70, -200, 800);
  TH2F *x_z_Hit_2_Bottom = new TH2F("xz2","xz2",70, -420, 420, 70, -200, 800);
  TH2F *x_y_AllHits_South = new TH2F("xy5","xy5",70, -360, 360, 70, -360, 360);
  TH2F *x_y_AllHits_North = new TH2F("xy6","xy6",70, -360, 360, 70, -360, 360);
  TH2F *x_z_AllHits_Bottom = new TH2F("xz3","xz3",70, -420, 420, 70, -200, 800);
  TH1F *Theta_xz_South =new TH1F("theta_xz_1", "theta_xz_1", 50, -50, 50);
  TH1F *Theta_yz_South =new TH1F("theta_yz_1", "theta_yz_1", 50, -50,50);
  TH1F *Theta_xz_North =new TH1F("theta_xz_2", "theta_xz_2", 50, -50,50);
  TH1F *Theta_yz_North =new TH1F("theta_yz_2", "theta_yz_2", 50, -50,50);
  TH1F *CRT_Theta_xz_South =new TH1F("crttheta_xz_1", "crttheta_xz_1", 50, -200, 200);
  TH1F *CRT_Theta_yz_South =new TH1F("crttheta_yz_1", "crttheta_yz_1", 50, -200,200);
  TH1F *CRT_Theta_xz_North =new TH1F("crttheta_xz_2", "crttheta_xz_2", 50, -200,200);
  TH1F *CRT_Theta_yz_North =new TH1F("crttheta_yz_2", "crttheta_yz_2", 50, -200,200);
  int nEntries = myTree->GetEntries(); // Get the number of entries in this tree
  for (int iEnt = 0; iEnt < nEntries; iEnt++) { //Big loop of entries
    myTree->GetEntry(iEnt); // Gets the next entry (filling the linked variables)
    for (size_t i=0; i < ct_x1->size(); ++i) {
      dx =ct_x2->at(i) - ct_x1->at(i);
      dy =ct_y2->at(i) - ct_y1->at(i);
      dz =ct_z2->at(i) - ct_z1->at(i);
      CRT_theta_xz= atan2(dx,dz)*(180/TMath::Pi());
      CRT_theta_yz= atan2(dy,dz)*(180/TMath::Pi());
      modified_CRT_theta_yz=CRT_theta_yz;
      CRT_Theta_xz_South->Fill(CRT_theta_xz);                                                                                                                                                               
      CRT_Theta_yz_South->Fill(CRT_theta_yz);
      if(dy<0.0) modified_CRT_theta_yz=-CRT_theta_yz;
      if(ct_z1->at(i)>-200 && ct_z1->at(i)<-150){
	if((ct_x1->at(i)>-300 && ct_x1->at(i)<-250)||(ct_x1->at(i)<300 && ct_x1->at(i)>250)){
	  x_y_Hit_1_South->Fill(ct_x1->at(i),ct_y1->at(i));
	}
      }
      if(ct_z2->at(i)>-200 && ct_z2->at(i)<-150 && ct_z1->at(i)>750 && ct_z1->at(i)<800 ){   
	if(ct_x2->at(i)>-300 && ct_x2->at(i)<-250){
	  if(ct_x1->at(i)>-300 && ct_x1->at(i)<-250){
	    x_y_Hit_2_South_TPC1->Fill(ct_x2->at(i),ct_y2->at(i));
	  }
	}
      }

      if(ct_z2->at(i)>-200 && ct_z2->at(i)<-150 && ct_z1->at(i)>750 && ct_z1->at(i)<800 ){
        if(ct_x2->at(i)>250 && ct_x2->at(i)<300){
          if(ct_x1->at(i)>250 && ct_x1->at(i)<300){
	    x_y_Hit_2_South_TPC2->Fill(ct_x2->at(i),ct_y2->at(i));
          }
        }
      }
      if(ct_z1->at(i)>750 && ct_z1->at(i)<800){
	if((ct_x1->at(i)>-300 && ct_x1->at(i)<-250)||(ct_x1->at(i)<300 && ct_x1->at(i)>250)){
	  x_y_Hit_1_North->Fill(ct_x1->at(i),ct_y1->at(i));
	}
      }
      if(ct_z2->at(i)>750 && ct_z2->at(i)<800 && ct_z1->at(i)>-200 && ct_z1->at(i)<-150){
	if(ct_x2->at(i)>-300 && ct_x2->at(i)<-250){
	  if(ct_x1->at(i)>-300 && ct_x1->at(i)<-250){
	    x_y_Hit_2_North_TPC1->Fill(ct_x2->at(i),ct_y2->at(i));
	  }
	}
      }   
      if(ct_z2->at(i)>750 && ct_z2->at(i)<800 && ct_z1->at(i)>-200 && ct_z1->at(i)<-150){
        if(ct_x2->at(i)>270 && ct_x2->at(i)<300){
          if(ct_x1->at(i)>270 && ct_x1->at(i)<300){
	    x_y_Hit_2_North_TPC2->Fill(ct_x2->at(i),ct_y2->at(i));
          }
        }
      }
      if(ct_y1->at(i)>-400 && ct_y1->at(i)<-350){
	x_z_Hit_1_Bottom->Fill(ct_x1->at(i),ct_z1->at(i));
      }
      if(ct_y2->at(i)>-400 && ct_y2->at(i)<-350){
	x_z_Hit_2_Bottom->Fill(ct_x2->at(i),ct_z2->at(i));
      }
    }


    //*************** Selection using chit variable start ***************// 

    // TPC loop 1 (used to define cuts)
    bool HasTrkType4= false;
    for(size_t j=0; j<muontrk_type->size(); ++j){ //Looping over TPC tracks
      if(muontrk_type->at(j)==4){//Selecting track type 4
	HasTrkType4=true;
      }// if
    }//for


    //CRT loop 1 (used to define cuts)
    bool IsInHitRegionSouth=false;
    bool IsInHitRegionNorth=false;
    for(size_t i=0; i < chit_x->size(); ++i){ //Looping over CRT hits
      if(chit_z->at(i)>-200 && chit_z->at(i)<-150 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360 ){ //Selectiong South CRT
	if((chit_x->at(i)>-300 && chit_x->at(i)<-250)||(chit_x->at(i)<300 && chit_x->at(i)>250)){ // Cut on x region
	  IsInHitRegionSouth=true;
	} //CRT if
      }//x cut if
      if(chit_z->at(i)>750 && chit_z->at(i)<800 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360){ //Selecting north CRT
	if((chit_x->at(i)>-300 && chit_x->at(i)<-250)||(chit_x->at(i)<300 && chit_x->at(i)>250)){ //Cut on x region
	  IsInHitRegionNorth=true;
	}// CRT if
      }// x cut if
    }//CRT hits for
     
    if(HasTrkType4 && (IsInHitRegionSouth || IsInHitRegionNorth)){  // Cuts required
	//TPC loop 2 (used to fill histograms)
	for(size_t j=0; j<muontrk_type->size(); ++j){ //Looping over TPC tracks 

	  //modifying TPC angles start     
	  double modified_theta_xz = muontrk_theta_xz->at(j);
	  double  modified_theta_yz =muontrk_theta_yz->at(j);

	  if ( muontrk_theta_xz->at(j) > 90 && muontrk_tpc->at(j) == 0 ) {
	    modified_theta_xz = muontrk_theta_xz->at(j) - 180.0;
	  }
	  else if ( muontrk_theta_xz->at(j) > 90 && muontrk_tpc->at(j) == 1 ) {
	    modified_theta_xz = 180.0 - muontrk_theta_xz->at(j);
	  }
	  else if ( muontrk_theta_xz->at(j) < 90 && muontrk_tpc->at(j) == 1 ) {
	    modified_theta_xz = -muontrk_theta_xz->at(j);
	  }

	  if ( muontrk_theta_yz->at(j) < -90 ) {
	    modified_theta_yz = muontrk_theta_yz->at(j) + 180.;
	  }
	  if ( muontrk_theta_yz->at(j) > 90 ) {
	    modified_theta_yz = muontrk_theta_yz->at(j) - 180.;
	  }
	  //modifying TPC angles end
	 
	  if(muontrk_type->at(j)==4){
	    Theta_xz_South->Fill(modified_theta_xz);
	    Theta_yz_South->Fill(modified_theta_yz);
	  }
	}//TPC trk loop end

	//CRT loop 2 (Used to fill histograms)
	for(size_t i=0; i < chit_x->size(); ++i){ //Looping over CRT hits
	  if(chit_z->at(i)>-200 && chit_z->at(i)<-150 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360 ){ //Selectiong South CRT 
	    if((chit_x->at(i)>-300 && chit_x->at(i)<-250)||(chit_x->at(i)<300 && chit_x->at(i)>250)){ //Cut on x region
	    x_y_AllHits_South->Fill(chit_x->at(i),chit_y->at(i));
	    } //South CRT
	  }// x region
	  if(chit_z->at(i)>750 && chit_z->at(i)<800 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360){ //Selecting north CRT  
	    if((chit_x->at(i)>-300 && chit_x->at(i)<-250)||(chit_x->at(i)<300 && chit_x->at(i)>250)){ //Cut on x region 
	    x_y_AllHits_North->Fill(chit_x->at(i),chit_y->at(i));
	    }// North CRT
	  }// x region
	}//CRT hit loop end
      }// end of required cuts if 

   //*************** Selection using chit variable end ***************//


  } //Big loop of entries
 
    //*************** Making histograms pretty ********************//
  /*
    TCanvas*  x_y_1 = new TCanvas ("xy1", "xy1", 900, 700);

    x_y_Hit_1_South->Draw("COLZ");
    x_y_Hit_1_South->SetTitle("South CRT Muon Hit 1");
    x_y_Hit_1_South->GetXaxis()->SetTitle("x Position (cm)");
    x_y_Hit_1_South->GetYaxis()->SetTitle("y Position (cm)");

    x_y_1->SaveAs(SaveDir + "/x_y_Hit_1_South.pdf");  

    TCanvas*  x_y_2 = new TCanvas ("xy2", "xy2", 900, 700);
    x_y_Hit_1_North->Draw("COLZ");
    x_y_Hit_1_North->SetTitle("North CRT Muon Hit 1");
    x_y_Hit_1_North->GetXaxis()->SetTitle("x Position (cm)");
    x_y_Hit_1_North->GetYaxis()->SetTitle("y Position (cm)");

    x_y_2->SaveAs(SaveDir + "/x_y_Hit_1_North.pdf");

    TCanvas*  x_y_3 = new TCanvas ("xy3", "xy3", 900, 700);

    x_y_Hit_2_South_TPC1->Draw("COLZ");
    x_y_Hit_2_South_TPC1->SetTitle("South CRT Muon Hit 2 TPC 1");
    x_y_Hit_2_South_TPC1->GetXaxis()->SetTitle("x Position (cm)");
    x_y_Hit_2_South_TPC1->GetYaxis()->SetTitle("y Position (cm)");

    x_y_3->SaveAs(SaveDir + "/x_y_Hit_2_South_TPC1.pdf");

    TCanvas*  x_y_4 = new TCanvas ("xy4", "xy4", 900, 700);
    x_y_Hit_2_North_TPC1->Draw("COLZ");
    x_y_Hit_2_North_TPC1->SetTitle("North CRT Muon Hit 2 TPC 1");
    x_y_Hit_2_North_TPC1->GetXaxis()->SetTitle("x Position (cm)");
    x_y_Hit_2_North_TPC1->GetYaxis()->SetTitle("y Position (cm)");

    x_y_4->SaveAs(SaveDir + "/x_y_Hit_2_North_TPC1.pdf");

    TCanvas*  x_y_7 = new TCanvas ("xy7", "xy7", 900, 700);

    x_y_Hit_2_South_TPC2->Draw("COLZ");
    x_y_Hit_2_South_TPC2->SetTitle("South CRT Muon Hit 2 TPC 2");
    x_y_Hit_2_South_TPC2->GetXaxis()->SetTitle("x Position (cm)");
    x_y_Hit_2_South_TPC2->GetYaxis()->SetTitle("y Position (cm)");

    x_y_7->SaveAs(SaveDir + "/x_y_Hit_2_South_TPC2.pdf");

    TCanvas*  x_y_8 = new TCanvas ("xy8", "xy8", 900, 700);
    x_y_Hit_2_North_TPC2->Draw("COLZ");
    x_y_Hit_2_North_TPC2->SetTitle("North CRT Muon Hit 2 TPC 2");
    x_y_Hit_2_North_TPC2->GetXaxis()->SetTitle("x Position (cm)");
    x_y_Hit_2_North_TPC2->GetYaxis()->SetTitle("y Position (cm)");

    x_y_8->SaveAs(SaveDir + "/x_y_Hit_2_North_TPC2.pdf");



    TCanvas*  x_z_1 = new TCanvas ("xz1", "xz1", 900, 700);
    x_z_Hit_1_Bottom->Draw("COLZ");
    x_z_Hit_1_Bottom->SetTitle("Bottom CRT Muon Hit 1");
    x_z_Hit_1_Bottom->GetXaxis()->SetTitle("x Position (cm)");
    x_z_Hit_1_Bottom->GetYaxis()->SetTitle("z Position (cm)");

    x_z_1->SaveAs(SaveDir + "/x_z_Hit_1_Bottom.pdf");

    TCanvas*  x_z_2 = new TCanvas ("xz2", "xz2", 900, 700);

    x_z_Hit_2_Bottom->Draw("COLZ");
    x_z_Hit_2_Bottom->SetTitle("Bottom CRT Muon Hit 2");
    x_z_Hit_2_Bottom->GetXaxis()->SetTitle("x Position (cm)");
    x_z_Hit_2_Bottom->GetYaxis()->SetTitle("z Position (cm)");

    x_z_2->SaveAs(SaveDir + "/x_z_Hit_2_Bottom.pdf");
  */ 

  
    TCanvas*  x_y_5 = new TCanvas ("xy5", "xy5", 900, 700);

    x_y_AllHits_South->Draw("COLZ");
    x_y_AllHits_South->SetTitle("South CRT Muon Hits");
    x_y_AllHits_South->GetXaxis()->SetTitle("x Position (cm)");
    x_y_AllHits_South->GetYaxis()->SetTitle("y Position (cm)");

    x_y_5->SaveAs(SaveDir + "/x_y_AllHits_South.pdf");

    TCanvas*  x_y_6 = new TCanvas ("xy6", "xy6", 900, 700);

    x_y_AllHits_North->Draw("COLZ");
    x_y_AllHits_North->SetTitle("North CRT Muon Hits");
    x_y_AllHits_North->GetXaxis()->SetTitle("x Position (cm)");
    x_y_AllHits_North->GetYaxis()->SetTitle("y Position (cm)");

    x_y_6->SaveAs(SaveDir + "/x_y_AllHits_North.pdf");
    /*
    TCanvas*  x_z_3 = new TCanvas ("xz3", "xz3", 900, 700);

    x_z_AllHits_Bottom->Draw("COLZ");
    x_z_AllHits_Bottom->SetTitle("Bottom CRT Muon Hits");
    x_z_AllHits_Bottom->GetXaxis()->SetTitle("x Position (cm)");
    x_z_AllHits_Bottom->GetYaxis()->SetTitle("z Position (cm)");

    x_z_3->SaveAs(SaveDir + "/x_z_AllHits_Bottom.pdf");
  */
  TCanvas*  theta_xz_South = new TCanvas ("theta_xz1", "theta_xz1", 900, 700);

  Theta_xz_South->Draw("");
  Theta_xz_South->SetTitle("#theta_{xz} South");
  Theta_xz_South->GetXaxis()->SetTitle("Angle (degrees)");
  Theta_xz_South->GetYaxis()->SetTitle("");
  Theta_xz_South->SetLineWidth(4);
  Theta_xz_South->SetLineColor(kRed-6);

  theta_xz_South->SaveAs(SaveDir + "/Theta_xz_south.pdf");
  /*
  TCanvas*  theta_xz_North = new TCanvas ("theta_xz2", "theta_xz2", 900, 700);

  Theta_xz_North->Draw("");
  Theta_xz_North->SetTitle("#theta_{xz} North");
  Theta_xz_North->GetXaxis()->SetTitle("Angle (degrees)");
  Theta_xz_North->GetYaxis()->SetTitle("");
  Theta_xz_North->SetLineWidth(4);
  Theta_xz_North->SetLineColor(kRed-6);

  theta_xz_North->SaveAs(SaveDir + "/Theta_xz_North.pdf");
  */
  TCanvas*  theta_yz_South = new TCanvas ("theta_xy1", "theta_xy1", 900, 700);

  Theta_yz_South->Draw("");
  Theta_yz_South->SetTitle("#theta_{yz} South");
  Theta_yz_South->GetXaxis()->SetTitle("Angle (degrees)");
  Theta_yz_South->GetYaxis()->SetTitle("");
  Theta_yz_South->SetLineWidth(4);
  Theta_yz_South->SetLineColor(kRed-6);

  theta_yz_South->SaveAs(SaveDir + "/Theta_yz_south.pdf");
  /*
  TCanvas*  theta_yz_North = new TCanvas ("theta_xy2", "theta_xy2", 900, 700);

  Theta_yz_North->Draw("");
  Theta_yz_North->SetTitle("#theta_{yz} North");
  Theta_yz_North->GetXaxis()->SetTitle("Angle (degrees)");
  Theta_yz_North->GetYaxis()->SetTitle("");
  Theta_yz_North->SetLineWidth(4);
  Theta_yz_North->SetLineColor(kRed-6);

  theta_yz_North->SaveAs(SaveDir + "/Theta_yz_North.pdf");
  */
  /*
  TCanvas*  CRT_theta_xz_South = new TCanvas ("crttheta_xz1", "crttheta_xz1", 900, 700);

  CRT_Theta_xz_South->Draw("");
  CRT_Theta_xz_South->SetTitle("#theta_{xz} South");
  CRT_Theta_xz_South->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_xz_South->GetYaxis()->SetTitle("");
  CRT_Theta_xz_South->SetLineWidth(4);
  CRT_Theta_xz_South->SetLineColor(kRed-6);

  CRT_theta_xz_South->SaveAs(SaveDir + "/CRT_Theta_xz_south.pdf");

  TCanvas*  CRT_theta_xz_North = new TCanvas ("crttheta_xz2", "crttheta_xz2", 900, 700);

  CRT_Theta_xz_North->Draw("");
  CRT_Theta_xz_North->SetTitle("#theta_{xz} North");
  CRT_Theta_xz_North->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_xz_North->GetYaxis()->SetTitle("");
  CRT_Theta_xz_North->SetLineWidth(4);
  CRT_Theta_xz_North->SetLineColor(kRed-6);

  CRT_theta_xz_North->SaveAs(SaveDir + "/CRT_Theta_xz_North.pdf");

  TCanvas*  CRT_theta_yz_South = new TCanvas ("crttheta_xy1", "crttheta_xy1", 900, 700);

  CRT_Theta_yz_South->Draw("");
  CRT_Theta_yz_South->SetTitle("#theta_{yz} South");
  CRT_Theta_yz_South->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_yz_South->GetYaxis()->SetTitle("");
  CRT_Theta_yz_South->SetLineWidth(4);
  CRT_Theta_yz_South->SetLineColor(kRed-6);

  CRT_theta_yz_South->SaveAs(SaveDir + "/CRT_Theta_yz_south.pdf");

  TCanvas*  CRT_theta_yz_North = new TCanvas ("crttheta_xy2", "crttheta_xy2", 900, 700);

  CRT_Theta_yz_North->Draw("");
  CRT_Theta_yz_North->SetTitle("#theta_{yz} North");
  CRT_Theta_yz_North->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_yz_North->GetYaxis()->SetTitle("");
  CRT_Theta_yz_North->SetLineWidth(4);
  CRT_Theta_yz_North->SetLineColor(kRed-6);

  CRT_theta_yz_North->SaveAs(SaveDir + "/CRT_Theta_yz_North.pdf");
  */

}
