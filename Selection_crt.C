#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

{
  gROOT->Reset();
  const TString SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";

  gStyle->SetOptStat(0); //Removing the stats box
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
  TH2F *x_y_Hit_1 = new TH2F("xy1","xy1",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_2 = new TH2F("xy3","xy3",70, -360, 360, 70, -360, 360);
  TH2F *x_y_AllHits_South = new TH2F("xy5","xy5",70, -360, 360, 70, -360, 360);
  TH2F *x_y_AllHits_North = new TH2F("xy6","xy6",70, -360, 360, 70, -360, 360);
  TH1F *TPC_Theta_xz =new TH1F("theta_xz_1", "theta_xz_1", 50, -90, 90);
  TH1F *TPC_Theta_yz =new TH1F("theta_yz_1", "theta_yz_1", 50, -90,90);
  TH1F *CRT_Theta_xz =new TH1F("crttheta_xz_1", "crttheta_xz_1", 50, -200, 200);
  TH1F *CRT_Theta_yz =new TH1F("crttheta_yz_1", "crttheta_yz_1", 50, -200,200);
  int nEntries = myTree->GetEntries(); // Get the number of entries in this tree
  for (int iEnt = 0; iEnt < nEntries; iEnt++) { //Big loop of entries
    myTree->GetEntry(iEnt); // Gets the next entry (filling the linked variables)
   
    //*************** Selection start ***************// 

    //CRT hit 1/2 loop (used to define cuts)
    bool HasHitInNorthandSouth=false;
    for (size_t c=0; c < ct_x1->size(); ++c) {
      if(((ct_z1->at(c)>-200 && ct_z1->at(c)<-150) && (ct_z2->at(c)>750 && ct_z2->at(c)<800)) || ((ct_z1->at(c)>750 && ct_z1->at(c)<800) && (ct_z2->at(c)>-200 && ct_z2->at(c)<-150))){                           
	if((ct_y1->at(c)>-360 && ct_y1->at(c)<360) && (ct_y2->at(c)>-360 && ct_y2->at(c)<360)){                                                                                                                   
	  if(((ct_x1->at(c)>-300 && ct_x1->at(c)<-250)&& (ct_x2->at(c)>-300 && ct_x2->at(c)<-250)) || ((ct_x1->at(c)>250 && ct_x1->at(c)<300)&& (ct_x2->at(c)>250 && ct_x2->at(c)<300))){ 
	    HasHitInNorthandSouth=true;
	  }
	}
      }
    }

    // TPC loop 1 (used to define cuts)
    bool HasTrkType4= false;
    for(size_t j=0; j<muontrk_type->size(); ++j){ //Looping over TPC tracks
      if(muontrk_type->at(j)==4 || muontrk_type->at(j)==5){//Selecting track type 4 + 5 
	HasTrkType4=true;
      }// if
    }//for


    //CRT loop 1 (used to define cuts)
    bool IsInHitRegionSouth=false;
    bool IsInHitRegionNorth=false;
    for(size_t i=0; i < chit_x->size(); ++i){ //Looping over CRT hits
      if(chit_z->at(i)>-200 && chit_z->at(i)<-150 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360 ){ //Selectiong South CRT
	if((chit_x->at(i)>-200 && chit_x->at(i)<-150)||(chit_x->at(i)<200 && chit_x->at(i)>150)){ // Cut on x region
	  IsInHitRegionSouth=true;
	} //CRT if
      }//x cut if
      if(chit_z->at(i)>750 && chit_z->at(i)<800 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360){ //Selecting north CRT
	if((chit_x->at(i)>-200 && chit_x->at(i)<-150)||(chit_x->at(i)<200 && chit_x->at(i)>150)){ //Cut on x region
	  IsInHitRegionNorth=true;
	}// CRT if
      }// x cut if
    }//CRT hits for
     
    if(HasTrkType4 && HasHitInNorthandSouth){  // Cuts required
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
	 
	if(muontrk_type->at(j)==4 || muontrk_type->at(j)==5){
	  TPC_Theta_xz->Fill(modified_theta_xz);
	  TPC_Theta_yz->Fill(modified_theta_yz);
	}
      }//TPC trk loop end

      //CRT loop 2 (Used to fill histograms)
      for(size_t i=0; i < chit_x->size(); ++i){ //Looping over CRT hits
	if(chit_z->at(i)>-200 && chit_z->at(i)<-150 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360 ){ //Selectiong South CRT 
	  if((chit_x->at(i)>-200 && chit_x->at(i)<-150)||(chit_x->at(i)<200 && chit_x->at(i)>150)){ //Cut on x region
	    x_y_AllHits_South->Fill(chit_x->at(i),chit_y->at(i));
	  } //South CRT
	}// x region
	if(chit_z->at(i)>750 && chit_z->at(i)<800 && chit_x->at(i)>-360 && chit_x->at(i)<360 && chit_y->at(i)>-360 && chit_y->at(i)<360){ //Selecting north CRT  
	  if((chit_x->at(i)>-200 && chit_x->at(i)<-150)||(chit_x->at(i)<200 && chit_x->at(i)>150)){ //Cut on x region 
	    x_y_AllHits_North->Fill(chit_x->at(i),chit_y->at(i));
	  }// North CRT
	}// x region
      }//CRT hit loop end

      //CRT hit 1/2 loop (used to fill histograms)
      for (size_t c=0; c < ct_x1->size(); ++c) {
	double dx;
	double dy;
	double dz;
	double CRT_theta_xz;
	double CRT_theta_yz;
	if(((ct_z1->at(c)>-200 && ct_z1->at(c)<-150) && (ct_z2->at(c)>750 && ct_z2->at(c)<800)) || ((ct_z1->at(c)>750 && ct_z1->at(c)<800) && (ct_z2->at(c)>-200 && ct_z2->at(c)<-150))){
	  if((ct_y1->at(c)>-360 && ct_y1->at(c)<360) && (ct_y2->at(c)>-360 && ct_y2->at(c)<360)){
	    if(((ct_x1->at(c)>-300 && ct_x1->at(c)<-250)&& (ct_x2->at(c)>-300 && ct_x2->at(c)<-250)) || ((ct_x1->at(c)>250 && ct_x1->at(c)<300)&& (ct_x2->at(c)>250 && ct_x2->at(c)<300))){ 
	      x_y_Hit_1->Fill(ct_x1->at(c),ct_y1->at(c));
	      x_y_Hit_2->Fill(ct_x2->at(c),ct_y2->at(c));

	      
	      dx =ct_x2->at(c) - ct_x1->at(c);
	      dy =ct_y2->at(c) - ct_y1->at(c);
	      dz =ct_z2->at(c) - ct_z1->at(c);
		
	      CRT_theta_xz= atan2(dx,dz)*(180/TMath::Pi());
	      /*
		if(CRT_theta_xz<-12 && CRT_theta_xz>-24){
		cout<<dx<< "   dx"<<endl;
		cout<<dz<< "   dz"<<endl;
		cout<<ct_x2->at(c)<< "    x2"<<endl;
		cout<<ct_x1->at(c)<< "    x1"<<endl;
		cout<<CRT_theta_xz<< "   theta_xz"<<endl;
		} */
	      CRT_theta_yz= atan2(dy,dz)*(180/TMath::Pi());
	      /*
		if(CRT_theta_yz<-12 &&CRT_theta_yz>-24){
		cout<<dy<< "   dy"<<endl; 
		cout<<dz<< "   dz"<<endl;
		cout<<ct_y2->at(c)<< "    y2"<<endl;
		cout<<ct_y1->at(c)<< "    y1"<<endl;
		cout<<CRT_theta_yz<< "   theta_yz"<<endl;
		}*/
	      CRT_Theta_xz->Fill(CRT_theta_xz);
	      CRT_Theta_yz->Fill(CRT_theta_yz);
	      
	    }//x if
	  }//y if 
	}//z if
      }// end of CRT hit 1/2 loop

	
    }// end of required cuts if 

    //*************** Selection end ***************//


  } //Big loop of entries
 
    //*************** Making pretty histograms ********************//
  
  TCanvas*  x_y_1 = new TCanvas ("xy1", "xy1", 900, 700);

  x_y_Hit_1->Draw("COLZ");
  x_y_Hit_1->SetTitle("CRT Muon Hit 1");
  x_y_Hit_1->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_1->GetYaxis()->SetTitle("y Position (cm)");

  x_y_1->SaveAs(SaveDir + "/x_y_Hit_1.pdf");  


  TCanvas*  x_y_3 = new TCanvas ("xy3", "xy3", 900, 700);

  x_y_Hit_2->Draw("COLZ");
  x_y_Hit_2->SetTitle("CRT Muon Hit 2");
  x_y_Hit_2->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_2->GetYaxis()->SetTitle("y Position (cm)");

  x_y_3->SaveAs(SaveDir + "/x_y_Hit_2.pdf");
 

  
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
  

  TCanvas*  TPC_theta_xz = new TCanvas ("theta_xz1", "theta_xz1", 900, 700);

  TPC_Theta_xz->Draw("");
  TPC_Theta_xz->SetTitle("#theta_{xz} TPC");
  TPC_Theta_xz->GetXaxis()->SetTitle("Angle (degrees)");
  TPC_Theta_xz->GetYaxis()->SetTitle("");
  TPC_Theta_xz->SetLineWidth(4);
  TPC_Theta_xz->SetLineColor(kRed-6);

  TPC_theta_xz->SaveAs(SaveDir + "/TPC_Theta_xz.pdf");
  
  TCanvas*  TPC_theta_yz = new TCanvas ("theta_xy1", "theta_xy1", 900, 700);

  TPC_Theta_yz->Draw("");
  TPC_Theta_yz->SetTitle("#theta_{yz} TPC");
  TPC_Theta_yz->GetXaxis()->SetTitle("Angle (degrees)");
  TPC_Theta_yz->GetYaxis()->SetTitle("");
  TPC_Theta_yz->SetLineWidth(4);
  TPC_Theta_yz->SetLineColor(kRed-6);

  TPC_theta_yz->SaveAs(SaveDir + "/TPC_Theta_yz.pdf");
 
  
  TCanvas*  CRT_theta_xz = new TCanvas ("crttheta_xz1", "crttheta_xz1", 900, 700);

  CRT_Theta_xz->Draw("");
  CRT_Theta_xz->SetTitle("#theta_{xz} CRT");
  CRT_Theta_xz->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_xz->GetYaxis()->SetTitle("");
  CRT_Theta_xz->SetLineWidth(4);
  CRT_Theta_xz->SetLineColor(kRed-6);

  CRT_theta_xz->SaveAs(SaveDir + "/CRT_Theta_xz.pdf");


  TCanvas*  CRT_theta_yz = new TCanvas ("crttheta_xy1", "crttheta_xy1", 900, 700);

  CRT_Theta_yz->Draw("");
  CRT_Theta_yz->SetTitle("#theta_{yz} CRT");
  CRT_Theta_yz->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_yz->GetYaxis()->SetTitle("");
  CRT_Theta_yz->SetLineWidth(4);
  CRT_Theta_yz->SetLineColor(kRed-6);

  CRT_theta_yz->SaveAs(SaveDir + "/CRT_Theta_yz.pdf");
  
}
