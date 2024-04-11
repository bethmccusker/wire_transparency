#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<TMath.h>
#include<TTree.h>
#include<TFile.h>
#include<TROOT.h>
#include<iostream>

void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour)
{
  c1->cd();
  TList *outline = new TList;
  TPolyLine3D *p1 = new TPolyLine3D(4);
  TPolyLine3D *p2 = new TPolyLine3D(4);
  TPolyLine3D *p3 = new TPolyLine3D(4);
  TPolyLine3D *p4 = new TPolyLine3D(4);
  p1->SetLineColor(colour);
  p1->SetLineWidth(2);
  p1->Copy(*p2);
  p1->Copy(*p3);
  p1->Copy(*p4);
  outline->Add(p1);
  outline->Add(p2);
  outline->Add(p3);
  outline->Add(p4); 
  TPolyLine3D::DrawOutlineCube(outline, rmin, rmax);
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();

}

void DrawTrack(TCanvas *c1, double* point1, double* point2, int colour)
{
  TPolyLine3D *track = new TPolyLine3D(2);
  track->SetLineWidth(2);
  track->SetLineColor(colour);
  track->SetPoint(0,point1[0], point1[1], point1[2]);
  track->SetPoint(1,point2[0],point2[1],point2[2]);
  track->Draw();
}


void event_display(double *track_start,double *track_end,double index){

  TCanvas *c1 = new TCanvas("c1","");
 
  double  rmin[3] {-200,-200,0}; 
  double  rmax[3] {-200,200,500};

  DrawCube(c1, rmin, rmax, kRed-6 );
 
  DrawTrack(c1,track_start,track_end,kBlue-6);

  double arrow_start_x = - 0.7;
  double arrow_start_y = -0.8;
  double arrow_length = 0.2; 
  TLatex coord_y;
  TLatex coord_z;
  coord_y.DrawLatex(arrow_start_x-0.015,arrow_start_y+0.2+0.03, "y");
  coord_z.DrawLatex(arrow_start_x+0.13 +0.02, arrow_start_y-0.02, "z");
  TArrow *arrow_y = new TArrow(arrow_start_x, arrow_start_y, arrow_start_x, arrow_start_y+0.2, 0.01, "|>");
  TArrow *arrow_z = new TArrow(arrow_start_x, arrow_start_y, arrow_start_x+0.13, arrow_start_y, 0.01, "|>");
  arrow_y->SetLineWidth(2);
  arrow_z->SetLineWidth(2);
  arrow_y->Draw();
  arrow_z->Draw();


  TView3D *view = (TView3D*) TView::CreateView(1);

  double c[3] = { 0, 0, 250 };
  double s[3] = { 800, -800, 800 };
  view->SetRange(-600, -600, -300, 600, 600, 900);
  view->DefineViewDirection(s, c,
			    1, 0,
			    0, 1,
			    0, 1,
			    view->GetTnorm(),
			    view->GetTback());

  const std::string SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";
   std::string event = SaveDir + ("/event_display_" + std::to_string(index) + ".pdf");
   c1->SaveAs(event.c_str());  
 //   c1->SaveAs(SaveDir + "/event_display.pdf");

   //  delete c1;
}

void event_display_test(){
  // double track_start[3] {-200,70,-100};                                                                                                                                                             
  //  double track_end[3] {-200,-70,600};  
  //  event_display(double * track_start, double *track_end);
}
