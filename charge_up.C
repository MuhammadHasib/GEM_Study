#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include<stdlib.h>
#ifdef __CINT__
#define ROOT_GuiTypes 0
#endif
using namespace std;



int charge_up()
{
gStyle->SetOptStat(0000);
const n=   2500;
Double_t a[n];

ifstream infile, in;
string line;
string input_filename, output_filename, output_graph;

int m =16, l, f=0, xy1=0; ;   //Change m for the no. of input files
string ext1 = "CUp_20072015_";
string ext2 =".mca";
string ext3 = ".pdf";
TCanvas *c1= new TCanvas("c1"," "); 		//Canvas definition
Double_t m_e[500],s_e[500], m_s[500], s_s[500];	
Double_t H[500], M[500], S[500]; 	

int input_files =1;

  for(int two =0; two<1; two++){  
    for(int three =0; three<3; three++) {		
      for (int i=0;i<10;i++) {
        for (int j=0;j<10;j++) {
	if(input_files<254 ) {   ///input files
	stringstream ssi;
	ssi<<two<<three<<i<<j;
	input_filename = ext1 + ssi.str() + ext2;
	output_graph = ssi.str() + ext3; 
	infile.open(input_filename.c_str(),ios::in|ios::out);	
	if (!infile)
	{
		cout << " ERROR: " <<input_filename.c_str()<<" file not found" << endl;
		cout << endl;

		}
	else
		{
		
		printf("\n\n\n\n");
		cout<<"Input_file_no. = "<<input_files<<"\n";
		++input_files;
		cout<<"File " <<input_filename.c_str()<<" opened successfully" << endl;
		TH1F *h1 = new TH1F("h1","Pulse Height Spectrum",2000,0,400);
		h1->GetXaxis()->SetTitle("Charge (ADC channels)");
		h1->GetYaxis()->SetTitle("Counts");
		h1->GetXaxis()->CenterTitle();
  	 	h1->GetYaxis()->CenterTitle();
		TSpectrum *s = new TSpectrum(2);
		Double_t par[3000];
		
		Double_t bounds[50], xlow1, xlow2, xup1, xup1;
		int count=1;
		int line_number=1;

		while (getline(infile, line))
       		{
		line_number++;

		if(line_number == 11) {
        	istringstream ss(line);
		char cg1, cg2;
		string gr1, gr2,gr3, gr4, g5;
		int val1, val2, val3;
		cout<<line<<endl;
		ss >> gr1>> gr2>> gr3 >>val1>>cg1>>val2>>cg2>>val3;
		H[xy1]=val1;
		M[xy1]=val2;
		S[xy1]=val3;
		cout<<"Timing information of "<<input_filename.c_str()<<" is "<<endl;
		cout<<"File No. "<<xy1+1<<endl;
		cout<<H[xy1]<<"h "<<M[xy1]<<"m "<<S[xy1]<<"s "<<endl;
		++xy1;
		
			}


		if (line_number>130) { //reading entries in the file 
        	istringstream ss(line);
        	string name;
        	int var1;
        	ss >>var1;
		h1->SetBinContent(count, var1);
		++count;            //counter
			} // Closing brace for line number
			
       		}

		c1->cd();
		h1->Rebin(10); 
		h1->Draw();

		//===========================================>

		Int_t nfound = s->Search(h1,4,"new"); //searching peaks
		Float_t *xpeaks = s->GetPositionX(); //assigning peak positions to array element
   		npeaks = 0;
   		Float_t *xpeaks = s->GetPositionX();

		for (int d=0; d<nfound; d++)

		{
		Float_t xy = xpeaks[d];
		bounds[d] = xy;

			}

		//================================================>

		xlow1 = (3.6*bounds[1] - bounds[0])/2.5; 
		xup1 = (1.5*bounds[1] + bounds[0])/2.5;
		xlow2 = (1 *bounds[0] + bounds[1])/2;
		xup2 = (3*bounds[0] - bounds[1])/1.5;

		//=================================================>

		TF1 *f1 = new TF1("f1","gaus",xlow1,xup1);		//escape fit function definition
		TF1 *f2 = new TF1("f2","gaus",xlow2,xup2);		//signal fit function definition
		TF1 *total = new TF1("total","gaus(0)+gaus(3)",xlow1,xup2);



   		for (p=0;p<nfound;p++) {
      		Float_t xp = xpeaks[p];
      		Int_t bin = h1->GetXaxis()->FindBin(xp);
      		Float_t yp = h1->GetBinContent(bin);
      		par[3*npeaks+0] = yp;
      		par[3*npeaks+1] = xp;
      		par[3*npeaks+2] = 2;
		f1->FixParameter(1, par[1]);
		f2->FixParameter(1, par[4]);
      		npeaks++;
   		}
		//===================================================>


		printf("Found %d useful peaks to fit \n\n",npeaks);
   		printf("Now fitting: Be patient \n\n");	

		//=================================================>

		 h1->Fit(f1,"R");
	         h1->Fit(f2,"R+");
		 f1->GetParameters(&par[0]);
   		 f2->GetParameters(&par[3]);
		 total->SetParameters(par);
		 h1->Fit(total,"R+");

		 string aaa = "45.txt";	
		if (aaa.compare( input_filename.c_str())) //this loop is just remove one bad point in the escape peak of file 45.txt
		{

		m_e[f] = total->GetParameter(1);  //saving mean of escape peak into the array
		s_e[f] = total->GetParameter(2);    //saving sigma of escape peak into the array
		m_s[f] = total->GetParameter(4);    //saving mean of sinal peak into the array
		s_s[f] = total->GetParameter(5); //saving sigma of the signal peak into the array

			}
		 else 

		{
		m_e[f] = f1->GetParameter(1);      //saving mean of escape peak into the array 
		s_e[f] = total->GetParameter(2);    //saving sigma of escape peak into the array
		m_s[f] = total->GetParameter(4);    //saving mean of sinal peak into the array
		s_s[f] = total->GetParameter(5);    //saving sigma of the signal peak into the array

			}

		 f++;

		//====================================================>
		c1->Update();
		c1->SaveAs(output_graph.c_str());
		h1->Reset();

		}	
	

	infile.close();
	
	}
	
	}
      }   
    }
   }
 
//==================================================================>Drawing first graph
//Calculations

Double_t H_M[500], Time[500]; 
for (int ab = 0; ab <input_files; ab++) { 

	H_M[ab] = H[ab]*60;
	Time[ab] = H_M[ab] + M[ab] + (S[ab]/60) ; //subtracted by 600
	//cout<<Time[ab]<<endl;

}
TCanvas *c2 = new TCanvas("c2"," Charge_up", 200 ,10,700,500);
TPad *pad = new TPad("pad","",0,0,1,1);
pad->SetFillColor(42);
pad->SetGrid();
pad->Draw();
pad->cd();

// draw a frame to define the range
   TH1F *hr1 = pad->DrawFrame(630, 140, 850, 170); //(-0.4,0,1.2,150);
   hr1->SetXTitle("Time (min)");
   hr1->SetYTitle("Peak Energy (ADC Counts)");
   hr1->GetXaxis()->CenterTitle();
   hr1->GetYaxis()->CenterTitle();
   pad->GetFrame()->SetFillColor(21);
   pad->GetFrame()->SetBorderSize(12);
TGraph *graph1 = new TGraph (input_files-2, Time, m_s);
graph1->SetLineColor(4);
graph1->Draw("LP");

//TMultiGraph *multigraph = new TMultiGraph();
//graph1->GetXaxis()->SetTitle("time (min)");
//graph1->GetYaxis()->SetTitle("Peak Energy");
//c2->SaveAs("Charge_Up.pdf");
//======================================================================>

string line1;
int xyz=0, line_number1 = 1;
Double_t H_2[2500], M_2[2500], Temp_2[2500], Pressure[2500], Humidity[2500];
infile.open("temp.txt");	
	if (!infile)
	{
		cout << " ERROR: file not found" << endl;
		cout << endl;

		}
	else
		{

		while (getline(infile, line1))
       		{
		line_number1++;
		if(line_number1 > 22 && line_number1< 278) {        //cut over here
        	istringstream ssk(line1);
		char cg1, cg2;
		int var1, var2, var3, var4, var5, var6, var7, var8;
		Double_t var9, var10, var11, var12, var13, var14, var15;
		ssk >> var1>> var2>>var3 >>var4>>var5>>var6>>var7>>var8>>var9>>var10>>var11>>var12>>var13>>var14>>var15;
		H_2[xyz]=var4;
		M_2[xyz]=var5;
		Temp_2[xyz]=var8;
		Pressure[xyz] = var14 ;
		cout<<"Pressure = "<<Pressure[xyz]<<endl;
		Humidity[xyz] = var15 ;
		cout<<"Humidity = "<<Humidity[xyz]<<endl;
		++xyz;
		}
						
		}

	}	
//=================================================================>
//Calculations

Double_t H_M2[2500], Time_2[2500], T_1[2500], Temp_1[2500], Temp_b_Press[2500]; 
for (int abc = 0; abc <xyz-1; abc++) { 

	H_M2[abc] = H_2[abc]*60;
	Time_2[abc] = H_M2[abc] + M_2[abc] ; //timeing information
	Temp_1[abc] =  (Temp_2[abc]/100);  //temperature divided by 100
	Temp_b_Press[abc] = Temp_1[abc]/Pressure[abc];
	cout<<"Temp_b_Press = "<<Temp_b_Press[abc]<<endl;

}

//graph2->GetXaxis()->SetRange(500,900);
c2->cd();
TPad *overlay = new TPad("overlay","",0,0,1,1);
   overlay->SetFillStyle(4000);
   overlay->SetFillColor(0);
   overlay->SetFrameFillStyle(4000);
   overlay->Draw();
   overlay->cd();
TGraph *graph2 = new TGraph (xyz-2, Time_2, Temp_1);
graph2->SetLineColor(2);

Double_t xmin = pad->GetUxmin();
   Double_t ymin = 18;
   Double_t xmax = pad->GetUxmax();
   Double_t ymax = 26;
   TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
   hframe->GetXaxis()->SetLabelOffset(99);
   hframe->GetYaxis()->SetLabelOffset(99);
   //hframe->GetYaxis()->CenterTitle();
   //hframe->SetYTitle("Temperature");
   graph2->Draw("LP");

//Draw an axis on the right side
   TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
   axis->SetLineColor(kRed);
   axis->SetTitle("Temperature (C)");
   axis->CenterTitle();   //modified
   axis->SetTitleColor(kRed);
   axis->SetLabelColor(kRed);
   axis->Draw();
c2->SaveAs("Charge_Up2.pdf");
//===================================================================================>
TCanvas *c3 = new TCanvas("c3"," Peak Energy versus Pressure", 200 ,10,700,500);
c3->cd();
TGraph *graph3 = new TGraph (xyz-2, Time_2, Pressure);  //Time_2, Pressure
graph3->SetLineColor(4);
graph3->Draw("ALP");
graph3->GetXaxis()->SetTitle("Pressure");
graph3->GetYaxis()->SetTitle("Peak Energy (ADC Counts)");
graph3->GetXaxis()->CenterTitle();
graph3->GetYaxis()->CenterTitle();
c3->SaveAs("Peak_Pressure.pdf");
//=================================================================================>

TCanvas *c4 = new TCanvas("c4","Peak Energy versus Humidity", 200 ,10,700,500);
c4->cd();
TGraph *graph4 = new TGraph (xyz-2, Time_2, Humidity);
graph4->SetLineColor(3);
graph4->Draw("ALP");
graph4->GetXaxis()->SetTitle("Humidity");
graph4->GetYaxis()->SetTitle("Peak Energy (ADC Counts)");
graph4->GetXaxis()->CenterTitle();
graph4->GetYaxis()->CenterTitle();
c4->SaveAs("Peak_Humidity.pdf");

//=====================================================================================>
TCanvas *c5 = new TCanvas("c5"," Peak_Energy versus Temperature/Pressure", 200 ,10,700,500);
c5->cd();
TGraph *graph5 = new TGraph (xyz-2, Time_2, Temp_b_Press);
graph5->SetLineColor(4);
graph5->Draw("ALP");
graph5->GetYaxis()->SetTitle("Temperature/Pressure");
graph5->GetYaxis()->SetTitle("Peak Energy (ADC Counts)");
graph5->GetXaxis()->CenterTitle();
graph5->GetYaxis()->CenterTitle();
c5->SaveAs("Peak_Temp_Press.pdf");


} // main program closing brace





