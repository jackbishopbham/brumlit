TGraph *gXS_0 = new TGraph("XS_0.txt");
TGraph *gXS_1 = new TGraph("XS_1.txt");
TGraph *gstopping = new TGraph("stopping_power.txt");
TGraph *gXS_t = new TGraph("XS_t.txt"); //Atomic Data and Nuclear Data Tables, Vol. 15, No. 1, January 1975
TGraph *gA0,*gA1,*gA2;
TGraph *gA0p,*gA1p,*gA2p;
TRandom3 *rndm = new TRandom3();
double d2r=atan(1.)/45.;
double pi=4.*atan(1.);
double mp=1.007825;
double mn=1.008665;
double mLi=7.016004;
double mBe=7.016929;
double amu=931.4941;
double Q=-1.64424;
double max_depth = 0.2;//mm
double InteractionE(double);
double InteractionTheta(double);
void makesumgraph();
double Eloss(double,double);
double ConvertCMtoLab(double Eint,double thetaCM);
double NeutronEnergy(double Eint, double thetalab, double thetaCM);
double Einitial=2.6;//MeV
double Eth=1.8804 ;//Threshold energy
double thetaCM,thetalab;
TH1F *h1 = new TH1F("h1","",800,1.8,2.6);
TH1F *hEn = new TH1F("hEn",";Neutron Energy [MeV];Yield",1000,0.0,1.0);
TH1F *hEn_source = new TH1F("hEn_source","",1000,0.0,1.0);
TH1F *hth = new TH1F("hth","",1800,0.0,180);
TH2F *hEth = new TH2F("hEth","",180,0,180,1000,0,1.);
TH2F *hEthCM = new TH2F("hEthCM","",180,0,180,1000,0,1.);
TH2F *hthth = new TH2F("hthth","COM vs lab",1800,0,180,1800,0,180);
bool secondks=false;
bool excited=false;
void brumlit() {
	gStyle->SetOptStat(0);
	const int nsamples = 1000000;
	makesumgraph();
	for(int n=0;n<nsamples;n++) {
		//Sample Ep_interaction
		double Eint = InteractionE(Einitial);
		h1->Fill(Eint);
		//Sample theta CM
		thetaCM = InteractionTheta(Eint);
		if(thetaCM==0) continue;//Skip this guy?
		thetalab = ConvertCMtoLab(Eint,thetaCM);
		hth->Fill((thetalab/d2r));
		double En = NeutronEnergy(Eint,thetalab,thetaCM);
		if(En==0) continue;//Skip this guy
//		if(!excited) continue;
		hthth->Fill(thetaCM/d2r,thetalab/d2r);
		hEn->Fill(En);
		hEth->Fill(thetalab/d2r,En);
//		if(!secondks)hEth->Fill(thetalab/d2r,En);
		hEthCM->Fill(thetaCM/d2r,En);
	}
	TCanvas *c = new TCanvas();
	hth->Draw();
	TCanvas *c2 = new TCanvas();
	hEn->Draw();
	TGraph *g1 = new TGraph("sourceEp2600.dat");
	for(int x=0;x<g1->GetN();x++) {
		hEn_source->SetBinContent(hEn_source->GetXaxis()->FindBin(1e-3*g1->GetX()[x]),g1->GetY()[x]);
	}
	hEn_source->Scale(hEn->Integral(0,1000)/hEn_source->Integral(0,1000));
	hEn_source->SetLineColor(kRed);
	hEn_source->Draw("HISTSAME");
	TLegend *leg = new TLegend(0.7,0.9,0.7,0.9);
	leg->AddEntry(hEn,"BrumLiT");
	leg->AddEntry(hEn_source,"Rendimiento");
	leg->Draw("SAME");
	TCanvas *c3 = new TCanvas();
	c3->Divide(2,1);
	c3->cd(1);
	hEth->Draw("COLZ");
	c3->cd(2);
	hEthCM->Draw("COLZ");
	TCanvas *c4 = new TCanvas();
	hthth->Draw("COLZ");
	TCanvas *c5 = new TCanvas();
//	cout<<gXS_0->Eval(1.94999)<<endl;
	gXS_t->Draw();
	TCanvas *c6 = new TCanvas();
	h1->Draw();
}

double InteractionE(double Ep) {
	//Rejection sampling for z_interaction based on XS(Ep) and stopping power
	//Then return Ep for z
	bool trials=true;
	double XSmax = gXS_t->Eval(Einitial);//Max XS is at the highest E
	XSmax = 1000;
	int loopcounter=0;
	while(trials) {//Randomly choose a depth into the target (rather than energy to keep same target thickness)
		double trial_z = max_depth*rndm->Rndm();
		//Calc XS at Z
		double energy = Eloss(Ep,trial_z);
		if(energy<=Eth) continue;
		double XS = gXS_t->Eval(energy);
		if(energy<1.925) {//NEED TO ADD IN THE TWO OPTIONS
			double C0=6;//From Lee and Zhou
			double A=164.913;//mb MeV/sr
			double x=C0*sqrt(1-Eth/energy);
			XS=4.*pi*(A*x)/(energy*(1+x)*(1+x));
		}
//		cout<<energy<<"\t"<<XS<<endl;
//		cout<<loopcounter<<"\t"<<energy<<"\t"<<XS<<"\t"<<trial_z<<"\t"<<XSmax<<endl;
//		cout<<"RATIO: "<<energy<<"\t"<<gXS_1->Eval(energy)/gXS_0->Eval(energy)<<endl;
		if(energy>2.5 && rndm->Rndm()<gXS_1->Eval(energy)/gXS_0->Eval(energy)) {
//			cout<<"Excited"<<endl;
			excited=true;
		}
		else excited=false;
		if(XS>XSmax) {cout<<"Max XS is exceeded!!"<<endl;}//
		if(XS>XSmax*rndm->Rndm()) {
			return energy;
		}
		loopcounter++;
		if(loopcounter>1000) {
			cout<<"Rejection sampling not working properly - 1000 samples taken"<<endl;
			return 0;
		}
	}
	return 0;
}

void makesumgraph() {
//	gXS_t = gXS_0;
	ifstream in;
	in.open("XS_ang0.txt");
	double inA0[13],inA1[13],inA2[13],inE[13];
	int count=0;
	while(in.good() && count<13) {
		in>>inE[count]>>inA0[count]>>inA1[count]>>inA2[count];
		count++;
	}
	gA0 = new TGraph(13,inE,inA0);
	gA1 = new TGraph(13,inE,inA1);
	gA2 = new TGraph(13,inE,inA2);
//
	ifstream in2;
	in2.open("XS_ang1.txt");
	double inA0p[2],inA1p[2],inA2p[2],inEp[2];
	count=0;
	while(in2.good() && count<2) {
		double readE,readA0,readA1,readA2;
//		in2>>readE>>readA0>>readA1>>readA2;
		in2>>inEp[count]>>inA0p[count]>>inA1p[count]>>inA2p[count];
//		cout<<readE<<"\t"<<readA0<<"\t"<<readA1<<"\t"<<readA2<<endl;
		count++;
	}
	gA0p = new TGraph(2,inEp,inA0p);
	gA1p = new TGraph(2,inEp,inA1p);
	gA2p = new TGraph(2,inEp,inA2p);
}

double Eloss(double Ep, double z) {//Energy of proton of initial energy Ep after z mm in nat Li
	double range=0.0237*Ep*Ep+0.0369*Ep-0.0102;//mm
//	cout<<"INitial range\t"<<range<<endl;
	if(z>range) return 0;//Fully stopped before Z
	range-=z;
	double Eout=-7.98*range*range+10.0761*range+0.6029;//MeV
//	cout<<Ep<<"\t"<<range<<"\t"<<z<<"\t"<<Eout<<endl;
	return Eout;
}

double InteractionTheta(double Ep) {//COM theta
//	return (100+10*rndm->Rndm())*d2r;
	double random_theta = acos(-1+2*rndm->Rndm());
	random_theta = 2.*pi*rndm->Rndm();
	double A[3]={1,0,0};
	bool sample=true;
	double XSmax=3;
	int counter=0;
	if(Ep<1.95) {//Below where we have data it is isotropic
		A[0]=1.;
		A[1]=0.;
		A[2]=0.;
	}
	else {
		if(!excited) {
			A[0]=gA0->Eval(Ep);
			A[1]=gA1->Eval(Ep);
			A[2]=gA2->Eval(Ep);
		}
		if(excited) {
			A[0]=gA0p->Eval(Ep);
			A[1]=gA1p->Eval(Ep);
			A[2]=gA2p->Eval(Ep);
		}
	}
//	cout<<Ep<<"\t"<<excited<<"\t"<<A[0]<<"\t"<<A[1]<<"\t"<<A[2]<<endl;
	while(sample) {
		double XS_sample = XSmax*rndm->Rndm();
		random_theta = 2.*pi*rndm->Rndm();
		double XS=sin(random_theta)*(A[0]*1+A[1]*cos(random_theta)+A[2]*0.5*(3*cos(random_theta)*cos(random_theta)-1));
		if(XS>XSmax) cout<<"Sampling error for angle"<<endl;
//		cout<<XS<<"\t"<<XSmax<<endl;
		if(XS_sample<XS) {
			return random_theta;//Isotropic for now
		}
		counter++;
		if(counter>1000) {
			cout<<"Over 1000 samples taken - no solution"<<endl;
			return 0;
		}
	}
}

double ConvertCMtoLab(double Ep, double thetaCM) {//Convert thetaCM to lab
	double Ex=0;
	if(excited) Ex=0.477;
	double gamma = sqrt(((mp*mn)/(mBe*mLi))*(Ep/(Ep+(Q-Ex)*(1.+mp/mLi))));//Calculate gamma for the conversion to lab COM
	gamma = gamma*((mp+mLi)/(mn+mBe));//Account for change in the CM velocity
	double thetalab = atan2(sin(thetaCM),(cos(thetaCM)+gamma));//
	secondks = false;
	if(gamma*cos(thetaCM)<=-1.) secondks=true;//Neutron is going bac
//	thetalab=(12*rndm->Rndm())*d2r;//TEST
//	if(rndm->Rndm()>0.5) secondks=true;//TEST
//	cout<<secondks<<"\t"<<thetaCM/d2r<<"\t"<<gamma*cos(thetaCM)<<"\t"<<gamma<<endl;
	return thetalab;
}

double NeutronEnergy(double Ep, double theta, double thetaCM) {//Get neutron energy for a given Ep and theta
/*
	double a,b,c;
	double p0=sqrt(2.*mp*Ep);
	double Etot=Ep+Q;
	a=1.+mBe/mn;
	b=2.*p0*cos(theta);
	c=p0*p0-2.*mBe*Etot;
	double x=(-b-sqrt(b*b-4.*a*c))/(2.*a);
	double thetacrit = acos(sqrt((a*c)/(p0*p0)));
	secondks=false;
	if(theta>thetacrit) secondks=true;
//	cout<<theta/d2r<<"\t"<<b*b-4.*a*c<<endl;
//	Catch second kinematic solution
//	hEth->Fill(thetacrit/d2r,(b*b/(4*a*a))/(2.*mn),1e5);
//	cout<<thetacrit/d2r<<endl;
	if(secondks) {
		x=(-b+sqrt(b*b-4.*a*c))/(2.*a);
	}
	double En = (x*x)/(2.*mn);
//	Krane pg 382 eq:
	double numer = sqrt(mp*mn*Ep)*cos(theta)+sqrt(mp*mn*Ep*cos(theta)*cos(theta)+(mLi+mn)*(mLi*Q+(mLi-mp)*Ep));
	if(secondks) numer = sqrt(mp*mn*Ep)*cos(theta)-sqrt(mp*mn*Ep*cos(theta)*cos(theta)+(mLi+mn)*(mLi*Q+(mLi-mp)*Ep));
	double denom = mLi+mn;
	En = TMath::Power(numer/denom,2);
//
*/
	double En;
	double Ex=0;
	if(excited) Ex=0.431;//MeV
	double ECM = mLi*Ep/(mLi+mp)+Q-Ex;
	if(ECM<0) {
		cout<<"Below threshold! Something is wrong!"<<endl;
		return 0;
	}
	double p_nT = sqrt(2.*ECM*mn*(mBe/(mn+mBe)));
	double phi=2*pi*rndm->Rndm();
	double pn[3]={p_nT*sin(thetaCM)*cos(phi),p_nT*sin(thetaCM)*sin(phi),p_nT*cos(thetaCM)};//COM neutron momentum
	pn[2]+=(mp/(mp+mLi))*sqrt(2.*mp*Ep);//Add the momentum boost for the lab frame
	double En_new = 0;
	for(int i=0;i<3;i++) En_new+=pn[i]*pn[i];
	En_new /= (2.*mn);
	En = En_new;
	theta = acos(pn[2]/p_nT);
	return En;
}
