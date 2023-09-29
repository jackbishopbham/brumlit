void angles() {
	TF1 *f1 = new TF1("f1","[0]+[1]*cos(x)+[2]*0.5*(3*cos(x)*cos(x)-1)",0,6);
	f1->SetParameter(0,1);
	f1->SetParameter(1,0);
	f1->SetParameter(2,0);
	f1->Draw();
}
