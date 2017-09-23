/* Author: Cesar Serna
 * Date: 31/oct/2015
 * 
 * In this file are all the general routines used for classical simulations, every program must use
 * this routines, in this way there will not exist different versions of the same.
 */
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define MAX_N 1000
using namespace std;
int previous_val_label_i,previous_val_label_j; //these are used to restart to a previous configuration of in a trial movement
int label_i_changed,label_j_changed;
double potential_energy(double);
double potential_energyOneSite(double);
double distance2OneSite(int,int);
/****************************************************
*			Declaration of all variables			*
*****************************************************/
int N, NMOVE, NSUB, load_prev;
double T, RHO, DISPL=0.5, DISPLONESITE=0.5,S,SS,maxdS,rd,rc,epsilon;
double RC2,EPSILON; // the size of the site two the square
double xs[MAX_N], ys[MAX_N], zs[MAX_N]; // save the direction of each site 
double x[MAX_N], y[MAX_N], z[MAX_N]; // Used when simulating sites (associating fluids)
int label[MAX_N]; // used to label when a siteis overlaped with another site
double DIST; //defines the distance between the site and the center of the particle
double PI;
double sqrt2PI,lb; // this variables are used for semiclassical MC
long seed,*pseed;
string filename_prev_conf;
	
/****************************************************
*		Declaration of the random generator			*
*****************************************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 379
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran2(long *idum){
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0) {
		if (-(*idum) < 1)
			*idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
				if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
/****************************************************
*		Declaration of general functions			*
*****************************************************/
void print_label(string label){
	cout<<endl<<"#"<<label;
}
void set_program_configuration(){
	S=pow(RHO/(double)(N),(1.0/3.0));//sigma=(rho*V/N)**1/3;
	SS=S*S;	
}
string load_settings_file(){
	/* this function returns only the name of the file were the previous configuration will be readed*/
	char file_name[100];
	FILE *fp;
	fp = fopen("settings.dat","r");
	if (fp == NULL){//if cannot open the file
		print_label("The file that set the configuration of the simulation can not be read");
		return "ERROR: CANT READ THE FILE settings.dat";
	}else{
		print_label("File settings.dat opened succesfully");
	}
	char car;
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &N);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &T);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &RHO);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &lb);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &DISPL);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &NMOVE);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &NSUB);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &load_prev);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%s", &(*file_name));
	fclose(fp);
	string str(file_name);
	set_program_configuration();
	return str;
}

string load_settings_fileOneSite(){
	/* this function returns only the name of the file were the previous configuration will be readed*/
	char file_name[100];
	FILE *fp;
	fp = fopen("settings.dat","r");
	if (fp == NULL){//if cannot open the file
		print_label("The file that set the configuration of the simulation can not be read");
		return "ERROR: CANT READ THE FILE settings.dat";
	}else{
		print_label("File settings.dat opened succesfully");
	}
	char car;
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &N);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &T);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &RHO);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &lb);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &DISPL);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &NMOVE);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &NSUB);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%d", &load_prev);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%s", &(*file_name));
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &DISPLONESITE);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &rc);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &rd);
	while((car=fgetc(fp))!=' '){
	}
	fscanf (fp, "%lf", &epsilon);
	fclose(fp);
//	printf("\n You oppened a previous configuration and the value of DISPLONESITE is: %lf",DISPLONESITE);
	string str(file_name);
	set_program_configuration();
	return str;
}

int set_initial_conf_file(string filename){
	/* Set the initial configuration for a system with N particles in a cubic box
	*  simulation, the length of the box is equal to 1 */
	print_label("Initial configuration from file: " + filename);
	FILE *fp;
	fp = fopen(filename.c_str(),"r");
	if (fp == NULL){//if cannot open the file
		print_label("The file can not be read");
		return 1;
	}
	// omit the two first lines of the file
	char car;
	while((car=fgetc(fp))!='\n'){//Omit the header of the file (the first line is the  header)
	}
	/*while((car=fgetc(fp))!='\n'){//Omit the header of the file (the first line is the  header)
	}*/
	for(int n=0;n<N;n++){ // setting the position of each particle
		fscanf(fp,"%lf",&x[n]); 
		fscanf(fp,"%lf",&y[n]);
		fscanf(fp,"%lf",&z[n]);
	}
	fclose(fp);
}
int set_initial_conf_fileOneSite(string filename){
	/* Set the initial configuration for a system with N particles in a cubic box
	*  each particle has a site, the coordinates are stored in the representation:
	*  xCoorPart yCoorPart zCoorPart xCoorSite yCoorSite zCoorSite
	*/
	print_label("Initial configuration from file: " + filename);
	FILE *fp;
	fp = fopen(filename.c_str(),"r");
	if (fp == NULL){//if cannot open the file
		print_label("The file can not be read");
		return 1;
	}
	// omit the two first lines of the file
	char car;
	while((car=fgetc(fp))!='\n'){//Omit the header of the file (the first line is the  header)
	}
	/*while((car=fgetc(fp))!='\n'){//Omit the header of the file (the first line is the  header)
	}*/
	for(int n=0;n<N;n++){ // setting the position of each particle
		fscanf(fp,"%lf",&x[n]); 
		fscanf(fp,"%lf",&y[n]);
		fscanf(fp,"%lf",&z[n]);
		fscanf(fp,"%lf",&xs[n]); 
		fscanf(fp,"%lf",&ys[n]);
		fscanf(fp,"%lf",&zs[n]);
		fscanf(fp,"%d",&label[n]);
	}
	fclose(fp);
}
void set_initial_conf_array(){
	/* Set the initial configuration for a system with N particles in a cubic box
	*  simulation, the length of the box is equal to 1 */
	// max particles per side in the simulation box
	int part_side = (int)pow(N,1.0/3.0);
	// make sure that the particles for side are enough to complete N particles
	if ((double)part_side < pow(N,1.0/3.0)){
		part_side += 1; 
	}
	double delta_pos = 1.0/part_side; // diference between adjacent particles
	double xi, ti, zi; // coords of each particle
	int n=0;
	print_label("Initial configuration from a predefined array");
	do{ // setting the position of each particle
		x[n] = delta_pos/2.0 + (double)(n%part_side)*delta_pos; 
		y[n] = delta_pos/2.0 + (double)((n/part_side)%part_side)*delta_pos;
		z[n] = delta_pos/2.0 + (double)((n/part_side/part_side)%part_side)*delta_pos;
		//printf("\n%d\t%E\t%E\t%E",n+1,x[n],y[n],z[n]);
		n++;
	}while(n<N);
}
void set_initial_conf_arrayOneSite(){
	/* Set the initial configuration for a system with N particles in a cubic box
	*  simulation, the length of the box is equal to 1 */
	// max particles per side in the simulation box
	int part_side = (int)pow(N,1.0/3.0);
	// make sure that the particles for side are enough to complete N particles
	if ((double)part_side < pow(N,1.0/3.0)){
		part_side += 1; 
	}
	double delta_pos = 1.0/part_side; // diference between adjacent particles
	double xi, ti, zi; // coords of each particle
	int n=0;
	print_label("Initial configuration from a predefined array");
	do{ // setting the position of each particle
		x[n] = delta_pos/2.0 + (double)(n%part_side)*delta_pos; 
		y[n] = delta_pos/2.0 + (double)((n/part_side)%part_side)*delta_pos;
		z[n] = delta_pos/2.0 + (double)((n/part_side/part_side)%part_side)*delta_pos;
		xs[n] = 1.0;	// Set the site alingned with its corresponding sphere
		ys[n] = 0.0;
		zs[n]=0.0;
		label[n] = -1; // the -1 indicates the site is not overlaped to anything else
		//printf("\n%d\t%E\t%E\t%E",n+1,x[n],y[n],z[n]);
		n++;
	}while(n<N);
}
void save_conf(string file_name){
	int n=0;
	FILE *fp;
	fp = fopen(file_name.c_str(),"w");
	print_label("Saving configuration file: " + file_name );
	fprintf(fp,"# x\t\ty\t\tz\t configuration saved for N: %d",N);
	do{ // setting the position of each particle
		fprintf(fp,"\n%E\t%E\t%E",x[n],y[n],z[n]);
		n++;
	}while(n<N);
	fclose(fp);
}
void save_confOneSite(string file_name,string msg){
	int n=0;
	FILE *fp;
	fp = fopen(file_name.c_str(),"w");
	if(msg!=""){
		print_label("Saving configuration file: " + file_name );
	}
	fprintf(fp,"# x\t\ty\t\tz\t configuration saved for N: %d",N);
	do{ // setting the position of each particle
		
		fprintf(fp,"\n%E\t%E\t%E\t%E\t%E\t%E\t%d",x[n],y[n],z[n],xs[n],ys[n],zs[n],label[n]);
		n++;
	}while(n<N);
	fclose(fp);
}
/****************************************************
*				Functions of MC scheme    			*
*****************************************************/
void boundary_conditions(double& x, double& y, double&z){
	if(x>1.0){
		x = x-1.0;
	}else if(x<0.0){
		x = x+1.0;
	}
	if(y>1.0){
		y = y-1.0;
	}else if(y<0.0){
		y = y+1.0;
	}
	if(z>1.0){
		z = z-1.0;
	}else if(z<0.0){
		z = z+1.0;
	}
}
void move_particle(int n){
	// The movement of each particle is in sigmas units
	double displx,disply,displz;
	double newx,newy,newz;
	displx=DISPL*(2.0*ran2(pseed)-1.0);
	disply=DISPL*(2.0*ran2(pseed)-1.0);
	displz=DISPL*(2.0*ran2(pseed)-1.0);
	newx = x[n]+displx;
	newy = y[n]+disply;
	newz = z[n]+displz;
	boundary_conditions(newx,newy,newz);
	x[n]=newx;
	y[n]=newy;
	z[n]=newz;
}

int DoSiteOverlap(int n,int step){
	/** This function check if site of parDoSiteOverlapicle n overlaps with another
		Care must be taken in this function, when the site n moves several scenarios are possible
		1.- The site n was not overlaped in the previous step, after movement
			* The site n overlaps with a not overlaped site
			* The site n overlaps with an overlaped site
			* The site n does not overlap any other site
		2.- The site n was previously overlaped in the previous step, after movement.
			* The site n overlaps a not overlaped site
			* The site n overlaps an overlaped site
			* The site n does not overlap any other site
	 **/
	
	double rr;
	label_i_changed = -1;
	label_j_changed = -1;
	if(label[n]==-1){// site n was not previosuly overlapped
		for(int i=0;i<N;i++){
			if(n!=i){
				rr = distance2OneSite(i,n);
				if(rr<RC2 && label[i]==-1){ // *The site n overlaps with a not overlaped site
					previous_val_label_i = label[i];
					label_i_changed = i;
					label[n] = i; // Change label to indicate this sites are overlaped
					label[i] = n;
					/*if(step>255867){
						printf("\n1: Site %d overlaps with site %d, therefore the labels corresponding are %d %d ",n,i,label[n],label[label[n]]);
					}//*/
					
// 					getchar();
 					return 0;
				}
			}
		}
	}else{// site n was previously overlapped
		rr = distance2OneSite(n,label[n]);
		if(rr<RC2){// if the site n remains overlaped with is previously overlapped site then do nothing
 			return 0;
		}
		/*if(step>255867){
	 		printf("\n Site %d and site %d has the labels %d %d",n,label[n],label[n],label[label[n]]);
		}//*/
		label_i_changed = label[n];
		previous_val_label_i = label[label[n]];
		label[label[n]] = -1; // since the distance among site n and site label[n] is greater than RC, then site label[n] now is not overlaped
		
		/*if(step>255867)
	 		printf("\n22: Site %d left the overlap with %d, therefore the labels corresponding are ",n,label[n]);
		//*/
		label[n]=-1;
		/*if(step>255867)
			printf("%d %d deben ser -1 -1",label[n],label[previous_val_label_i]);
		//*/
// 		getchar();
		//label[n]=-1; // asumme in the movement site n change from to be overlaped with site label[n] to be not overlaped to any other site
		/** Now check if site n do end overlaped to a different site **/
		for(int j=0;j<N;j++){
			if(n!=j and j!=label[n]){
				rr = distance2OneSite(j,n);
				if(rr<RC2  and label[j]==-1 ){ // *The site n overlaps with a not overlaped site
					previous_val_label_j = label[j];
					label_j_changed = j;
					label[n] = j; // Change label to indicate this sites are overlaped
					label[j] = n;
					/*if(step>255867)
						printf("\n3: Site %d overlaps with site %d, therefore the labels corresponding are %d %d ",n,j,label[n],label[j]);
					//*/
// 					getchar();
					return 0;
				}
			}
		}	
	}
}
void move_particleOneSite(int n,int moveSiteAndParticle,int step){
	/*** This function moves particle an rotate site, the two extra argunments are used to save the
	 * previous state of the configuration, it because if the movements is rejected the configuration
	 * must return to the initial state.
	 * if moveSiteAndParticle =0 then only move the particle but not the site
	 * ***/
	double displx,disply,displz;
	double newx,newy,newz;
	displx=DISPL*(2.0*ran2(pseed)-1.0);
	disply=DISPL*(2.0*ran2(pseed)-1.0);
	displz=DISPL*(2.0*ran2(pseed)-1.0);
	newx = x[n]+displx;
	newy = y[n]+disply;
	newz = z[n]+displz;
	boundary_conditions(newx,newy,newz);
	x[n]=newx;
	y[n]=newy;
	z[n]=newz;
	// Now rotate the particle
	if(moveSiteAndParticle=='1'){
		newx = xs[n]+DISPLONESITE*(2.0*ran2(pseed)-1.0); //DISPLONESITE is like an angle
		newy = ys[n]+DISPLONESITE*(2.0*ran2(pseed)-1.0);
		newz = zs[n]+DISPLONESITE*(2.0*ran2(pseed)-1.0);
		double norma = sqrt(newx*newx+newy*newy+newz*newz);
		xs[n] = newx/norma;
		ys[n] = newy/norma;
		zs[n] = newz/norma;
		/* Check if site n overlaps with another site*/
		DoSiteOverlap(n,step); // check that wituation happens, and upgrade the labels coordinate
	}
}
double distance2(int i, int j){
	//return the distance between particle i and j to the square
	double xx,yy,zz;
	xx = x[i]-x[j];
	yy = y[i]-y[j];
	zz = z[i]-z[j];
	if(xx>0.5){
		xx=1.0-xx;
	}else if(xx<-0.5){
		xx=1.0+xx;
	}
	if(yy>0.5){
		yy=1.0-yy;
	}else if(yy<-0.5){
		yy=1.0+yy;
	}
	if(zz>0.5){
		zz=1.0-zz;
	}else if(zz<-0.5){
		zz=1.0+zz;
	}
	xx = xx*xx;
	yy = yy*yy;
	zz = zz*zz;
	return (xx+yy+zz);
}
double distance2OneSite(int i, int j){
	//return the distance between site of particle i and j to the square
	double xx,yy,zz,xsi,xsj,ysi,ysj,zsi,zsj;
	xsi = x[i]+xs[i]*DIST; // DIST is given in sigmas
	xsj = x[j]+xs[j]*DIST;
	ysi = y[i]+ys[i]*DIST;
	ysj = y[j]+ys[j]*DIST;
	zsi = z[i]+zs[i]*DIST;
	zsj = z[j]+zs[j]*DIST;
	xx = xsi-xsj;
	yy = ysi-ysj;
	zz = zsi-zsj;
	if(xx>0.5){
		xx=1.0-xx;
	}else if(xx<-0.5){
		xx=1.0+xx;
	}
	if(yy>0.5){
		yy=1.0-yy;
	}else if(yy<-0.5){
		yy=1.0+yy;
	}
	if(zz>0.5){
		zz=1.0-zz;
	}else if(zz<-0.5){
		zz=1.0+zz;
	}
	xx = xx*xx;
	yy = yy*yy;
	zz = zz*zz;
	return (xx+yy+zz);
}
double energy_due_particle(int n){
	double rr,energy=0.0;
	for(int i=0;i<N;i++){
		if(n!=i){
			rr = distance2(i,n);
			energy += potential_energy(rr);
		}
	}
	return energy;
}


int checkAnotherOverlap(int n,int j){
	// this function verifies if the site n has an overlap without taking into account the site j
	double rr;
	for(int i=0;i<N;i++){
		if(n!=i && j!=i){
			rr = distance2OneSite(i,n);	//Potential energy due to site-site interaction
			if(rr<RC2){
				return 1; //the site n has an overlap with site i and maybe with j
			}
			
		}
	}
	return 0;
}
double energy_due_particleOneSiteOLD(int n){
	if(label[n]!=-1){ 
		return potential_energyOneSite(0.0); // if label != -1 it is overlapped
	}
	return 0.0; // if label==-1 it is not overlpped
}
double energy_due_particleOneSiteNEW(int n){
	if(label[n]!=-1){ 
		return potential_energyOneSite(0.0); // if label != -1 it is overlapped
	}
	return 0.0; // if label==-1 it is not overlpped
}
double energy_due_particleOneSite1(int n){
	double rr,energy=0.0;
	int part=-1;
	for(int i=0;i<N;i++){
		if(n!=i){
			//rr = distance2(i,n);	//Potential energy due to main interaction
			//energy += potential_energy(rr);
			rr = distance2OneSite(i,n);	//Potential energy due to site-site interaction
			energy += potential_energyOneSite(rr);
			if(energy==-1.0){
				part = i;
			}
			if(energy<-1.0){
			//	printf("\nOverlap with particle %d and %d of particle %d",part,i,n);
				printf("\nsigma = %lf",rc*S);	
				printf("\nCoordinates of %d %lf %lf %lf %lf %lf %lf %lf",n,x[n],y[n],z[n],x[n]+xs[n]*DIST,y[n]+ys[n]*DIST,z[n]+zs[n]*DIST,sqrt(distance2OneSite(n,i)));
				printf("\nCoordinates of %d %lf %lf %lf %lf %lf %lf %lf",i,x[i],y[i],z[i],x[i]+xs[i]*DIST,y[i]+ys[i]*DIST,z[i]+zs[i]*DIST,sqrt(distance2OneSite(i,part)));
				printf("\nCoordinates of %d %lf %lf %lf %lf %lf %lf %lf",part,x[part],y[part],z[part],x[part]+xs[part]*DIST,y[part]+ys[part]*DIST,z[part]+zs[part]*DIST,sqrt(distance2OneSite(part,n)));
				//getchar();
				// La cuestion que uede pasar es que un sitio i se acerque a un m y disminuya la energia, de igual manera un sitio j se acerca a este mismo sitio m pero sin traslapar con el sitio
				// ientonces el sitio m tiene un doble traslape y llego a esa situacion sin que este se tuviera que mover m, por lo tanto para mover una sitio se tiene que revisar que
				// si se se traslapa con otro nuevo r, este r no tenga ya algun traslape previo.
				return 10000000000000000.0;
			}
		}
	}
	return energy;
}

double total_energy(){
	double energy=0.0,rr;
	for(int i=0; i<(N-1);i++){
		for(int j=i+1;j<N;j++){
			rr = distance2(i,j);
			energy += potential_energy(rr);
		}
	}
	return energy;
}
double total_energyOneSite(){
	/* Compute the initial energy to the overlap of sites
	 *
	 * This function runs over the "coordinate" label[i], if it is equal to -1 then it is not overlaped
	 * but if label[i]!=-1 then it is overlaped and counts for potential energy of the system.
	 */
	double energy=0.0,rr;
	for(int i=0; i<N;i++){
		if(label[i]!=-1){
			energy += potential_energyOneSite(0.0);
		}
	}
	return energy/2.0;
}
double slaterSum(double r_ij){
	double xi = sqrt2PI*(r_ij/S-1.0)/lb;
	/*if(r_ij<S){
		return 0.000000000000001;
	}*/
	//printf("\nXI %E dist %E",xi,1.0-exp(-xi*xi));
	return (1.0-exp(-xi*xi));
}
double energyDueParticlePotential(int n){
	// it only work for hard sphere potential
	double rr;
	for(int i=0;i<N;i++){
		
		if(n!=i){
			rr = distance2(i,n);
			if(rr<SS){
				//printf("\nCollision");
				return 10000.0;
			}//*/
		}
	}
	//printf("\n\nprod %E \t log(prod) %lf\n",prod,-log(prod));
	//getchar();
	return 0;
}
/*double freeEnergyDueParticleOneSite(int n){
	// 16/april/16 I dont know how to modify this to take into account the quantum effects according to Scheraga scheme
 	double rr,prod=1.0;
	for(int i=0;i<N;i++){
		
		if(n!=i){
			rr = distance2(i,n);
			if(rr<SS){
				//printf("\nCollision");
				return 10000.0;
			}//
			prod *= slaterSum(sqrt(rr));
		}
	}
	//printf("\n\nprod %E \t log(prod) %lf\n",prod,-log(prod));
	//getchar();
	return -log(prod);
}*/
double freeEnergyDueParticleCavityFunction(int n){
	double rr,prod=1.0;
	int i=0;
	if(n<=1){ // the particles that do not interact are the ones with index 0 and 1
		i = 2;// if the moving particle is one of these, then do not take into account the
	}		  // interaction between them	
	for(;i<N;i++){
		
		if(n!=i){
			rr = distance2(i,n);
			if(rr<SS){
				//printf("\nCollision");
				return 10000.0;
			}//*/
			prod *= slaterSum(sqrt(rr));
		}
	}
	//printf("\n\nprod %E \t log(prod) %lf\n",prod,-log(prod));
	//getchar();
	return -log(prod);
}
double countMonomers(){
	double cnt=0.0;
	for(int n=0;n<N;n++){
		if(label[n]==-1){
			cnt++;
		}
	}
	return cnt;
}



/****************************************************
*			  Automatization functions   			*
*****************************************************/
double betterDISPL(double actualDISPL,int totalMovs,double accMovs){
	/* This code will choose the better DISPL
	 * 
	 * The base of this code is the bisection method, here it is assumed the box
	 * is of length 1.0, therefore the maximum DISPL must be 1.0 and the minimum 0.0
	 */
	static double minimum = 0.0, maximum = 1.0;
	static double prevMovs = 0.0, prevAccMovs = 0.0;
	double currentAcc,partialMovs,partialAccMovs;
	partialMovs = totalMovs-prevMovs;
	partialAccMovs = accMovs-prevAccMovs;
	currentAcc = partialAccMovs/partialMovs;
	if(currentAcc>0.4){
		minimum = actualDISPL;
	}else{
		maximum = actualDISPL;
	}
	prevMovs = totalMovs;
	prevAccMovs = accMovs;
	//printf("\n Partial acceptance %lf",currentAcc);
	return (minimum+maximum)/2.0;
}
double betterDISPLOneSite(double actualDISPL,int totalMovs,double accMovs,double accep){
	/* This code will choose the better DISPL
	 * 
	 * The base of this code is the bisection method, here it is assumed the box
	 * is of length 1.0, therefore the maximum DISPL must be 1.0 and the minimum 0.0
	 */
	static double minimum = 0.0, maximum = 1.0;
	static double prevMovs = 0.0, prevAccMovs = 0.0;
	double currentAcc,partialMovs,partialAccMovs;
	partialMovs = totalMovs-prevMovs;
	partialAccMovs = accMovs-prevAccMovs;
	currentAcc = partialAccMovs/partialMovs;
	if(currentAcc>accep){
		minimum = actualDISPL;
	}else{
		maximum = actualDISPL;
	}
	prevMovs = totalMovs;
	prevAccMovs = accMovs;
	//printf("\n Partial acceptance %lf",currentAcc);
	return (minimum+maximum)/2.0;
}
double betterDISPLOneSite1(double actualDISPL,int totalMovs,double accMovs,double accep){
	/* This code will choose the better DISPL
	 * 
	 * The base of this code is the bisection method, here it is assumed the box
	 * is of length 1.0, therefore the maximum DISPL must be 1.0 and the minimum 0.0
	 */
	static double minimum = 0.0, maximum = 1.0;
	static double prevMovs = 0.0, prevAccMovs = 0.0;
	double currentAcc,partialMovs,partialAccMovs;
	partialMovs = totalMovs-prevMovs;
	partialAccMovs = accMovs-prevAccMovs;
	currentAcc = partialAccMovs/partialMovs;
	if(currentAcc>accep){
		minimum = actualDISPL;
	}else{
		maximum = actualDISPL;
	}
	prevMovs = totalMovs;
	prevAccMovs = accMovs;
	//printf("\n Partial acceptance %lf",currentAcc);
	return (minimum+maximum)/2.0;
}
void updateDISPL(double DISPL1){
	/* edit the file config.dat with the appropiate data for the simulation
	 * after it is going to be copied in the corresponding directory 
	*/
	print_label("Updating the value of DISPL to an optimized one in settings.dat");
	FILE *fp;
	fp=fopen("settings.dat","w");
	fprintf(fp,"Particles(N): %d",N);
	fprintf(fp,"\nT*: %lf",T);
	fprintf(fp,"\nRHO: %lf",RHO);
	fprintf(fp,"\nlb: %E",lb);
	fprintf(fp,"\nDISPL: %E",DISPL1);
	fprintf(fp,"\nNMOVE: %d",NMOVE);
	fprintf(fp,"\nNSUB: %d",NSUB);
	fprintf(fp,"\nLoad_prev_conf: %d",1);
	fprintf(fp,"\nSource_file: %s",filename_prev_conf.c_str());
	fclose(fp);
}
void updateDISPLOneSite(double DISPL1,double DISPL1ONESITE){
	/* edit the file config.dat with the appropiate data for the simulation
	 * after it is going to be copied in the corresponding directory 
	*/
	print_label("Updating the value of DISPL to an optimized one in settings.dat");
	FILE *fp;
	fp=fopen("settings.dat","w");
	fprintf(fp,"Particles(N): %d",N);
	fprintf(fp,"\nT*: %lf",T);
	fprintf(fp,"\nRHO: %lf",RHO);
	fprintf(fp,"\nlb: %E",lb);
	fprintf(fp,"\nDISPL: %E",DISPL1);
	fprintf(fp,"\nNMOVE: %d",NMOVE);
	fprintf(fp,"\nNSUB: %d",NSUB);
	fprintf(fp,"\nLoad_prev_conf: %d",1);
	fprintf(fp,"\nSource_file: %s",filename_prev_conf.c_str());
	fprintf(fp,"\nDISPL: %E",DISPL1ONESITE);
	fprintf(fp,"\nrc: %E",rc);
	fprintf(fp,"\nrd: %E",rd);
	fprintf(fp,"\nepsilon: %E",epsilon);
	fclose(fp);
}
void updateNMOVE(char *MOVES_1){
	/* edit the file config.dat with the appropiate data for the simulation
	 * after it is going to be copied in the corresponding directory 
	*/
	print_label("Updating the value of DISPL to an optimized one in settings.dat");
	FILE *fp;
	fp=fopen("settings.dat","w");
	fprintf(fp,"Particles(N): %d",N);
	fprintf(fp,"\nT*: %lf",T);
	fprintf(fp,"\nRHO: %lf",RHO);
	fprintf(fp,"\nlb: %E",lb);
	fprintf(fp,"\nDISPL: %E",DISPL);
	fprintf(fp,"\nNMOVE: %s",MOVES_1);
	fprintf(fp,"\nNSUB: %d",NSUB);
	fprintf(fp,"\nLoad_prev_conf: %d",1);
	fprintf(fp,"\nSource_file: %s",filename_prev_conf.c_str());
	fclose(fp);
}
void updateNMOVEOneSite(char *MOVES_1){
	/* edit the file config.dat with the appropiate data for the simulation
	 * after it is going to be copied in the corresponding directory 
	*/
	print_label("Updating the value of DISPL to an optimized one in settings.dat");
	FILE *fp;
	fp=fopen("settings.dat","w");
	fprintf(fp,"Particles(N): %d",N);
	fprintf(fp,"\nT*: %lf",T);
	fprintf(fp,"\nRHO: %lf",RHO);
	fprintf(fp,"\nlb: %E",lb);
	fprintf(fp,"\nDISPL: %E",DISPL);
	fprintf(fp,"\nNMOVE: %s",MOVES_1);
	fprintf(fp,"\nNSUB: %d",NSUB);
	fprintf(fp,"\nLoad_prev_conf: %d",1);
	fprintf(fp,"\nSource_file: %s",filename_prev_conf.c_str());
	fprintf(fp,"\nDISPLONESITE: %E",DISPLONESITE);
	fprintf(fp,"\nrc: %E",rc);
	fprintf(fp,"\nrd: %E",rd);
	fprintf(fp,"\nepsilon: %E",epsilon);
	fclose(fp);
}
/****************************************************
*			  		other functions 	  			*
*****************************************************/
struct tm *time_info;
double seconds;
time_t start_t,end_t; //saves the initial and final time of the program 


/****************************************************
*		  		definition of classes 	  			*
*****************************************************/
/*#define MAX_GR_BINS 200
class gr{
	private:
		int N; // the number of particles to compute the radial distribution function
		double timeComputed = 0.0;
		double *x;
		double *y;
		double *z;
		double gr_bin[MAX_GR_BINS];
		double maximum; // the maximum distance over which the gr will be computed
		double delta_r; // size of the bin.
		//double grFactor;
		string fileName="default_gr.dat";
	public:
		// member fucntions
		void compute_gr(); 
		void save_gr();
}
void gr(n,double& xcoords,double& ycoords,double& zcoords,double max,string outputFile){
	// Constructor for object of type gr
	N = n;
	maximum = max;
	delta_r = maximum/(double)MAX_GR_BINS;
	//grFactor = 1.0/(4.0*PI*RHO*delta_r*delta_r*delta_r);
	x = &xcoords; // pointing to the coordinates of each particle on the system
	y = &ycoords;
	z = &zcoords;
	fileName = outputFile;
}
void gr::compute_gr(){
	// This function counts the number of neighbors around each one of the N particles, 
	// because it uses the function square root, then it must no be used so extensively
	double dist;
	int bin_index;
	for(int i=0;i<(N-1);i++){
		for(int j=i+1;j<N;j++){
			double dx,dy,dz;
			dx = x[i]-x[j];
			dy = y[i]-y[j];
			dz = z[i]-z[j];
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if(dist>0.5) //only consider distances no longer than the half length of the box
				continue;
			bin_index=(int)(dist/delta_r);
			gr_bin[bin_index]+=2;//there is a particle corresponding to the distance bin_index times delta_r	
		}
	}
}
void gr::save_gr(double number_measures){
	FILE *fp,*fp2;
	double cnt;
	double gr_factor,ri,r,grVal;
	fp=fopen(fileName.c_str(),"w");
	fprintf(fp,"#This file contains the radial distribution function for a system with %d particles, rho* %lf, temp %lf, LB %lf and xlam %lf, it is result of %1.1lf average measures",N,RHO,T,XLAM,timesComputed);
	for(int n=0;n<MAX_GR_BINS;n++){
		ri = n*delta_r;
		gr_factor = 1.0/(4.0*PI*RHO*delta_r*(ri*(ri+delta_r)+delta_r*delta_r/3.0));
		r = delta_r*((double)n+0.5);
		grVal = gr_factor*gr_bin[n]/(timesComputed*(double)N)
		fprintf(fp,"\n%lf\t%lf",r,grVal);
	}
	fclose(fp);	
}	
/****************************************************
*	Definition of classical radial distribution fucntion 	  			*
*****************************************************/
#define MAX_GR_BINS 200
double delta_r=0.5/MAX_GR_BINS;
int grBin[MAX_GR_BINS];
double fact_gr;	
void compute_gr(){
	// This function counts the number of neighbors around each one of the N particles, 
	// because it uses the function square root, then it must no be used so extensively
	double dist;
	int binIndex;
	for(int i=0;i<(N-1);i++){
		for(int j=i+1;j<N;j++){
			dist=sqrt(distance2(i,j)); // compute the distance between particle i and j
			if(dist>0.5) //only consider distances no longer than the half length of the box
				continue;
			binIndex=(int)(dist/delta_r); // Count into the corresponding bin.
			grBin[binIndex]+=2;//there is a particle corresponding to the distance bin_index times delta_r	
		}
	}
}
void computeCavity_gr(){
	double dist;
	int binIndex;
	dist=sqrt(distance2(0,1)); // compute the distance between particle i and j
	binIndex=(int)(dist/delta_r); // Count into the corresponding bin.
	grBin[binIndex]+=2;//there is a particle corresponding to the distance bin_index times delta_r		
}
void saveCavity_gr(double numberMeasures){
	FILE *fp;
	double cnt;
	fp=fopen("radialCavity.dat","w");
	//fprintf(fp,"#This file contains the radial distribution function for a system with %d particles, rho* %lf, temp %lf, LB %lf and xlam %lf, it is result of %1.1lf average measures",N,RHO,TEMP,L_B,XLAM,number_measures	);
	for(int n=0;n<MAX_GR_BINS;n++){
		cnt=(double)(n*n+n)+ 0.33333333; // it is like r^2 a normalization factor for the g(r)
		fprintf(fp,"%lf\t%lf\n",delta_r*((double)n+0.5)/S,fact_gr*(double)grBin[n]/(2.0*cnt*numberMeasures));
	}
	fclose(fp);	
}
void save_gr(double numberMeasures){
	FILE *fp;
	double cnt;
	fp=fopen("radial.dat","w");
	//fprintf(fp,"#This file contains the radial distribution function for a system with %d particles, rho* %lf, temp %lf, LB %lf and xlam %lf, it is result of %1.1lf average measures",N,RHO,TEMP,L_B,XLAM,number_measures	);
	for(int n=0;n<MAX_GR_BINS;n++){
		cnt=(double)(n*n+n)+ 0.33333333; // it is like r^2 a normalization factor for the g(r)
		fprintf(fp,"%lf\t%lf\n",delta_r*((double)n+0.5)/S,fact_gr*(double)grBin[n]/(cnt*numberMeasures*(double)N));
	}
	fclose(fp);	
}
void saveStatus(int done, double energy,double seconds){
	FILE *fp;
	double cnt;
	fp=fopen("status.dat","w");
	//fprintf(fp,"#This file contains the radial distribution function for a system with %d particles, rho* %lf, temp %lf, LB %lf and xlam %lf, it is result of %1.1lf average measures",N,RHO,TEMP,L_B,XLAM,number_measures	);
	fprintf(fp,"#done\tEnergy\tSeconds\n");
	fprintf(fp,"%d\t%lf\t%lf",done,energy,seconds);
	fclose(fp);	
}
