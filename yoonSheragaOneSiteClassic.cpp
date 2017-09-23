/* Author: Cesar Serna
 * Date: 31/oct/2015
 *
 * the coordinates of each particle are 6, 3 of the main particle and other 3 of the site
 * in this case it is used another extra "coordinate" which label if a site is overlaped 
 * with another site, then, it has to be take care of update this label, when this overlap
 * and when it goes from be overlaped to be free.
 */
#include "rutinasOneSite.h"
#include <iostream>
using namespace std;
/****************************************************
*		Declaration of the potential energy		*
*****************************************************/
void checkEverthingOK(int step){
	for(int i=0;i<N;i++){
		if(label[i]!=-1){
			if(i!=label[label[i]]){
				printf("\nAqui ocurrio un error, en el paso %d",step);
				printf("\n particle i: %d label[label[i]]: %d ",i,label[label[i]]);
				//printf("label[]")
				save_confOneSite("final_conf.dat","");
				getchar();
			}
		}
	}
}
void printConfiguration(){
	for(int i=0;i<N;i++){
		printf("\n%lf\t%lf\t%lf\t%d",x[i],y[i],z[i],label[i]);
	}
	printf("\n");
}
double potential_energy(double rr){
	//Here is the definition of the interation potential
	double energy=0.0;
	if(rr<SS){ // Inside the core of the hard particle
		energy = 1000000.0;
	}/*else if(rr<=XLAMS2){
		energy-=1.0;
	}//*/
	/*if(rr>=16*SS)//The potential is truncated to 4sigmas
		return (0.0);
	energy = 4.0*(pow(SS/rr,6)-pow(SS/rr,3));//Lennard-Jones potential 
	//energy = energy + 0.00097632; // potential shift*/
	return energy;
}
char nmovs[10];
double potential_energyOneSite(double rr){
	//Here is the definition of the interation potential
	double energy=0.0;
	//if(rr<0.5326*0.5326*S*S){ // Inside the core of the hard particle
	//	energy -= 7;
	//if(rr<SS*0.2837){ // Inside the core of the hard particle
	if(rr<RC2){ // Inside the core of the hard particle
		energy -= EPSILON;
	}/*else if(rr<=XLAMS2){
		energy-=1.0;
	}//*/
	/*if(rr>=16*SS)//The potential is truncated to 4sigmas
		return (0.0);
	energy = 4.0*(pow(SS/rr,6)-pow(SS/rr,3));//Lennard-Jones potential 
	//energy = energy + 0.00097632; // potential shift*/
	return energy;
}
void printConf(){
	for(int n=0;n<N;n++){
		printf("\n%E\t%E\t%E\t%E\t%E\t%E\t%d",x[n],y[n],z[n],xs[n],ys[n],zs[n],label[n]);
	}
}
			
int main(int argc, char* argv[]){
	PI = acos(-1.0);
	sqrt2PI=sqrt(2.0*PI); // used in semiclassical simulation
	time(&start_t);//start counting time
	double oldx,oldy,oldz,oldxs,oldys,oldzs,oldE,newE,oldPotentialE,newPotentialE,movsAcc,prevMovsAcc=0,tot_E=0.0,nMeasures=0.0,perc=0.0;
	double currentAcc,partialAcc, energyCounter=0.0,prevEnergy=0.0,currentEnergy,freeEnergyCounter;
	int n;
    char moveSiteAndParticle = '0';// if it is zero then only move the particle to reach a 50% of acceptance
								// in order to find a value for DISPL and DISPLSITE for which the global acceptance
								// reach a value of 40%, 10^6 moves are done to find DISPL at 50% and 10^6 moves
								// are done and DISPLSITE are correcte to fit a 40% global accceptance
	pseed=&seed;
	seed=-123456789;
	filename_prev_conf = load_settings_fileOneSite(); // Load the settings of the simulation
	if (argc!=1){ // If have an argument in the simulation then update the settings.dat file
		// if an argument is received then it is assumed that DISPL and DISPLONESITE is already adjusted to
		// fit a global acceptance of 40%
		updateNMOVEOneSite(argv[1]);
		filename_prev_conf = load_settings_fileOneSite(); // Load the settings of the simulation
	}
	printf("\n# N:\t%d\n# T*:\t%E\n# RHO*\t%E\n# lb*:\t%E\n# DISPL:\t%E\n# maxdL:\t%E\n# NMOVE:\t%d\n# NSUB:\t%d\n# load prev config:\t%d\n#file name:\t%s\n#DISPLONESITE:\t%E\n#rc:\t%E\n#rd:\t%E\n#epsilon:\t%E",N,T,RHO,lb,DISPL,maxdS,NMOVE,NSUB,load_prev,filename_prev_conf.c_str(),DISPLONESITE,rc,rd,epsilon);	
	// Load initial configuration
	if (load_prev==1){
    		moveSiteAndParticle = '1';// if it is zero then only move the particle to reach a 50% of acceptance
		set_initial_conf_fileOneSite(filename_prev_conf);
	}else{
		set_initial_conf_arrayOneSite(); 
		DISPL = 0.5;
		DISPLONESITE = 0.5;
	}
	fact_gr=S*S*S/(4.0*PI*RHO*delta_r*delta_r*delta_r);// once the system variables are loaded the fact_gr is computed
	printf("\nSigma en unidades de la caja %lf",S);	   // this is used in computing the radial distribution function
	DIST = rd*S;	//Used for sites simulationsm defines the distance between the center of the site and thje center
	RC2 = rc*rc*S*S;
	EPSILON = epsilon;
	tot_E = total_energyOneSite();
	prevEnergy = tot_E;
	printf("\nInitial energy: %lf",tot_E/N);
	save_confOneSite("initial_conf.dat","");
	movsAcc=0.0;
	// Printing label current progress
	printf("\nDone\tEnergy\t\tPartial Acc\tAcceptance\tDISPL\t\tDISPONESITE\tSecs\tmean_x\tcurrent_X");
	double auxCnt1=0;
				// of the hard sphere, DIST define the distance in sigmas.
	double cntMonomer=0.0,currentMonomers=0.0;
	float seconds,temp;
	int oldLabel;
	label_j_changed = -1;
	label_i_changed = -1;
	int stepIndex = 255867;
	for(int i=0;i<NMOVE;i++){
		
		n = (int)(ran2(pseed)*(N)); // choose a particle randomly
		/*if(i>stepIndex){
			printf("\nThe particle moving is %d, old coordinates %lf %lf %lf %d",n,x[n],y[n],z[n],label[n]);
		}//*/
		oldx = x[n];//save old coordinates
		oldy = y[n];
		oldz = z[n];
		oldxs = xs[n]; // Save coordinates of the sites
		oldys = ys[n];
		oldzs = zs[n];
		oldLabel = label[n];
		oldE = energy_due_particleOneSiteOLD(n);
		oldPotentialE = energyDueParticlePotential(n)/(2.0);
		/*if(i>255867){
			printf("\nAlgo raro paso aquiiiii%d %lf",i,oldE);
			printf("\nCoordinates %d %lf %lf %lf",i,x[i],y[i],z[i]);
			printf("\nCoordinates of %d %lf %lf %lf %lf %lf %lf %lf",i,x[i],y[i],z[i],x[i]+xs[i]*DIST,y[i]+ys[i]*DIST,z[i]+zs[i]*DIST,oldE);
			temp = energy_due_particleOneSite1(n);
			getchar();
		}//*/
		
		//newE=energy_due_particleOneSite(n);			
		move_particleOneSite(n,moveSiteAndParticle,i);
		/*if(i>stepIndex){
			printf("\nThe particle moving is %d, NEW coordinates %lf %lf %lf %d",n,x[n],y[n],z[n],label[n]);
		}//*/
		newPotentialE = energyDueParticlePotential(n)/(2.0);			
		newE = energy_due_particleOneSiteNEW(n);			
		/*if(n==93 && newE>0){
			printf("\nAlgo raro paso aqui1 %d %lf",i,newE);
			temp = energy_due_particleOneSite1(n);
		}*/

		//printf("\nNEW %E\tOLD%E",newPotentialE,oldPotentialE);
		if((newE-oldE)+(newPotentialE-oldPotentialE)<0.0){
			energyCounter += newE-oldE; // count the increment of energy
			freeEnergyCounter += newPotentialE-oldPotentialE; // count the increment of free energy
			movsAcc++;
			/*if(i>255869){
				printf("\nacepted1 %lf %lf %d %lf %lf %d",newPotentialE-oldPotentialE, newE-oldE,i,newE,oldE,n);
				printf("\nLlamando a la funcio para ver que sitio se traslapo");
				temp = energy_due_particleOneSite1(n);
				getchar();
			}//*/
				
		}else{
			
			if(ran2(pseed)<exp(-(newPotentialE-oldPotentialE + newE-oldE)/T)){
				movsAcc++;
				freeEnergyCounter += newPotentialE-oldPotentialE; // count the increment of energy
				energyCounter += newE-oldE; // count the increment of energy
				/*if(i>255869){
					printf("\nacepted2 %lf %lf %d %lf %lf",newPotentialE-oldPotentialE, newE-oldE,i,newE,oldE);
					getchar();
				}///*/	
		//		printf("\nACCEPTED2");
			}else{
				/*if(i>=110444){
					printf("\nrejected %lf %lf %d %lf %lf",newPotentialE-oldPotentialE, newE-oldE,i,newE,oldE);
					getchar();
				}*/
// 				printf(" REJECTED");
				x[n]=oldx;
				y[n]=oldy;
				z[n]=oldz;
				xs[n] = oldxs;
				ys[n] = oldys;
				zs[n] = oldzs;
				label[n] = oldLabel;
				if(label_i_changed!=-1){
// 					printf(" i: label %d has the value %d and will be restored to the value %d",label_i_changed, label[label_i_changed], previous_val_label_i);
					label[label_i_changed] = previous_val_label_i;
				}
				
				if(label_j_changed!=-1){
// 					printf(" j: label %d has the value %d and will be restored to the value %d",label_j_changed, label[label_j_changed], previous_val_label_j);
					label[label_j_changed] = previous_val_label_j;
				}
// 				printf("\n-- Site %d have the label %d ",n,label[n]);
			}
		}
		/*if(i>stepIndex){
			printConfiguration();
		}
		checkEverthingOK(i);
		//*/
		//getchar();
		if(i%NSUB==0){
			// Taking averages
			currentEnergy = prevEnergy + energyCounter;
			//tot_E += total_energyOneSite();
			tot_E += currentEnergy;
			prevEnergy = currentEnergy;
			energyCounter=0.0;
			nMeasures++;
			currentMonomers = countMonomers(); 
			cntMonomer += currentMonomers;
			//compute_gr();
			//save_gr(nMeasures);
		}
		if((i+1)%(NMOVE/100)==0){
			perc++;
			if (perc>50){
				moveSiteAndParticle='1';
// 				getchar();
			}
			/*if (perc>=55)
				getchar();*/
			
			currentAcc = (double)movsAcc/(double)(i+1);
			time(&end_t);
			seconds=difftime(end_t,start_t); 
			partialAcc = (movsAcc-prevMovsAcc)/(i-auxCnt1);
			printf("\n%1.1f\t%E\t%E\t%E\t%E\t%E\t%1.2f\t%E\t%E",perc,tot_E/((double)N*nMeasures),partialAcc,currentAcc,DISPL,DISPLONESITE,seconds,cntMonomer/((double)N*nMeasures),currentMonomers/N);
			if((int)perc%5	== 0 && load_prev == 0 ){ 
				if (moveSiteAndParticle=='0'){
					DISPL = betterDISPLOneSite(DISPL,i+1,movsAcc,0.63);
				}else{
					DISPLONESITE = betterDISPLOneSite1(DISPLONESITE,i+1,movsAcc,0.4);
				}
			}//*/
			auxCnt1 = i;
			prevMovsAcc = movsAcc; 
			saveStatus(perc,tot_E/((double)N*nMeasures),seconds);
			save_confOneSite("final_conf.dat","");
		}
		
	}
	if(load_prev==0){
		updateDISPLOneSite(DISPL,DISPLONESITE); // change the settings.dat file, assume the next run of the program will use the 
	}						// optimized value of DISPL and change load_prev_conf flag value to 1 in case of a 
	cout<<endl;         	// pre-programmed run is scheduled	
	save_confOneSite("final_conf.dat","Salvando configuracion final");
	return 0;			
}
