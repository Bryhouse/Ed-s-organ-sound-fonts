
/*	ELA_synthesis_532: Variable Timbre Generation	 */
/*	This program was developed and written independently in 
	Feb 2007 by
	Edwin L. Althouse, elalthouse@comcast.net.	*/

/*	Modified 10/1/07 to use real (non-integer) values for the partial_indices  */

/*	Modified 10/1/07 to use phase offsets determined from curve fitting
	to the phase data, but using calculated frequencies instead of the continuous
	phase data   */

/*	*************************   */
/*	Modified 9/14/07 to scale timbre by imperical expression
	delta_y(dB) = -2.24*(H-1)*ln(ff/freq0), where delta_y(dB) is the difference in
	power between the spectral powers for each harmonic (H) component at frequency "f"
	relative to that specified at frequency = freq0.  The dB values computed will be
	converted to a numerical fraction that will be used as a multiplier to modify the
	envelope data provided at frequency freq0.		*/

/*	modified 9/30/07 to scale time axis by
	exp=log10(freq0/ff)/3.0;
	Fa=pow(10.0,exp);  i.e., Fa=cube root of (freq0/ff)	*/
/*	************************* */

/*	Modified 8 Sept 07 to extend mas waveform time to 1.5 sec  */

/*	This version includes a capability to shift frequencies sharp or flat in attack phase
	to create more interesting sound  */
/*	This version scales the chiff by fourth-root-of(freq0/ff)	*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void Read_myconfig532(char *fin);
void Read_f_offsets(char *fin);
void Read_notes(char *fin);
void Read_partial_indicies_real(char *fin);
void Read_time(char *fin);
void Read_envelope_names(void);
void Read_partial_indicies_real(char *fin);
void Write_output(char *fout);
void Modify_chiff(void);
void Interp_env_across_time(int k);
void Read_envelope(char fin);
void compute_signal(float freq,int ilast, int kk);
void Normalize_signal(void);
void Read_phase_offsets(char *fin2);

float interpolate2(float yyy,int nnn,int nn_end);
float maxftn(float a, float b, float c);

float pi = 3.14159271f,degtorad=.017453293;
float w, freq0,f0,freq1,scale_param,chiff_ref_freq,chiff_mult,maxtime;
float t[5000],ff;
float tmod[5000], shift[103], shft_ref;
float y[5000], yy[47000],AA[103],BB[103];
float output[47000];
float freq[55],timb_scale, dB_scale;
float temp2,temp3;
float time_offset_end, offset, offset_m,f_offset[103],partial_index[103];

char envelope_name[103][30],note[58][30];
char filename0[]="f_offsets.txt";
char filename1[]="myconfig532.txt";
char filename3[]="notes.txt";
char filename5[]="time.txt";
char filename6[]="envelope_names.txt";
char filename4[]="partial_indicies_real.txt";
char filename7[]="phase_offsets.txt";
char *fin,*fin2,*fout,*ptr5,*ptr6,*ptr9;

int n_offset_end,chiff_index,chiff_end,bypass;
int nsamp,maxpts,maxnotes=54,num_partials;
int start_note, stop_note;


/*	----------------------------------------------------------------------------------------------------	*/
/* 	Read_myconfig532 function definition	*/
void Read_myconfig532(char *fin)
{
	
	FILE *fptr1;
	if ((fptr1 = fopen(fin, "r")) == NULL) {
	printf("Cannot open myconfig532.txt for reading.\n");
	exit(1); }
	else {
		
		fscanf(fptr1,"%d",&nsamp);// number of samples in the envelope  array 
		fscanf(fptr1,"%d",&num_partials);// number of partial frequencies used in the synthesis
		fscanf(fptr1,"%f",&freq0); //freq at which envelopes were determined
		fscanf(fptr1,"%f",&freq1); //Reference freq for timbre conversion (near low end of notes)
		fscanf(fptr1,"%f",&scale_param); //sets slope or timbre conversion formula (=2.42 for Principal pipe)
		fscanf(fptr1,"%f",&maxtime); //max synthesis signal time; max is 1.5 sec
		fscanf(fptr1,"%d",&start_note); // first note to be synthesized, e.g. C2
		fscanf(fptr1,"%d",&stop_note); // last note to be synthesized, e.g. C7
		fscanf(fptr1,"%f",&time_offset_end); // time when freq-offset effect dissapears, e.g. 0.4 sec
		fscanf(fptr1,"%d",&bypass); // bypass =1 eliminates special chiff processing
		fscanf(fptr1,"%d",&chiff_index);//specifies the partial index that carries the chiff
		fscanf(fptr1,"%d",&chiff_end); // envelope sample point where chiff ends
		fscanf(fptr1,"%f",&chiff_ref_freq);//Reference frequency for chiff processing
		fscanf(fptr1,"%f",&chiff_mult);//Freq-independent multiplier for the chiff
		fclose(fptr1);
		
		printf("nsamp = %d\n", nsamp);
		printf("num_partials = %d\n", num_partials);
		printf("freq0 = %10.4f\n", freq0);
		printf("freq1 = %10.4f\n", freq1);
		printf("scale_param = %10.4f\n", scale_param);
		printf("maxtime = %5.3f\n", maxtime);
		printf("start_note = %d\n", start_note);
		printf("stop_note = %d\n", stop_note);
		printf("time_offset_end = %f\n", time_offset_end);
		printf("bypass = %d\n", bypass);
		printf("chiff_index = %d\n", chiff_index);
		printf("chiff_end= %d\n", chiff_end);
		printf("chiff_ref_freq = %5.3f\n", chiff_ref_freq);
		printf("chiff_mult = %E\n", chiff_mult);
		}
}

/*	---------------------------------------------------------------------	 */
void Read_f_offsets(char *fin)
{
	int i;
	FILE *fptr0;
	fptr0 = fopen(fin,"r");
	if (fptr0 == NULL) {
	printf("Cannot open f_offsets.txt for reading.\n");
	exit(0); }
	else
	{
		for (i=0;i<=num_partials-1;i++) fscanf(fptr0,"%f", &f_offset[i]);
		fclose(fptr0);

		printf("Some frequency offsets are:\n");
		 for (i=0;i<=3;i++) printf("f_offset [%d]= %f\n",i,&f_offset[i]);
	}
}

/* ------------------------------------------------------------------- */
//*	Read_notes function definition	*/
void Read_notes(char *fin)
{
	int p;
	FILE *fptr3;
	fptr3 = fopen(fin, "r");
	if (fptr3 == NULL) {
	printf("Cannot open notes.txt for reading.\n");
	exit(3); }
	else 
	{
		for (p=0;p<=maxnotes-1;p++) fscanf(fptr3, "%s", &note[p][0]);
		fclose(fptr3); 

		 printf("A few notes read are:\n");
		 for (p=0;p<=9;p++) printf("Note [%d]= %s\n",p,&note[p][0]); 
		
	}
				
}
				 
/*	-------------------------------------------------------------------	*/
/*	Read_time function definition	*/
void Read_time(char *fin)
{
	FILE *fptr5;
	int k;
	if ((fptr5 = fopen(fin, "r")) == NULL) {
	printf("Cannot open time.txt for reading.\n");
	exit(5); }
	else {
		
		for(k=0;k<=nsamp-1;k++) fscanf(fptr5,"%E",&t[k]);		
		fclose(fptr5);
//		printf("The time data are:\n");
//		for (k=0;k<=nsamp-1; k++) printf("t[%3d] = %7.2f\n",k,t[k]);
			
		}			
}		

/*	---------------------------------------------------------------------	*/
/* Compute_tmod_array() function defintion      */
void Compute_tmod_array(void)
{		
	int n;
	double Fa,exp;
	
	//old Fa was Fa = (0.05*sqrt(freq0/ff) + 0.03f)/0.080f;
	exp=log10(freq0/ff)/3.0;
	Fa=pow(10.0,exp);
	for (n = 0;n<= nsamp-2; n++) tmod[n] = t[n]*Fa;
	tmod[nsamp-1] =maxftn(tmod[nsamp-2]+.1f,2.0f,0.0f); /* force the last element of tmod[] to stay at the extreme end of 
	the time scale  */
}

/*------------------------------------------------------------------------*/
/* Modify_chiff function definition   */
void Modify_chiff(void)
{
	int j,jend;
	float temp;
	/*	Now we want to modify the amplitude of the chiff as the frequency changes relative
		to a reference frequency (chiff_ref_freq).  We use the scaling law
		new_amplitude = old_amplitude * (chiff_ref_freq/ff), where ff is the current note frequency.	*/
		
		jend=chiff_end;
		temp = pow(10.0,log10(chiff_ref_freq/ff)/4.0);
		for (j = 0; j<=jend; j++)  y[j] = y[j]*temp*chiff_mult;

}
/* ---------------------------------------------------------------------------------*/
/* Read_partial_indicies_real function definition    */
void Read_partial_indicies_real(char *fin)
{
	FILE *fptr6;	
	int i;
	if ((fptr6 = fopen(fin, "r")) == NULL) {
	printf("Cannot open partial_indicies_real.txt for reading.\n");
	exit(6); }
	else {
	
		for(i=0;i<=num_partials-1;i++) fscanf(fptr6,"%f",&partial_index[i]);	
		fclose(fptr6);
		

		}
}
/*	-----------------------------------------------------------------------	*/	
/*	Read_envelope_names function definition	*/
void Read_envelope_names(void)
{
	//printf("called Read_envelope_names \n");
	
	int k;
	fin2=filename6;
	FILE *fptr3;

	if ((fptr3 = fopen(fin2, "r")) == NULL) {
	printf("Cannot open envelope_names for reading.\n");
	exit(3); }
	else {
		
		for(k=0;k<=num_partials-1;k++) fscanf(fptr3,"%s",&envelope_name[k][0]);
			
		fclose(fptr3);

		printf("A few envelope_names are:\n");
		
		ptr5=&envelope_name[0][0];
		printf("%s\n",ptr5);
		if(num_partials>1)	
				{
				ptr6=&envelope_name[1][0];
				printf("%s\n",ptr6);
				}
			
		}
}
/* ----------------------------------------------------------------------- */

/*	Read_envelope function definition	*/
void Read_envelope(char *fin)
{
	//printf("called Read_envelope_file\n");

	FILE *fptr5;
	int k;
	if ((fptr5 = fopen(fin, "r")) == NULL) {
	printf("Cannot open envelope for reading.\n");
	exit(5); }
	else {
		
		for(k=0;k<=nsamp-1;k++) fscanf(fptr5,"%E",&y[k]);
			
		fclose(fptr5);
		
		//printf("A few values of y[] are:\n");
		//for (k=100;k<=105; k++) printf("y[%3d] = %E\n",k,y[k]);
		}
}
/*	----------------------------------------------------------------------	*/
/* Interp_env_across_time function definition  */
void Interp_env_across_time(int k)
{
	int mend,m,j;
	float t1,tm,tm1;
	double a,b,g;
	mend=nsamp-1;
	
	for (j=0;j<=maxpts-1;j++) 
	{t1=j/44100.0;
	
	if (t1<= 0.0) yy[j]=0.0;
								
				else for (m =0; m<= mend-1; m++)
				{	tm=tmod[m];
					tm1=tmod[m+1];
					while ((t1 > tm) && (t1 <= tm1)) 
					{
						a=(y[m+1]-y[m])/(tmod[m+1]-tmod[m]);
						b = (y[m]*tmod[m+1]-y[m+1]*tmod[m])/(tmod[m+1]-tmod[m]);
						g = a*t1+b;
						yy[j]=g;
						break;	
					}
				}
	}
}
/* -------------------------------------------------------------------- */
/*	compute_signal definition	*/
void compute_signal(float freq,int nn, int kk)
{ 
	int jj;
	float gg, offset,timb_scale, dB_scale,phase_offset;
	float index,index1;
	offset=0.0;
	index=partial_index[kk-1];
	index1=index-1;

	//	Scale partials with frequency to roll down timbre at higher freqs
	dB_scale = -scale_param*index1*log(freq/freq1);
	timb_scale=pow(10,dB_scale/20.0);

	for (jj=0;jj<=maxpts-1;jj++) 
	{		if (index*freq < 18000.) gg=1.0;
			else gg=0.0;
	
		phase_offset=BB[kk-1];

			phase_offset=phase_offset*shift[nn];
			output[jj]=output[jj] + timb_scale*gg*yy[jj]*sin(freq*index*w*jj+phase_offset);
		
	}
	//printf("output[2000]= %E\n", output[2000]);
	
}

/* -------------------------------------------------- */

/* ------------------------------------------------*/

/* Normalize_signal() function defintion  */
void Normalize_signal(void)
{
/*	Now normalize the maximum amplitude of the output[] array to 3E4.	*/
	int n;
	float vmax = 0.0f;
	float vmin = 0.0f;

	for (n=0;n<=maxpts-1;n++) if (output[n] > vmax) vmax = output[n];
	for (n=0;n<=maxpts-1;n++) if (output[n] < vmin) vmin = output[n];
	if (fabs(vmin) > vmax) vmax = fabs(vmin);

	printf("vmax = %e\n",vmax);
	for (n=0;n<=maxpts-1;n++) output[n] = 30000.0f*output[n]/vmax;
//	for (n=0;n<=maxpts-1; n=n+200) printf("output[%d] = %f\n", n,output[n]);
}
/* ------------------------------------------------------- */

/*	Defintion of Write_output()	*/
void Write_output(char *fout)
{

	FILE *fptr7;
	int k;
	printf("opening: %s\n", fout);
	if ((fptr7 = fopen(fout, "w")) == NULL) {
	printf("Cannot open output file for writing.\n");
	exit(7);} 
	else {
		for(k=0;k<=maxpts-1;k++) fprintf(fptr7,"%e\n",output[k]);			
		fclose(fptr7); }
}

/*	---------------------------------------------------------------  */
/*	Definition of function maxftn(a,b,c)	*/

float maxftn(float a, float b, float c)
{
float temp3;

	if (a>b) temp3 =a;
	else temp3 = b;
	if(c>temp3) temp3=c;

	return temp3;
}
/* ------------------------------------------------------- */
/* --- Definition of function interpolate2(yyy,nnn,nn_end) ---- */
float interpolate2(float yyy,int nnn,int nn_end)
{
	float temp4;
	double zz;
	zz=pi*nnn/(nn_end*2.0);	
	temp4=yyy*(1.0-sqrt(sin(zz)));
		
	return temp4;
}
/* ----------------------------------------------------------------------- */



/*	-----------------------------------------------------------------------	*/	
/* Read_phase_offsets(fin2) function definition  */
void Read_phase_offsets(char *fin2)
{
	printf("called phase_offsets.txt\n");

	FILE *fptr6;
	int k;
	if ((fptr6 = fopen(fin2, "r")) == NULL) {
	printf("Cannot open phase_offsets.txt for reading.\n");
	exit(6); }
	else {
		
		for(k=0;k<=num_partials-1;k++) fscanf(fptr6,"%E",&BB[k]);
			
		fclose(fptr6);
		
		printf("A few values of AA[] and BB[] are:\n");
		for (k=0;k<=num_partials-1; k++) printf("BB[%4d]  = %E\n",k,BB[k]);
		}
}	

/*	-----------------------------------------------------------------------	*/

/* ---------------Begin Main -----------------------------------*/
/* ------------------------------------------------------------- */	
/* ------------------------------------------------------------- */		
int main( void )
{
/*	nsamp is the number of time samples used to define each envelope and also for the 
time_profile.  They all must have the same number of samples and each envelope sample
must apply at the corresponding time_profile value.    */ 

/*	freq0 is the reference frequency for which the envelope was derived.
e.g., 130.8 Hz if the envelopes were derived from a C3 note	*/

/*	maxtime is the maximum amount of waveform time sample to calculate 
	(not to exceed 1 sec).  */

/*	maxnotes is the maximum number of keygroups to generate.  Must be <=54.	*/
		
/*	Input description of note keygroups (2 keys), starting with a 16' pipe at lowest key 
	(e.g., C1, D1, E1, F#1, G#1, A#1, ..... and ending at A#9.  This is 9 octaves or 
	108/2 = 54 keys since I will be generating only every other note	*/


int i,k,n,r,ilast;

w=2.0*pi/44100.0;
freq[1] = 32.7031;

/* 	Input configuation information via "myconfig532.txt"	*/
fin=filename1;
Read_myconfig532(fin);

maxpts=maxtime*44100-1;

fin=filename3;
Read_notes(fin);

fin=filename4;
Read_partial_indicies_real(fin); /* 1 for fundamental, 2 for 2nd harmonic, etc  */

/*	Read in the times associated with each point of the ADSR envelope profile */
fin = filename5;
Read_time(fin);

/*  Read in ADSR envelope names.     */
Read_envelope_names();

/*	Read in Phase names.	*/
fin2=filename7;
Read_phase_offsets(fin2);

/*	Compute shift factor for multiplying phase data */
f0=32.70317;
shft_ref=1.122462048;
shift[1] = f0/freq0;
for (k=2;k<=maxnotes;k++) shift[k]=shift[k-1]*shft_ref;  

/*	Calculate the approximate note frequencies (in keygroups of two) starting at C1 (= 32.7031 Hz)	*/
for (r = 2; r<= maxnotes; r++) freq[r] = freq[r-1]*1.122462048f;

for (i =start_note; i<= stop_note ; i++) 
{ // start notes or "i" loop 

	ff=freq[i];
	ilast=i;

/*	Compute the tmod[n] array, which corresponds to the time-shifted t[n] array where n 
	references increasing ADSR-profile samples (1 to nsamp).  The time shift must be done
	to account for the change in risetime of the pipe waveform as the frequency changes.
	The shift is toward the origin of the time axis as the frequency increases.  The following
	imperical law is currently being used:  
	Scaling factor = Fa = (0.25*sqrt(freq0/ff) + .03)/0.280  */
Compute_tmod_array();		


		for (n=0;n<=maxpts-1;n++) output[n]=0.0; // initialize output array for each new note
		
		for (k=1; k<= num_partials; k++) 		
		{// start partials or "k" loop
			
			fin=&envelope_name[k-1][0];
			Read_envelope(fin); // builds y[] array for current partial
			
			if ((bypass!=1)&&(k==chiff_index)) Modify_chiff(); 
			// modifies y[] array containing chiff
							
			Interp_env_across_time(k); // Builds yy[] array at 44100 sample rate.
					
			compute_signal(ff,ilast,k);

		}// end of partials or "k" loop

		Normalize_signal();

		fout=&note[ilast-1][0];
		Write_output(fout);

} //end of notes or "i" loop
			
return 0;

/*	this ends the outermost loop of execution. Now we go back to the beginning of the 
	i-loop, change the note and note frequency, and iterate through again.  */
//#endif

}