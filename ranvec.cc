// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/ranvec/ranvec.cc,v 1.6 2008-08-11 12:19:47 weigel Exp $

#include <machine_specific.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#ifndef FAST
#include "ranvec.h"
#endif

//==============================================================================
// initializations

const int r250::nrand = 1000;
const int r250::biginteger = 2147483647;
const double r250::bigfloat = 2147483647.0;
const int r250::nwarm = 10000;
const double r250::multiply = 16807.0;
const double r250::factor = 4.6566128752457969e-10;
const int r250::bigmagic1 = 250;
const int r250::smallmagic1 = 103;
const int r250::bigmagic2 = 521;
const int r250::smallmagic2 = 168;
const int r250::nbit = 32;

//==============================================================================

r250::~r250()
{
  delete[] rand_w_array1;
  delete[] rand_w_array2;
  delete[] random_numbers;
}

//==============================================================================

void r250::setup(int iseed)
{
  double rmod;
  int i, ihlp, imask1, imask2;
  int icyc, ncyc, nrest, ibas1, ibas2, ibas3;
 
  count=nrand;

  if(iseed <=0 || iseed >= biginteger)
    {
      cerr<<"Message from random number initialization:\n";
      cerr<<"Please specify a seed smaller than "<<biginteger<<"\n";
      exit(0);
    }
  
  rmod = static_cast<double>(iseed);

  /* Warm up the congruential generator */

  for(i = 0; i < nwarm; ++i)
    {
      rmod = multiply * rmod;
      rmod = rmod - (static_cast<double>(static_cast<int>(rmod * factor) ) ) * bigfloat;
      ihlp = static_cast<int>(rmod + 0.1);      /* This is done to get rid of */
      rmod = static_cast<double>(ihlp);         /* possible roundoff errors   */
    }

  rand_w_array1 = new int[bigmagic1+nrand];
  for(i=0; i < bigmagic1+nrand; rand_w_array1[i++]=0);
  rand_w_array2 = new int[bigmagic2+nrand];
  for(i=0; i < bigmagic2+nrand; rand_w_array2[i++]=0);
  
  /* Put congruential random numbers onto the working arrays */

  for(i = 0; i < bigmagic1; ++i)
    {
      rmod = multiply * rmod;
      rmod = rmod - ( static_cast<double>( static_cast<int>(rmod * factor) ) ) * bigfloat;
      ihlp = static_cast<int>(rmod + 0.1);
      rmod = static_cast<double>(ihlp);
      rand_w_array1[i] = ihlp;
    }

  for(i = 0; i < bigmagic2; ++i)
    {
      rmod = multiply * rmod;
      rmod = rmod - ( static_cast<double>( static_cast<int>(rmod * factor) ) ) * bigfloat;
      ihlp = static_cast<int>(rmod + 0.1);
      rmod = static_cast<double>(ihlp);
      rand_w_array2[i] = ihlp;
    }

  /* Linear independence of the bit columns for both generators. */
  /* Put ones on the main diagonal, and zeroes above.            */
  /* & is the bitwise AND                                        */
  /* | is the bitwise OR                                         */
  /* ^ is the bitwise XOR                                        */
  
  imask1 = 1;
  imask2 = biginteger;
  for(i = nbit - 2; i > 0; --i)
    {
      rand_w_array1[i] = ( rand_w_array1[i] | imask1 ) & imask2;
      rand_w_array2[i] = ( rand_w_array2[i] | imask1 ) & imask2;
      imask2 = imask2 ^ imask1;
      imask1 = imask1 * 2;
    }
  rand_w_array1[0] = imask1;    /* This last element is treated separately */
  rand_w_array2[0] = imask1;    /* in order to avoid overflow in imask1    */

  /* Warm up. Same structure as in vector_random_generator.      */
  /* Double loop structure to enable vectorization of inner loop */

  /* First, generator one */

  ncyc  = nrand / smallmagic1;
  nrest = nrand - smallmagic1 * ncyc;

  ibas3 = bigmagic1;                 /* position of first new random number */
  ibas2 = bigmagic1 - smallmagic1;   /* position of first input for this    */
  ibas1 = 0;                         /* position of second input for this   */
  
  for(icyc = 0; icyc < ncyc; ++icyc) {
      for(i = 0; i < smallmagic1; ++i)
	rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] ^ rand_w_array1[ibas2 + i];
      ibas1 = ibas1 + smallmagic1;
      ibas2 = ibas2 + smallmagic1;
      ibas3 = ibas3 + smallmagic1;
  }
  
  if(nrest > 0) {
    for(i = 0; i < nrest; ++i)
      rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] ^ rand_w_array1[ibas2 + i];
  }

  /* Put last elements to the beginning */

  for(i = 0; i < bigmagic1; ++i)
    rand_w_array1[i] = rand_w_array1[nrand + i];
  
  /* Now the same for the second generator */

  ncyc  = nrand / smallmagic2;
  nrest = nrand - smallmagic2 * ncyc;
  
  ibas3 = bigmagic2;                 /* position of first new random number */
  ibas2 = bigmagic2 - smallmagic2;   /* position of first input for this    */
  ibas1 = 0;                         /* position of second input for this   */
  
  for(icyc = 0; icyc < ncyc; ++icyc) {
      for(i = 0; i < smallmagic2; ++i)
	  rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] ^ rand_w_array2[ibas2 + i];
      ibas1 = ibas1 + smallmagic2;
      ibas2 = ibas2 + smallmagic2;
      ibas3 = ibas3 + smallmagic2;
  }
  
  if(nrest > 0) {
    for(i = 0; i < nrest; ++i)
      rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] ^ rand_w_array2[ibas2 + i];
  }
  
  /* Put last elements to the beginning */

  for(i = 0; i < bigmagic2; ++i)
    rand_w_array2[i] = rand_w_array2[nrand + i];

  random_numbers = new double[nrand];
  for(i=0; i < nrand; random_numbers[i++]=0);
}

//==============================================================================

void r250::write(const char* file)
{
  FILE *fp;

  fp = fopen(file, "w");
  fwrite(rand_w_array1,sizeof(int),bigmagic1,fp);
  fwrite(rand_w_array2,sizeof(int),bigmagic2,fp);
  fclose(fp);
}

//==============================================================================

void r250::read(const char *file)
{
  FILE *fp;

  fp = fopen(file ,"r");
  fread(rand_w_array1,sizeof(int),bigmagic1,fp);
  fread(rand_w_array2,sizeof(int),bigmagic2,fp);
  fclose(fp);
}

//==============================================================================

void r250::fill_arrays()
{
  int i, icyc, ncyc, nrest, ibas1, ibas2, ibas3;

  /* First, run generator one */

  ncyc  = nrand / smallmagic1;
  nrest = nrand - smallmagic1 * ncyc;

  ibas3 = bigmagic1;                 /* position of first new random number */
  ibas2 = bigmagic1 - smallmagic1;   /* position of first input for this    */
  ibas1 = 0;                         /* position of second input for this   */
  
  for(icyc = 0; icyc < ncyc; ++icyc) {
      for(i = 0; i < smallmagic1; ++i)
	rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] ^ rand_w_array1[ibas2 + i];
      ibas1 = ibas1 + smallmagic1;
      ibas2 = ibas2 + smallmagic1;
      ibas3 = ibas3 + smallmagic1;
  }
  
  if(nrest > 0) {
      for(i = 0; i < nrest; ++i)
	  rand_w_array1[ibas3 + i] = rand_w_array1[ibas1 + i] ^ rand_w_array1[ibas2 + i];
  }

  /* Put last elements to the beginning */
  
  for(i = 0; i < bigmagic1; ++i)
    rand_w_array1[i] = rand_w_array1[nrand + i];

  /* Now the same for the second generator */

  ncyc  = nrand / smallmagic2;
  nrest = nrand - smallmagic2 * ncyc;
  
  ibas3 = bigmagic2;                 /* position of first new random number */
  ibas2 = bigmagic2 - smallmagic2;   /* position of first input for this    */
  ibas1 = 0;                         /* position of second input for this   */
  
  for(icyc = 0; icyc < ncyc; ++icyc) {
    for(i = 0; i < smallmagic2; ++i)
      rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] ^ rand_w_array2[ibas2 + i];
    ibas1 = ibas1 + smallmagic2;
    ibas2 = ibas2 + smallmagic2;
    ibas3 = ibas3 + smallmagic2;
  }
  
  if(nrest > 0) {
    for(i = 0; i < nrest; ++i)
      rand_w_array2[ibas3 + i] = rand_w_array2[ibas1 + i] ^ rand_w_array2[ibas2 + i];
  }
  
  /* Put last elements to the beginning */
  
  for(i = 0; i < bigmagic2; ++i)
    rand_w_array2[i] = rand_w_array2[nrand + i];

  /* Generate normalized random numbers:                        */
  /* Take output from generator one and combine it with         */
  /* that from generator two, via a simple XOR                  */
  
  for(i = 0; i < nrand; ++i)
    random_numbers[i] = factor * (rand_w_array1[i + bigmagic1] ^ rand_w_array2[i + bigmagic2]);
}

//==============================================================================
// initializations
 
const unsigned long mersenne::MATRIX_A = 0x9908b0dfUL;
const unsigned long mersenne::UPPER_MASK = 0x80000000UL;
const unsigned long mersenne::LOWER_MASK = 0x7fffffffUL;

void mersenne::init_genrand(unsigned long s)
{
  mt[0]= s & 0xffffffffUL;
  for (mti=1; mti<N; mti++) {
    mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
    mt[mti] &= 0xffffffffUL;
  }
}

void mersenne::setup(unsigned long init_key[], int key_length)
{
  int i, j, k;
  init_genrand(19650218UL);
  i=1; j=0;
  k = (N>key_length ? N : key_length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))	+ init_key[j] + j;
    mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    i++; j++;
    if (i>=N) { mt[0] = mt[N-1]; i=1; }
    if (j>=key_length) j=0;
  }
  for (k=N-1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
    mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    i++;
    if (i>=N) { mt[0] = mt[N-1]; i=1; }
  }
  
  mt[0] = 0x80000000UL; 
}

unsigned long mersenne::genrand_int32()
{
  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MATRIX_A};

  if (mti >= N) { /* generate N words at one time */
    int kk;
    
    if (mti == N+1)  { /* if init_genrand() has not been called, */
      cout<<"mersenne: using default seed"<<endl;
      init_genrand(5489UL); /* a default initial seed is used */
    }
      
    for (kk=0;kk<N-M;kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<N-1;kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
    
    mti = 0;
  }
  
  y = mt[mti++];
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return y;
}
//==============================================================================

void mersenne::write(const char* file)
{
  FILE *fp;

  fp = fopen(file, "w");
  fwrite(mt,sizeof(unsigned long),N,fp);
  fclose(fp);
}

//==============================================================================

void mersenne::read(const char *file)
{
  FILE *fp;

  fp = fopen(file ,"r");
  fread(mt,sizeof(unsigned long),N,fp);
  mti = N;
  fclose(fp);
}

//==============================================================================
// methods for Poisson distribution

int ranvec::poisson_large(double lambda)
{
    // "Rejection method PA" from "The Computer Generation of Poisson Random Variables" by A. C. Atkinson
    // Journal of the Royal Statistical Society Series C (Applied Statistics) Vol. 28, No. 1. (1979)
    // The article is on pages 29-35. The algorithm given here is on page 32.

    double c = 0.767 - 3.36/lambda;
    double beta = M_PI/sqrt(3.0*lambda);
    double alpha = beta*lambda;
    double k = log(c) - lambda - log(beta);

    for(;;)
    {
        double u = floating();
        double x = (alpha - log((1.0 - u)/u))/beta;
        int n = (int) floor(x + 0.5);
        if (n < 0) continue;
        double v = floating();
        double y = alpha - beta*x;
        double temp = 1.0 + exp(y);
        double lhs = y + log(v/(temp*temp));
        double rhs = k + n*log(lambda) - logfactorial(n);
        if (lhs <= rhs) return n;
    }
}

double ranvec::logfactorial(int n)
{
    if (n > 254)
    {
        double x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*M_PI) + 1.0/(12.0*x);
    }
    else
    {
        double lf[] =        
        {
            0.000000000000000,
            0.000000000000000,
            0.693147180559945,
            1.791759469228055,
            3.178053830347946,
            4.787491742782046,
            6.579251212010101,
            8.525161361065415,
            10.604602902745251,
            12.801827480081469,
            15.104412573075516,
            17.502307845873887,
            19.987214495661885,
            22.552163853123421,
            25.191221182738683,
            27.899271383840894,
            30.671860106080675,
            33.505073450136891,
            36.395445208033053,
            39.339884187199495,
            42.335616460753485,
            45.380138898476908,
            48.471181351835227,
            51.606675567764377,
            54.784729398112319,
            58.003605222980518,
            61.261701761002001,
            64.557538627006323,
            67.889743137181526,
            71.257038967168000,
            74.658236348830158,
            78.092223553315307,
            81.557959456115029,
            85.054467017581516,
            88.580827542197682,
            92.136175603687079,
            95.719694542143202,
            99.330612454787428,
            102.968198614513810,
            106.631760260643450,
            110.320639714757390,
            114.034211781461690,
            117.771881399745060,
            121.533081515438640,
            125.317271149356880,
            129.123933639127240,
            132.952575035616290,
            136.802722637326350,
            140.673923648234250,
            144.565743946344900,
            148.477766951773020,
            152.409592584497350,
            156.360836303078800,
            160.331128216630930,
            164.320112263195170,
            168.327445448427650,
            172.352797139162820,
            176.395848406997370,
            180.456291417543780,
            184.533828861449510,
            188.628173423671600,
            192.739047287844900,
            196.866181672889980,
            201.009316399281570,
            205.168199482641200,
            209.342586752536820,
            213.532241494563270,
            217.736934113954250,
            221.956441819130360,
            226.190548323727570,
            230.439043565776930,
            234.701723442818260,
            238.978389561834350,
            243.268849002982730,
            247.572914096186910,
            251.890402209723190,
            256.221135550009480,
            260.564940971863220,
            264.921649798552780,
            269.291097651019810,
            273.673124285693690,
            278.067573440366120,
            282.474292687630400,
            286.893133295426990,
            291.323950094270290,
            295.766601350760600,
            300.220948647014100,
            304.686856765668720,
            309.164193580146900,
            313.652829949878990,
            318.152639620209300,
            322.663499126726210,
            327.185287703775200,
            331.717887196928470,
            336.261181979198450,
            340.815058870798960,
            345.379407062266860,
            349.954118040770250,
            354.539085519440790,
            359.134205369575340,
            363.739375555563470,
            368.354496072404690,
            372.979468885689020,
            377.614197873918670,
            382.258588773060010,
            386.912549123217560,
            391.575988217329610,
            396.248817051791490,
            400.930948278915760,
            405.622296161144900,
            410.322776526937280,
            415.032306728249580,
            419.750805599544780,
            424.478193418257090,
            429.214391866651570,
            433.959323995014870,
            438.712914186121170,
            443.475088120918940,
            448.245772745384610,
            453.024896238496130,
            457.812387981278110,
            462.608178526874890,
            467.412199571608080,
            472.224383926980520,
            477.044665492585580,
            481.872979229887900,
            486.709261136839360,
            491.553448223298010,
            496.405478487217580,
            501.265290891579240,
            506.132825342034830,
            511.008022665236070,
            515.890824587822520,
            520.781173716044240,
            525.679013515995050,
            530.584288294433580,
            535.496943180169520,
            540.416924105997740,
            545.344177791154950,
            550.278651724285620,
            555.220294146894960,
            560.169054037273100,
            565.124881094874350,
            570.087725725134190,
            575.057539024710200,
            580.034272767130800,
            585.017879388839220,
            590.008311975617860,
            595.005524249382010,
            600.009470555327430,
            605.020105849423770,
            610.037385686238740,
            615.061266207084940,
            620.091704128477430,
            625.128656730891070,
            630.172081847810200,
            635.221937855059760,
            640.278183660408100,
            645.340778693435030,
            650.409682895655240,
            655.484856710889060,
            660.566261075873510,
            665.653857411105950,
            670.747607611912710,
            675.847474039736880,
            680.953419513637530,
            686.065407301994010,
            691.183401114410800,
            696.307365093814040,
            701.437263808737160,
            706.573062245787470,
            711.714725802289990,
            716.862220279103440,
            722.015511873601330,
            727.174567172815840,
            732.339353146739310,
            737.509837141777440,
            742.685986874351220,
            747.867770424643370,
            753.055156230484160,
            758.248113081374300,
            763.446610112640200,
            768.650616799717000,
            773.860102952558460,
            779.075038710167410,
            784.295394535245690,
            789.521141208958970,
            794.752249825813460,
            799.988691788643450,
            805.230438803703120,
            810.477462875863580,
            815.729736303910160,
            820.987231675937890,
            826.249921864842800,
            831.517780023906310,
            836.790779582469900,
            842.068894241700490,
            847.352097970438420,
            852.640365001133090,
            857.933669825857460,
            863.231987192405430,
            868.535292100464630,
            873.843559797865740,
            879.156765776907600,
            884.474885770751830,
            889.797895749890240,
            895.125771918679900,
            900.458490711945270,
            905.796028791646340,
            911.138363043611210,
            916.485470574328820,
            921.837328707804890,
            927.193914982476710,
            932.555207148186240,
            937.921183163208070,
            943.291821191335660,
            948.667099599019820,
            954.046996952560450,
            959.431492015349480,
            964.820563745165940,
            970.214191291518320,
            975.612353993036210,
            981.015031374908400,
            986.422203146368590,
            991.833849198223450,
            997.249949600427840,
            1002.670484599700300,
            1008.095434617181700,
            1013.524780246136200,
            1018.958502249690200,
            1024.396581558613400,
            1029.838999269135500,
            1035.285736640801600,
            1040.736775094367400,
            1046.192096209724900,
            1051.651681723869200,
            1057.115513528895000,
            1062.583573670030100,
            1068.055844343701400,
            1073.532307895632800,
            1079.012946818975000,
            1084.497743752465600,
            1089.986681478622400,
            1095.479742921962700,
            1100.976911147256000,
            1106.478169357800900,
            1111.983500893733000,
            1117.492889230361000,
            1123.006317976526100,
            1128.523770872990800,
            1134.045231790853000,
            1139.570684729984800,
            1145.100113817496100,
            1150.633503306223700,
            1156.170837573242400,
        };
        return lf[n];
    }
}
