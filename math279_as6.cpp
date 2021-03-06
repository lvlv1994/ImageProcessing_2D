#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

//FUNCTION HEADERS HERE
//   FUNCTIONS FOR READING AND WRITING
double **readpicture(ifstream &infile, int &width, int &height);
char readchar(ifstream& infile, int &num);
int readint(ifstream& infile, int numbytes, int &num);
int char2int(char thechar);
int char2int(char *thechar, int numbytes);
void readrgb(ifstream& infile, double *rgb, int &num);
void writepicture(ofstream& outfile, double **f, int width, int height);
char int2char(unsigned long theint);
void int2char(char *thechar, int theint);
void writechar(ofstream &outfile, char thechar);
void writeint(ofstream &outfile, int theint, int numbytes);
void writergb(ofstream& outfile, double r, double g, double b);
//   FUNCTIONS FOR CREATING AND REMOVING MEMORY
double **matrix(int width, int height);
void free_matrix(double **f, int width, int height);
//   RANDOM FUNCTIONS
double ran2(int &idum);
double gasdev(int &idum);
//   YOUR FUNCTIONS HERE

//MAIN BODY HERE
int main()
{
   int width;
   int height;
   double **f;
   double x, y;
   int i, j;
   int idum = -1;

// if reading in picture
   ifstream infile("bread.bmp", ios::in | ios::binary);
   f = readpicture(infile,width,height);
   infile.close();
    for(i = 0; i < width;i++)
        for (j = 0; j < height ; j++) {
            f[i][j]+=0.1*gasdev(idum);
    }
    
   /* for(i = 0; i < width;i++)
        for (j = 0; j < height ; j++) {
            f[i][j]+=10*gasdev(idum);
        }
   */
    /*for(i = 0; i < width;i++)
        for (j = 0; j < height ; j++) {
            if (fabs(gasdev(idum)) < 0.9)
                f[i][j] = 1;
        }
  */
   /*for(i = 0; i < width;i++)
        for (j = 0; j < height ; j++) {
            if (fabs(gasdev(idum)) < 0.9) f[i][j] = 0;
        }*/
    // if creating picture
//   width = 100;
//   height = 100;
//   f = matrix(width,height);

// writing picture
/*
    double h = 0.05;
    double eps = h;
    double s = pow(h,3)/4;
    double fx,fy,fxx,fyy,fxy;
    for(int k = 0; k < 200; k++) {
        //cout<<"k= "<<k<<endl;
      for (i = 0;i<width;i++)
        for(j = 0;j<height;j++) {
            fx = (f[min(i+1,width-1)][j]-f[max(i-1,0)][j])/(2*eps);
            fy = (f[i][min(j+1,height-1)]-f[i][max(j-1,0)])/(2*eps);
            fxx = (f[min(i+1,width-1)][j]-2*f[i][j]+f[max(i-1,0)][j])/(eps*eps);
            fyy = (f[i][min(j+1,height-1)]-2*f[i][j]+f[i][max(j-1,0)])/(eps*eps);
            fxy = (f[min(i+1,width-1)][min(j+1,height-1)]-f[min(i+1,width-1)][max(j-1,0)]-f[max(i-1,0)]
                   [min(j+1,height-1)]+f[max(i-1,0)][max(j-1,0)])/(4*eps*eps);
            f[i][j] = f[i][j]+s*((fxx+fyy)/sqrt(fx*fx+fy*fy+eps*eps)-(fx*(fxx*fx+fxy*fy)+fy*(fxy*fx+fyy*fy))/
                                 sqrt(pow(fx*fx+fy*fy+eps*eps,3))); }
    }
 */
    double h = 0.05;
    double eps;
    double s;
    double fx,fy,fxx,fyy,fxy;
    eps = pow(10,-9);
    s = h*h/4;
    for(int k = 0; k < 400; k++) {
        //cout<<"k= "<<k<<endl;
        for (i = 0;i<width;i++)
        for(j = 0;j<height;j++) {
            fx = (f[min(i+1,width-1)][j]-f[max(i-1,0)][j])/(2*h);
            fy = (f[i][min(j+1,height-1)]-f[i][max(j-1,0)])/(2*h);
            fxx = (f[min(i+1,width-1)][j]-2*f[i][j]+f[max(i-1,0)][j])/(h*h);
            fyy = (f[i][min(j+1,height-1)]-2*f[i][j]+f[i][max(j-1,0)])/(h*h);
            fxy = (f[min(i+1,width-1)][min(j+1,height-1)]-f[min(i+1,width-1)][max(j-1,0)]-f[max(i-1,0)]
                   [min(j+1,height-1)]+f[max(i-1,0)][max(j-1,0)])/(4*h*h);
            f[i][j] = f[i][j]+s*((fxx+fyy)-(fx*(fxx*fx+fxy*fy)+fy*(fxy*fx+fyy*fy))/(fx*fx+fy*fy+eps*eps));
        } }
    cout<<f[1][1]<<endl;
   ofstream outfile("result.bmp", ios::out | ios::binary);
   writepicture(outfile,f,width,height);
   outfile.close();

   free_matrix(f,width,height);
}

//FUNCTION CONTENT HERE

//**************************************************************************
//FUNCTIONS FOR READING AND WRITING
//**************************************************************************
double **readpicture(ifstream &infile, int &width, int &height)
{
   int i, j, k;
   char junk, theformat[2];
   double rgb[3], **f;
   int num = 0;
   int upsidedown = 0;
   cout << "reading picture" << endl;

   for (i = 0; i < 2; i++)
      theformat[i] = readchar(infile,num);
   int filesize = readint(infile,4,num);
   cout << "   numbytes = " << filesize << endl;
   for (i = 0; i < 4; i++)
      junk = readchar(infile,num);
   int offset = readint(infile,4,num);
   int headerstart = num;
   int headersize = readint(infile,4,num);
   if (headersize == 12)
   {
      width = readint(infile,2,num);
      height = readint(infile,2,num);
   }
   else
   {
      width = readint(infile,4,num);
      height = readint(infile,4,num);
   }
   if (height < 0)
   {
      height = abs(height);
      upsidedown = 1;
   }
   cout << "   width = " << width << endl;
   cout << "   height = " << height << endl;
   int numcolorpanes = readint(infile,2,num);
   int bitsperpixel = readint(infile,2,num);
   if (bitsperpixel != 24)
   {
      cout << "ERROR: this program is only built for 1 byte per rgb, not a total of " 
           << bitsperpixel << " bits" << endl;
      exit(1);
   }
   for (i = num-headerstart; i < headersize; i++)
      junk = readchar(infile,num);

   f = matrix(width,height);

   if (upsidedown)
      for (j = height-1; j >= 0; j--)
      {
         for (i = 0; i < width; i++)
         {
            readrgb(infile,rgb,num);
            f[i][j] = 0.0;
            for (k = 0; k < 3; k++)
               f[i][j] += rgb[k]*rgb[k];
            f[i][j] = sqrt(f[i][j])/sqrt(3.0);
         }
         for (i = 0; i < (4-(3*width)%4)%4; i++)
            junk = readchar(infile,num);
      }
   else
      for (j = 0; j < height; j++)
      {
         for (i = 0; i < width; i++)
         {
            readrgb(infile,rgb,num);
            f[i][j] = 0.0;
            for (k = 0; k < 3; k++)
               f[i][j] += rgb[k]*rgb[k];
            f[i][j] = sqrt(f[i][j])/sqrt(3.0);
         }
         for (i = 0; i < (4-(3*width)%4)%4; i++)
            junk = readchar(infile,num);
      }

   return f;
}

int readint(ifstream& infile, int numbytes, int &num)
{
   char *temp;
   int i;
   int theint;

   temp = new char[numbytes];

   for (i = 0; i < numbytes; i++)
      infile.read(reinterpret_cast<char *>(&(temp[i])),sizeof(char));
   theint = char2int(temp,numbytes);

   delete[] temp;
   num += numbytes;

   return theint;
}

int char2int(char thechar)
{
   int i, theint, imask;
   char cmask;

   cmask = 1;
   imask = 1;
   theint = 0;
   for (i = 0; i < 8; i++)
   {
      if (thechar & cmask)
      {
         theint += imask;
      }
      cmask = cmask << 1;
      imask *= 2;
   }

   return theint;
}

int char2int(char *thechar, int numbytes)
{
   int i;
   int theint, power;

   power = 1;
   theint = 0;
   for (i = 0; i < numbytes; i++)
   {
      theint += char2int(thechar[i])*power;
      power = power*256;
   }

   return theint;
}

char readchar(ifstream& infile, int &num)
{
   char thechar;

   infile.read(reinterpret_cast<char *>(&thechar),sizeof(char));
   num++;

   return thechar;
}

void readrgb(ifstream& infile, double *rgb, int &num)
{
   char ctemp;
   int i, itemp, ijunk;

   for (i = 0; i < 3; i++)
   {
//      ctemp = readchar(infile);
//      itemp = char2int(ctemp);
      itemp = readint(infile,1,ijunk);
      rgb[i] = static_cast<double>(itemp)/255.0;
      num++;
   }
}

void writepicture(ofstream& outfile, double **f, int width, int height)
{
   int numbytes;
   int i, j;
   cout << "writing picture" << endl;

   writechar(outfile,66);
   writechar(outfile,77);

// number of bytes
   numbytes = 54+height*(3*width+(4-(3*width)%4)%4);
   cout << "   numbytes = " << numbytes << endl;
   writeint(outfile,numbytes,4);

   writeint(outfile,0,2);
   writeint(outfile,0,2);
   writeint(outfile,54,4);
   writeint(outfile,40,4);

// width
   cout << "   width = " << width << endl;
   writeint(outfile,width,4);
// height
   cout << "   height = " << height << endl;
   writeint(outfile,height,4);

   writeint(outfile,1,2);
   writeint(outfile,24,2);
   writeint(outfile,0,4);
   writeint(outfile,16,4);
   writeint(outfile,2835,4);
   writeint(outfile,2835,4);
   writeint(outfile,0,4);
   writeint(outfile,0,4);

   for (j = 0; j < height; j++)
   {
      for (i = 0; i < width; i++)
         writergb(outfile,f[i][j],f[i][j],f[i][j]);
      for (i = 0; i < (4-(3*width)%4)%4; i++)
         writechar(outfile,0);
   }
}

char int2char(unsigned long theint)
{
   char thechar, mask;
   int i;

   mask = 1;
   thechar = 0;
   for (i = 0; i < 8; i++)
   {
      if (theint%2)
         thechar += mask;
      theint /= 2;
      mask = mask << 1;
   }

   return thechar;
}

void int2char(char *thechar, int theint, int numbytes)
{
   int temp = theint;
   int i;

   for (i = 0; i < numbytes; i++)
   {
      thechar[i] = int2char(temp%256);
      temp = temp/256;
   }
}

void writechar(ofstream& outfile, char thechar)
{
   outfile.write(reinterpret_cast<char *>(&thechar),sizeof(char));
}

void writeint(ofstream &outfile, int theint, int numbytes)
{
   char *temp = new char[numbytes];
   int i;

   int2char(temp,theint,numbytes);
   for (i = 0; i < numbytes; i++)
      outfile.write(reinterpret_cast<char *>(&(temp[i])),sizeof(char));

   delete[] temp;
}

void writergb(ofstream& outfile, double r, double g, double b)
{
   int i, irgb;
   char temp;
   char a;
   double rgb[3];

   rgb[0] = r;
   rgb[1] = g;
   rgb[2] = b;
   for (i = 0; i < 3; i++)
   {
      irgb = static_cast<int>(floor(255.0*rgb[i]+0.5));
      if (irgb < 0)
         irgb = 0;
      else if (irgb > 255)
         irgb = 255;

      temp = int2char(irgb);
      outfile.write(reinterpret_cast<char *>(&temp),sizeof(char));
   }
}

//**************************************************************************
//FUNCTIONS FOR CREATING AND REMOVING MEMORY
//**************************************************************************

double **matrix(int width, int height)
{
   double **f;
   int i;

   f = new double*[width];
   for (i = 0; i < width; i++)
      f[i] = new double[height];

   return f;
}

void free_matrix(double **f, int width, int height)
{
   int i;

   for (i = 0; i < width; i++)
      delete f[i];
   delete f;
}

//**************************************************************************
//RANDOM FUNCTIONS
//**************************************************************************

double ran2(int &idum)
{
   const int IM1 = 2147483563;
   const int IM2 = 2147483399;
   const double AM = (1.0/IM1);
   const int IMM1 = (IM1-1);
   const int IA1 = 40014;
   const int IA2 = 40692;
   const int IQ1 = 53668;
   const int IQ2 = 52774;
   const int IR1 = 12211;
   const int IR2 = 3791;
   const int NTAB = 32;
   const int NDIV = (1+IMM1/NTAB);
   const double EPS = 1.2e-7;
   const double RNMX = (1.0-EPS);
   int j;
   int k;
   static int idum2 = 123456789;
   static int iy = 0;
   static int iv[NTAB];
   double temp;


   if (idum <= 0)
   {
      if (-idum < 1)
         idum = 1;
      else
         idum = -idum;
      idum2 = idum;
      for (j = NTAB+7; j >= 0; j--)
      {
         k = idum/IQ1;
         idum = IA1*(idum-k*IQ1)-k*IR1;
         if (idum < 0)
            idum += IM1;
         if (j < NTAB)
            iv[j] = idum;
      }
      iy = iv[0];
   }
   k = idum/IQ1;
   idum = IA1*(idum-k*IQ1)-k*IR1;
   if (idum < 0)
      idum += IM1;
   k = idum2/IQ2;
   idum2 = IA2*(idum2-k*IQ2)-k*IR2;
   if (idum2 < 0)
      idum2 += IM2;
   j = iy/NDIV;
   iy = iv[j]-idum2;
   iv[j] = idum;
   if (iy < 1)
      iy += IMM1;
   if ((temp = AM*iy) > RNMX)
      return RNMX;
   else
      return temp;
}

// normally distributed random numbers with mean 0, standard deviation 1.
double gasdev(int &idum)
{
   int iset = 0;
   double gset;
   double fac, rsq, v1, v2;

   if (iset == 0)
   {
      do
      {
         v1 = 2.0*ran2(idum)-1.0;
         v2 = 2.0*ran2(idum)-1.0;
         rsq = v1*v1+v2*v2;
      }
      while (rsq >= 1.0 || rsq == 0.0);
      fac = sqrt(-2.0*log(rsq)/rsq);
      gset = v1*fac;
      iset = 1;
      return v2*fac;
   }
   else
   {
      iset = 0;
      return gset;
   }
}

//**************************************************************************
//YOUR FUNCTIONS HERE
//**************************************************************************
