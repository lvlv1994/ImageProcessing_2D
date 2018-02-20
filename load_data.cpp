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
//   YOUR FUNCTIONS HERE

//MAIN BODY HERE
int main()
{
   int width;
   int height;
   double **f;
   double x, y;
   int i, j;

// if reading in picture
   ifstream infile("deadflowers.bmp", ios::in | ios::binary);
   f = readpicture(infile,width,height);
   infile.close();
// if creating picture
//   width = 100;
//   height = 100;
//   f = matrix(width,height);

// writing picture
/*   int a= 0,b =0;
   for (i = 0; i< width; i++) for( j = 0; j< height; j++) {
       if (f[i][j] < 0.5) {
           f[i][j] = 0;
           a++;
       }
       else {
           f[i][j] = 1;
           b++;
       }
   }
   cout << "The number of pixels with intensity of 1 is " << a << endl;
   cout << "The number of pixels with intensity of 0 is " << b << endl;
*/
/*    int counter[8] = { };
    float h = 20;
    double e = 7.5 * h;


    for (i = 0; i < width; i++)
        for (j = 0; j < height; j++) {
            if (f[i][j] <= 3.5 / h) {
                f[i][j] = 0;
                counter[0]++;
            }
            else if (f[i][j] <= (float)7.5 / h && f[i][j] > (float)3.5/h) {
                f[i][j] = 7*h;
                counter[1]++;
            }
            else if (f[i][j] <= (float)8.5 / h && f[i][j] > (float)7.5/h) {
                f[i][j] = 8*h;
                counter[2]++;
            }
            else if (f[i][j] <= (float)10 / h && f[i][j] > (float)8.5/h) {
                f[i][j] = 9*h;
                counter[3]++;
            }
            else if (f[i][j] <= (float)11.5 / h && f[i][j] > (float)10/h) {
                f[i][j] = 11*h;
                counter[4]++;
            }
            else if (f[i][j] <= (float)12.5 / h && f[i][j] > (float)11.5/h) {
                f[i][j] = 12*h;
                counter[5]++;
            }
            else if (f[i][j] <= (float)16.5 / h && f[i][j] > (float)12.5/h) {
                f[i][j] = 13*h;
                counter[6]++;
            }
            else {
                f[i][j] = 1;
                counter[7]++;
            }
    }
    cout << "The number of pixels with intensity of 0 is " << counter[0] << endl;
    cout << "The number of pixels with intensity of 7h is " << counter[1] << endl;
    cout << "The number of pixels with intensity of 8h is " << counter[2] << endl;
    cout << "The number of pixels with intensity of 9h is " << counter[3] << endl;
    cout << "The number of pixels with intensity of 11h is " << counter[4] << endl;
    cout << "The number of pixels with intensity of 12h is " << counter[5] << endl;
    cout << "The number of pixels with intensity of 13h is " << counter[6] << endl;
    cout << "The number of pixels with intensity of 1 is " << counter[7] << endl;

    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            for (int n = 0; n < 256; n++) {
                if(f[i][j]==(n)/255.0)
                f[i][j] = n/255.0;
            }
            f[i][j] = 0.45*f[i][j]+0.55 ;
        }
    }
 
    double max = 0;
    double min = 0;
    for (i = 0; i < width; i++)
        for (j = 0; j < height; j++)
        {
            if (f[i][j] < min) {
                min = f[i][j];
            }
        
            if (f[i][j] > max) {
                 max = f[i][j];
                
            }
       
    }
    cout << "b is "<<max << endl;
    cout << "a is "<<min << endl;
   */
    for (i = 0; i < width; i++)
        for (j = 0; j < height; j++)
            {
                f[i][j]*= 1/0.458824;
            }
    // writing picture
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
//YOUR FUNCTIONS HERE
//**************************************************************************
