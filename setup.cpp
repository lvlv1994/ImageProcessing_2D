#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
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
void myround(double **f,int width,int height);
void round3a(double **f, int width, int height);
void round3ab(double **f, int width, int height);
void round4a(double **f, int width, int height);
double g(double c, double **f, int width, int height);
double CH(double **f, int width, int height, double c);
double hist(double **f, int width, int height, double c);
void roundhw1p2(double **f, int width, int height);
double cubicInterpolate(double y1,double y2,double y3,double y4,double x);
void problem1(double **f, int width, int height, int N, double percentage);
//MAIN BODY HERE
int main()
{
   int width;
   int height;
   double **f , **newf, **mask;
   double h;
   int i, j;

// if reading in picture
   ifstream infile("seal.bmp", ios::in | ios::binary);
   f = readpicture(infile,width,height);
   infile.close();
   //myround(f,width,height);
   //roundhw1p2(f,width,height);
   //round3a(f,width,height);
   //round3ab(f,width,height);
   //g(1.0,f,width,height);
  /* mask = matrix(width,height);
   newf = matrix(width,height);
   for (i = 0;i<width;i++)
 for(j = 0;j<height;j++)
 {
    if ((double(i)/4 - i/4 == 0)||
        (double(j)/9 - j/9 == 0)||
        (double(j+1)/9 - (j+1)/9 == 0)||
        (double(j+2)/9 - (j+2)/9 == 0)||
        (i<=50&&i>=10)||(j>=10&&j<=30))
     mask[i][j] = 1;
    else
     mask[i][j] = 0;
 }

 for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
        if (mask[i][j] == 1)
          f[i][j] = 0.5;

*/
















// if creating picture
 /*  width = 500;
   height = 500;
   f = matrix(width,height);
   mask = matrix(width,height);
   newf = matrix(width,height);
   double m;

   m = 1.0/500;
   for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
   {
         if ((0.5+i)*m*2-1 <= -0.5)
              f[i][j] = 0.25;
         if ((0.5+i)*m*2-1 > -0.5 && (0.5+i)*m*2-1 < 0.5)
              f[i][j] = 0.75;
         if ((0.5+i)*m*2-1 > 0.5)
              f[i][j] = 0.4;

 }

   for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
        {
          if (((0.5+i)*m*2-1 >= -0.2 &&
                (0.5+i)*m*2-1 <=0.2)&&
              ((0.5+j)*m*2-1 >= -0.01 &&
               (0.5+j)*m*2-1 <= 0.01))

             mask[i][j] = 1;

         else
             mask[i][j] = 0;
 }

 for(i = 0; i < width; i++)
    for(j = 0; j < height; j++)
        if (mask[i][j] == 1)
          f[i][j] = 0.5;
*/

    //at k = h^2/4 and using Von Newmman condition,
    //this is a routine heat flow

 /*  int n, maxstep = 10000;
    newf = matrix(width,height);

    for(n = 1; n <= maxstep; n++)
    {
        cout<<"step = "<<n<<endl;

        for(i = 0; i < width; i++)
         for(j = 0; j < height; j++)
            if (mask[i][j] == 1)
                newf[i][j] = (f[max(i-1,0)][j]+f[min(i+1,0)][j]+f[i][max(j-1,0)]+
                              f[i][min(j+1,height-1)])/4.0;
        for(i = 0; i < width; i++)
         for(j = 0; j < height; j++)
            if (mask[i][j] == 1)
            f[i][j] = newf[i][j];
    }
*/
    // this is the routine TV minimization flow
    int n, r, s, maxstep = 0;
    double RHS, **df2, df[2], gradf2, mess, k;
    df2 = matrix(2,2);
    h = 0.1;
    k = h*h*h/4.0;
    for( n = 1; n <= maxstep; n++)
    {
        cout<<"step = "<< n<<endl;
        for(i = 0; i < width; i++)
         for(j = 0; j < height; j++)
            if (mask[i][j] == 1)
            {
                df[0] = (f[min(i+1,width-1)][j]-f[max(i-1,0)][j])/(2.0*h);
                df[1] = (f[i][min(j+1,height-1)]-f[i][max(j-1,0)])/(2.0*h);
                df2[0][0] = (f[min(i+1,width-1)][j]-2.0*f[i][j]+f[max(i-1,0)][j])/(h*h);
                df2[1][1] = (f[i][min(j+1,height-1)]-2.0*f[i][j]+f[i][max(j-1,0)])/(h*h);
                df2[0][1] = (f[min(i+1,width-1)][min(j+1,height-1)]-
                             f[min(i+1,width-1)][max(j-1,0)]-
                             f[max(i-1,0)][min(j+1,height-1)]+
                             f[max(i-1,0)][max(j-1,0)])/(4.0*h*h);
                df2[1][0] = df2[0][1];
                gradf2 = df[0]*df[0]+df[1]*df[1];
                mess = 0.0;
                for(r = 0; r < 2; r++)
                  for(s = 0; s < 2; s++)
                    mess += df[r]*df2[r][s]*df[s];

                RHS = (df2[0][0]+df2[1][1]-mess/(gradf2+h*h))/sqrt(gradf2+h*h);
                newf[i][j] = f[i][j]+k*RHS;

                }

        for(i = 0; i < width; i++)
         for(j = 0; j < height; j++)
            if (mask[i][j] == 1 )
            f[i][j] = newf[i][j];


    }


   //int newheight , newwidth;






// writing picture
   ofstream outfile("result.bmp", ios::out | ios::binary);
   writepicture(outfile,f,width,height);
   outfile.close();

   free_matrix(f,width,height);
   /* Windows only */


system("pause");
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

void myround(double **f,int width,int height)
{

    int i,j, countblack = 0, countwhite = 0;
    for (i = 0; i < width; i++ ){
        for(j = 0; j < height; j++){
        if(f[i][j]<=0.5){
                f[i][j]=0;
                countblack++;
        }
        else{
                f[i][j]=1;
                countwhite++;
        }
}
    }
    cout << "number of black pixels    " <<countblack;
    cout <<"   number of white pixels    "<<countwhite;
}

void roundhw1p2(double **f, int width, int height){
    int i,j,ch[8];
    double h = 0.05;

 for (i = 0; i < width; i++)
 for (j = 0; j < height; j++)
 {
 if (f[i][j] <= 3.5*h)
 {
 f[i][j] = 0;
 ch[0]++;
 }
 else if (f[i][j] <= 7.5*h && f[i][j] > 3.5*h)
 {
 f[i][j] = 7*h;
 ch[1]++;
 }
 else if (f[i][j] <= 8.5*h && f[i][j] > 7.5*h)
 {
 f[i][j] = 8*h;
 ch[2]++;
 }
 else if (f[i][j] <= 10*h && f[i][j] > 8.5*h)
 {
 f[i][j] = 9*h;
 ch[3]++;
 }
 else if (f[i][j] <= 11.5*h && f[i][j] > 10*h)
 {
 f[i][j] = 11*h;
 ch[4]++;
 }
 else if (f[i][j] <= 12.5*h && f[i][j] > 11.5*h)
 {
 f[i][j] = 12*h;
 ch[5]++;
 }
 else if (f[i][j] <= 16.5*h && f[i][j] > 12.5*h)
 {
 f[i][j] = 13*h;
 ch[6]++;
 }
 else
 {
 f[i][j] = 1;
 ch[7]++;
 }
 }

 cout << "The number of pixels with intensity of 0 is " << ch[0] << endl;

 cout << "The number of pixels with intensity of 7h is " << ch[1] << endl;

 cout << "The number of pixels with intensity of 8h is " << ch[2] << endl;

 cout << "The number of pixels with intensity of 9h is " << ch[3] << endl;

 cout << "The number of pixels with intensity of 11h is " << ch[4] << endl;

 cout << "The number of pixels with intensity of 12h is " << ch[5] << endl;

 cout << "The number of pixels with intensity of 13h is " << ch[6] << endl;

 cout << "The number of pixels with intensity of 1 is " << ch[7] << endl;}

double cubicInterpolate(double y1,double y2,double y3,double y4,double x)
{
 return  -y1*(x-1.5)*(x-2.5)*(x-3.5)/6 + y2*(x-0.5)*(x-2.5)*(x-3.5)/2 - y3*(x-0.5)*(x-1.5)*(x-3.5)/2 +
y4*(x-0.5)*(x-1.5)*(x-2.5)/6;
}


void round3a(double **f, int width, int height)
{
    int i,j;
    for(i = 0; i < width; i++){
        for(j = 0; j < height; j++){
            f[i][j]=0.45*f[i][j]+0.55;
        }
    }
}

void round3ab(double **f, int width, int height)
{
    int i,j;
    for(i = 0; i < width; i++){
        for(j = 0; j < height; j++){
            f[i][j]=0.45*f[i][j];
        }
    }
}

void round4a(double **f, int width, int height)
{
    double min,max;
    int i,j;
    min = 1.0;
    max = 0.0;
    //loop to find extrema
    for(i = 0; i < width; i++){
        for(j = 0; j < height; j++){
            if (f[i][j] < min){min = f[i][j];}
            if (f[i][j] > max){max = f[i][j];}
        }
    }
        for(i = 0; i < width; i++){
        for(j = 0; j < height; j++){
            f[i][j] = (f[i][j]-min)/(max-min);
        }
    }
}

double g(double c, double **f, int width, int height)
{
    double min,max;
    int i,j;
    min = 1.0;
    max = 0.0;
    //loop to find extrema
    for(i = 0; i < width; i++){
        for(j = 0; j < height; j++){
            if (f[i][j] < min){min = f[i][j];}
            if (f[i][j] > max){max = f[i][j];}
        }
    }
    //cout <<"the max is " << max << endl;
    //cout <<"the min is " << min << endl;
    return (CH(f,width,height,c)-CH(f,width,height,min))/
            (CH(f,width,height,max)-CH(f,width,height,min));
}

double CH(double **f, int width, int height, double c)
{
    int i,j,count = 0;
    for(i = 0; i < width; i++){
        for(j = 0; j < height; j++){
            if (f[i][j] <= c){count++;}
        }
    }
    return static_cast<double>(count);
}
double hist(double **f, int width, int height, double c)
{
    int i,j,count=0;
    double tol = 1.0e-10;
    for(i=0;i<width;i++){
        for(j=0;j<height;j++){
            if(f[i][j]>c-tol && f[i][j]<c+tol){count++;}
    }

    return static_cast<double> (count);
}
}
void problem1(double **f, int width, int height, int N, double percentage)
{
    int *p = (int *)malloc(sizeof(int)*N);
    int i,j;
    for(i = 0; i < N; i++)
        p[i] = 0;

    double h = 1.0/(N-1);
    for(j = 0; j < width; j++)
        for(i = 0; i < height; i++)
    {
        p[(int)(f[j][i]/h)]++;

    }
    int sum = 0;
    for(i = 0; i < N; i++)
    {
        sum += p[i];
        p[i] = sum;
    }
    double *dp = (double *)malloc(sizeof(double)*N);
    for(i=0;i<N;i++)
        dp[i] = (double)(p[i])/(width*height);

    double a = 0.0, b = 1.0, t;
    for(i=0;i<N;i++)
    {
        if(dp[i]>=percentage)
        break;
    }

    if(i>0)
        t = (percentage - dp[i-1])/(dp[i]-dp[i-1]);

    a = (i-1+t)*h;
    printf("a=%f/n",a);

    for(i = 0; i < N; i++)
    {
        if(dp[i]>=1-percentage)
        break;
    }
    if(i>0)
        t = (1-percentage-dp[i-1])/(dp[i]-dp[i-1]);

    b = (i-1+t)*h;
    printf("b=%f/n",b);

    for(j = 0; j < width; j++)
    for(i = 0; i < height; i++)
    {
        if(f[j][i]<a) f[j][i] = 0.0;
    else if(f[j][i]>0) f[j][i] = 1.0;
    else f[j][i] = (f[j][i]-a)/(b-a);
    }
    free(p);
    free(dp);

}

