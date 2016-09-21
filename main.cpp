/* 
 * File:   main.cpp
 * Author: jcoelho
 *
 * Created on September 11, 2016, 2:18 PM
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Image.h"
#include "Cluster.h"

using namespace std;



/**
 * TODO: converte uma imagem rgb em uma imagem Lab.
 * @param rgb - imagem rgb de entrada.
 * @return - imagem em Lab.
 */
Image convertImageFromRGB2Lab( const Image& rgb )
{
    printf("Converting to Lab\n");
    int w = rgb.getW();
    int h = rgb.getH();
    Image labImg = Image(w,h);
    for(int i=0;i<w;i++)
    {
        for(int j=0;j<h;j++)
        {
            Pixel p = rgb.getPixel(i,j);
            //printf("%d%d\n",i,j);
            p=rgb.rgbToXYZ(p);
            p=rgb.XYZToLab(p);
            int ind = rgb.computePosition( i, j );
            //printf("p: %f %f %f ",p[0],p[1],p[2]);
            labImg.setPixel( ind, p );
            //printf("ind %d\n",ind );
        }
        //printf("%d of %d\n",i,w);
    }

    return labImg;

}



/**
 * TODO: converte uma imagem Lab em uma imagem em rgb.
 * @param Lab - imagem Lab de entrada.
 * @return - imagem em RGB.
 */
Image convertImageFromLAB2RGB( const Image& Lab )
{
    printf("Re-Converting to RGB\n");
    int w = Lab.getW();
    int h = Lab.getH();
    Image rgbImg = Image(w,h);
    for(int i=0;i<w;i++)
    {
        for(int j=0;j<h;j++)
        {
            Pixel p = Lab.getPixel(i,j);
            //printf("%d%d\n",i,j);
            p=Lab.LabToXYZ(p);
            p=Lab.XYZTorgb(p);
            //p[2]*=1.1;
            //p[0]*=1.1;
            //p[1]*=1.1;
            //printf("p: %f %f %f ",p[0],p[1],p[2]);
            int ind = Lab.computePosition( i, j );
            rgbImg.setPixel( ind, p );
            //printf("ind %d\n",ind );
        }
        //printf("%d of %d\n",i,w);
    }

    return rgbImg;
}



/**
 * TODO: inicializa os clusters.
 * @param clusters - clusters que DEVEM ser alocados e inicializados.
 * @param Lab - imagem em Lab.
 * @param k - numero desejado de superpixels.
 * @return - numero de superpixels.
 */
int initializeClusters( Cluster*& clusters, Image& Lab, int k )
{
    int nk = k;
    clusters = new Cluster[k];
    int w = Lab.getW();
    int h = Lab.getH();
    int size = sqrt((w * h)/k);

    float x=0.0;
    float y=0.0;
    for(int i =0 ; i< k; i++)
    {
        clusters[i] = Cluster();
        Pixel p = Lab.getPixel((int)x,(int)y);
        clusters[i].set(p,x,y);
        x += size;
        if(x>=w)
        {
            y+=size;
            x=0.0;
        }
    }


    return nk;
}


double distanceFromSuperPixels(Cluster sp,Pixel p, int x, int y,double mc, double ms)
{
    Pixel c = sp.getPixel();
    double dc = sqrt((c[0] - p[0])*(c[0] - p[0])+(c[1] - p[1])*(c[1] - p[1])+(c[2] - p[2])*(c[2] - p[2]));
    double ds = sqrt((sp.getX() - x)*(sp.getX() - x)+(sp.getY() - y)*(sp.getY() - y));

    double dt = sqrt( (dc/mc)*(dc/mc) + (ds/ms)* (ds/ms)) ;

    return dt;
}

/**
 * TODO: realiza o algoritmo de superpixels.
 * @param Lab - Imagem em lab.
 * @param clusters - clusters inicializados.
 * @param labels - labels inicializadas.
 * @param k - numero de superpixels.
 * @param M - compacidade.
 */
void performSuperPixelsAlgorithm( Image& Lab, Cluster* clusters, int *labels, int k, double M )
{



}



void drawContoursAroundSegments( Image& rgb, int* labels, Pixel borderColor = Pixel( ) )
{
    int w = rgb.getW( );
    int h = rgb.getH( );

    const int dx8[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

    int sz = w * h;
    vector<bool> istaken( sz, false );
    vector<int> contourx( sz );
    vector<int> contoury( sz );
    int mainindex( 0 );
    int cind( 0 );
    for (int j = 0; j < h; j++)
    {
        for (int k = 0; k < w; k++)
        {
            int np( 0 );
            for (int i = 0; i < 8; i++)
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if (( x >= 0 && x < w ) && ( y >= 0 && y < h ))
                {
                    int index = y * w + x;

                    if (false == istaken[index])//comment this to obtain internal contours
                    {
                        if (labels[mainindex] != labels[index]) np++;
                    }
                }
            }
            if (np > 1)
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                //img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    int numboundpix = cind; //int(contourx.size());
    for (int j = 0; j < numboundpix; j++)
    {
        for (int n = 0; n < 8; n++)
        {
            int x = contourx[j] + dx8[n];
            int y = contoury[j] + dy8[n];
            if (( x >= 0 && x < w ) && ( y >= 0 && y < h ))
            {
                int ind = rgb.computePosition( x, y );
                if (!istaken[ind])
                {
                    rgb.setPixel( ind, borderColor );
                }
            }
        }
    }
}



void enforceLabelConnectivity( const int* labels, //input labels that need to be corrected to remove stray labels
                               const int width,
                               const int height,
                               int*& nlabels, //new labels
                               int& numlabels, //the number of labels changes in the end if segments are removed
                               const int& K ) //the number of superpixels desired by the user
{
    const int dx4[4] = { -1, 0, 1, 0 };
    const int dy4[4] = { 0, -1, 0, 1 };

    const int sz = width * height;
    const int SUPSZ = sz / K;

    for (int i = 0; i < sz; i++) nlabels[i] = -1;
    int label( 0 );
    int* xvec = new int[sz];
    int* yvec = new int[sz];
    int oindex( 0 );
    int adjlabel( 0 ); //adjacent label
    for (int j = 0; j < height; j++)
    {
        for (int k = 0; k < width; k++)
        {
            if (0 > nlabels[oindex])
            {
                nlabels[oindex] = label;
                //--------------------
                // Start a new segment
                //--------------------
                xvec[0] = k;
                yvec[0] = j;
                //-------------------------------------------------------
                // Quickly find an adjacent label for use later if needed
                //-------------------------------------------------------
                {
                    for (int n = 0; n < 4; n++)
                    {
                        int x = xvec[0] + dx4[n];
                        int y = yvec[0] + dy4[n];
                        if (( x >= 0 && x < width ) && ( y >= 0 && y < height ))
                        {
                            int nindex = y * width + x;
                            if (nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                        }
                    }
                }

                int count( 1 );
                for (int c = 0; c < count; c++)
                {
                    for (int n = 0; n < 4; n++)
                    {
                        int x = xvec[c] + dx4[n];
                        int y = yvec[c] + dy4[n];

                        if (( x >= 0 && x < width ) && ( y >= 0 && y < height ))
                        {
                            int nindex = y * width + x;

                            if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }

                    }
                }
                //-------------------------------------------------------
                // If segment size is less then a limit, assign an
                // adjacent label found before, and decrement label count.
                //-------------------------------------------------------
                if (count <= SUPSZ >> 2)
                {
                    for (int c = 0; c < count; c++)
                    {
                        int ind = yvec[c] * width + xvec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }
    numlabels = label;

    if (xvec) delete [] xvec;
    if (yvec) delete [] yvec;
}



void SuperPixels( Image& rgb, int k, double M )
{
    //TODO: Converte a imagem para LAb.
    Image Lab =convertImageFromRGB2Lab(rgb);

    //TODO: Calcula o numero de pixels cada superpixel.
    int w = rgb.getW();
    int h = rgb.getH();
    int size = sqrt((w * h)/k);
    //Todo: Inicializa os os clusters.
    Cluster* clusters;
    
    int nk = initializeClusters( clusters, Lab, k );


    //TODO: aloca e inicializa labels.

    int* labels = (int*)malloc(sizeof(int)*w*h); //index = pixel, value = cluster associated 
    for(int i=0;i<(w*h);i++)
    {
        labels[i] = -1;
    }


    //TODO: Executa o algoritmo.

   // performSuperPixelsAlgorithm(Lab, clusters, labels, k, M);

    //    int* nlabels = new int[size];
    //    enforceLabelConnectivity( labels, w, h, nlabels, k, double(size ) / double( s * s ) );
    //    for (int i = 0; i < size; i++)
    //        labels[i] = nlabels[i];

    //if (nlabels)
    //    delete [] nlabels;

    //TODO: define as novas cores dos pixels.


    //TODO: Converte a imagem de volta.
    rgb = convertImageFromLAB2RGB(Lab);
    //Desenha os contornos. Deve passar a imagem em rgb e o vetor de labels.
    //drawContoursAroundSegments( rgb, labels );
}



/*
 * 
 */
int main( int argc, char** argv )
{

    Image l;
    if (l.readBMP( "AB_ufv_0675.bmp" ))
    {
        printf( "Leitura executada com sucesso\n" );
    }
    
    SuperPixels( l, 512, 20 );
    
    if (l.writeBMP( "AB_ufv_06752.bmp" ))
    {
        printf( "Escrita executada com sucesso\n" );
    }

    return 0;
}

