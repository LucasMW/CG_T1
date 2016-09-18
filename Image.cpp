/* 
 * File:   Image.cpp
 * Author: jcoelho
 * 
 * Created on September 11, 2016, 2:18 PM
 */

#include "Image.h"
#include "ImageLoader.h"
#include "Pixel.h"
#include <cstring>
#include <cassert>
#include <cmath>



Image::Image( )
{
    _w = _h = 0;
    _data = 0;
}



Image::Image( int w, int h )
{
    _w = w;
    _h = h;
    _data = new Pixel[w * h];
}



Image::Image( const Image& orig )
{
    delete[] _data;

    //Copia as informacoes da imagem.
    _w = orig._w;
    _h = orig._h;

    //Aloca memoria.
    _data = new Pixel[3 * _w * _h];

    //Copia informacoes.
    memcpy( _data, orig._data, _w * _h * sizeof (Pixel ) );
}



Image::~Image( )
{
    delete _data;
}



bool Image::readBMP( const char *path )
{
    float *rgb;

    //Le a imagem.
    ImageLoader l;
    l.readBMP( rgb, _w, _h, path );

    //Copia as informacoes para a imagem.
    delete _data;
    _data = new Pixel[_w * _h];

    for (int j = 0; j < _h; j++)
    {
        for (int i = 0; i < _w; i++)
        {
            int pos = computePosition( i, j );
            for (int k = 0; k < 3; k++)
            {
                _data[pos][k] = rgb[3 * pos + k];
            }
        }
    }
    delete rgb;
}



bool Image::writeBMP( const char *path )
{
    float *rgb = new float[3 * _w * _h];
    for (int j = 0; j < _h; j++)
    {
        for (int i = 0; i < _w; i++)
        {
            int pos = computePosition( i, j );
            for (int k = 0; k < 3; k++)
            {
                rgb[3 * pos + k] = _data[pos][k];
            }
        }
    }

    //Le a imagem.
    ImageLoader l;
    l.writeBMP( rgb, _w, _h, path );
    delete rgb;
}



Pixel Image::getPixel( int i, int j ) const
{
    int pos = computePosition( i, j );
    return _data[pos];
}



int Image::getW( ) const
{
    return _w;
}



int Image::getH( ) const
{
    return _h;
}



void Image::setPixel( int i, int j, Pixel p )
{
    int pos = computePosition( i, j );
    _data[pos] = p;
}



void Image::setPixel( int position, Pixel p )
{
    _data[position] = p;
}



Pixel Image::rgbToXYZ( Pixel rgb )
{
    if (rgb[0] > 0.04045)
    {
        rgb[0] = pow( ( rgb[0] + 0.055 ) / 1.055, 2.4 );
    }
    else
    {
        rgb[0] = rgb[0] / 12.92;
    }

    if (rgb[1] > 0.04045)
    {
        rgb[1] = pow( ( rgb[1] + 0.055 ) / 1.055, 2.4 );
    }
    else
    {
        rgb[1] = rgb[1] / 12.92;
    }

    if (rgb[2] > 0.04045)
    {
        rgb[2] = pow( ( rgb[2] + 0.055 ) / 1.055, 2.4 );
    }
    else
    {
        rgb[2] = rgb[2] / 12.92;
    }

    rgb[0] = rgb[0] * 100.0;
    rgb[1] = rgb[1] * 100.0;
    rgb[2] = rgb[2] * 100.0;

    Pixel XYZ;

    //Observer. = 2°, Illuminant = D65
    XYZ[0] = rgb[0] * 0.4124 + rgb[1] * 0.3576 + rgb[2] * 0.1805;
    XYZ[1] = rgb[0] * 0.2126 + rgb[1] * 0.7152 + rgb[2] * 0.0722;
    XYZ[2] = rgb[0] * 0.0193 + rgb[1] * 0.1192 + rgb[2] * 0.9505;

    return XYZ;
}



Pixel Image::XYZToLab( Pixel XYZ )
{
    XYZ[0] /= 95.047; //ref_X =  95.047   Observer= 2°, Illuminant= D65
    XYZ[1] /= 100.000; //ref_Y = 100.000
    XYZ[2] /= 108.883; //ref_Z = 108.883

    if (XYZ[0] > 0.008856)
    {
        XYZ[0] = pow( XYZ[0], ( 1.0 / 3.0 ) );
    }
    else
    {
        XYZ[0] = ( 7.787 * XYZ[0] ) + ( 16.0 / 116.0 );
    }

    if (XYZ[1] > 0.008856)
    {
        XYZ[1] = pow( XYZ[1], ( 1.0 / 3.0 ) );
    }
    else
    {
        XYZ[1] = ( 7.787 * XYZ[1] ) + ( 16.0 / 116.0 );
    }

    if (XYZ[2] > 0.008856)
    {
        XYZ[2] = pow( XYZ[2], ( 1.0 / 3.0 ) );
    }
    else
    {
        XYZ[2] = ( 7.787 * XYZ[2] ) + ( 16.0 / 116.0 );
    }

    Pixel Lab;
    Lab[0] = ( 116 * XYZ[1] ) - 16;
    Lab[1] = 500 * ( XYZ[0] - XYZ[1] );
    Lab[2] = 200 * ( XYZ[1] - XYZ[2] );

    return Lab;
}



Pixel Image::XYZTorgb( Pixel XYZ )
{
    XYZ[0] /= 100; //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
    XYZ[1] /= 100; //Y from 0 to 100.000
    XYZ[2] /= 100; //Z from 0 to 108.883

    Pixel rgb;

    rgb[0] = XYZ[0] * 3.2406 + XYZ[1] * -1.5372 + XYZ[2] * -0.4986;
    rgb[1] = XYZ[0] * -0.9689 + XYZ[1] * 1.8758 + XYZ[2] * 0.0415;
    rgb[2] = XYZ[0] * 0.0557 + XYZ[1] * -0.2040 + XYZ[2] * 1.0570;

    if (rgb[0] > 0.0031308)
    {
        rgb[0] = 1.055 * pow( rgb[0], ( 1 / 2.4 ) ) - 0.055;
    }
    else
    {
        rgb[0] = 12.92 * rgb[0];
    }

    if (rgb[1] > 0.0031308)
    {
        rgb[1] = 1.055 * pow( rgb[1], ( 1 / 2.4 ) ) - 0.055;
    }
    else
    {
        rgb[1] = 12.92 * rgb[1];
    }

    if (rgb[2] > 0.0031308)
    {
        rgb[2] = 1.055 * pow( rgb[2], ( 1 / 2.4 ) ) - 0.055;
    }
    else
    {
        rgb[2] = 12.92 * rgb[2];
    }
    return rgb;
}



Pixel Image::LabToXYZ( Pixel Lab )
{
    Pixel XYZ;
    XYZ[1] = ( Lab[0] + 16 ) / 116;
    XYZ[0] = Lab[1] / 500 + XYZ[1];
    XYZ[2] = XYZ[1] - Lab[2] / 200;

    if (pow( XYZ[1], 3 ) > 0.008856)
    {
        XYZ[1] = pow( XYZ[1], 3 );
    }
    else
    {
        XYZ[1] = ( XYZ[1] - 16.0 / 116.0 ) / 7.787;
    }

    if (pow( XYZ[0], 3 ) > 0.008856)
    {
        XYZ[0] = pow( XYZ[0], 3 );
    }
    else
    {
        XYZ[0] = ( XYZ[0] - 16.0 / 116.0 ) / 7.787;
    }
    if (pow( XYZ[2], 3 ) > 0.008856)
    {
        XYZ[2] = pow( XYZ[2], 3 );
    }
    else
    {
        XYZ[2] = ( XYZ[2] - 16.0 / 116.0 ) / 7.787;
    }

    XYZ[0] = 95.047 * XYZ[0]; //ref_X =  95.047     Observer= 2°, Illuminant= D65
    XYZ[1] = 100.000 * XYZ[1]; //ref_Y = 100.000
    XYZ[2] = 108.883 * XYZ[2]; //ref_Z = 108.883

    return XYZ;
}



int Image::computePosition( int i, int j ) const
{
    return j * _w + i;
}



Image& Image::operator=( const Image& img )
{
    delete [] _data;
    //Copia as informacoes da imagem.
    _w = img._w;
    _h = img._h;

    //Aloca memoria.
    _data = new Pixel[3 * _w * _h];

    //Copia informacoes.
    memcpy( _data, img._data, _w * _h * sizeof (Pixel ) );

    return *this;
}