/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Cluster.h
 * Author: jcoelho
 *
 * Created on September 12, 2016, 10:16 AM
 */

#ifndef CLUSTER_H
#define CLUSTER_H
#include "Pixel.h"

class Cluster
{
public:
    /**
     * Construtor para o cluster.
     * @param lab - componentes Lab do cluster.
     * @param (x, y) - componentes da posicao do cluster.
     */
    Cluster( Pixel lab = Pixel( ), float x = 0.0, float y = 0.0 );

    /**
     * 
     * @return 
     */
    Pixel getPixel( ) const;

    /**
     * 
     * @return 
     */
    float getX( ) const;

    /**
     * 
     * @return 
     */
    float getY( ) const;

    /**
     * 
     * @param lab
     * @param x
     * @param y
     */
    void set( Pixel lab, float x, float y );

    /**
     * 
     * @param lab
     */
    void setPixel( Pixel lab );

    /**
     * 
     * @param x
     * @param y
     */
    void setPosition( float x, float y );

private:
    /**
     * Componente do Lab.
     */
    Pixel _lab;

    /**
     * Componentes da posicao.
     */
    float _x, _y;
};

#endif /* CLUSTER_H */

