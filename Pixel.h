/* 
 * File:   Pixel.h
 * Author: jcoelho
 *
 * Created on September 11, 2016, 6:03 PM
 */

#ifndef PIXEL_H
#define PIXEL_H

class Pixel
{
public:
    /**
     * Construtor default.
     */
    Pixel( );

    /**
     * Sobrecarga do operador [] para acessar as informacoes dos pixels.
     * @param i - posicao a ser acessada.
     * @return - referencia para a posicao.
     */
    float operator[]( unsigned int i ) const;

    /**
     * Sobrecarga do operador [] para acessar as informacoes dos pixels.
     * @param i - posicao a ser acessada.
     * @return - referencia para a posicao.
     */
    float& operator[]( unsigned int i );

    /**
     * Realiza o produto escalar entre dois pixels.
     * @param p - pixel a ser multiplicado pelo pixel corrente.
     * @return - resulatado do produto escalar.
     */
    float operator*( const Pixel& p ) const;
    
    /**
     * Sobrecarga do operador -.
     * @param p - pixel a ser subtraido.
     * @return - resultado da operacao.
     */
    Pixel operator-( const Pixel& p ) const;


private:
    float _val[3];
};

#endif /* PIXEL_H */

