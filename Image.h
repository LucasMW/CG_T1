/* 
 * File:   Image.h
 * Author: jcoelho
 *
 * Created on September 11, 2016, 2:18 PM
 */

#ifndef IMAGE_H
#define IMAGE_H

#include "Pixel.h"

class Image
{
public:
    /**
     * Construtor default.
     */
    Image( );

    /**
     * Construtor que recebe a dimensao da imagem.
     * @param w - largura da imagem.
     * @param h - altura da imagem.
     */
    Image( int w, int h );

    /**
     * Construtor de copia.
     * @param orig - imagem a ser usada para construir uma nova.
     */
    Image( const Image& orig );

    /**
     * Destrutor.
     */
    virtual ~Image( );

    /**
     * Le uma imagem em formato BMP.
     * @param path - caminho do arquivo.
     * @return - true se o arquivo pode ser aberto com sucesso e false
     * caso contrario.
     */
    bool readBMP( const char *path );

    /**
     * Escreve a imagem em formato BMP.
     * @param path - caminho do arquivo a ser escrito.
     * @return - true se o arquivo pode ser escrito com sucesso false contrario.
     */
    bool writeBMP( const char *path );

    /**
     * Obtem o pixel da posicao (i, j);
     * @param (i, j) - posicao da imagem que deve ser obtido o novo pixel.
     * @return - novo pixel.
     */
    Pixel getPixel( int i, int j ) const;

    /**
     * Obtem a largura da imagem.
     * @return - largura da imagem.
     */
    int getW( ) const;

    /**
     * Obtem a altura da imagem.
     * @return - altura da imagem.
     */
    int getH( ) const;

    /**
     * Define um novo pixel na imagem.
     * @param (i, j) - posicao na imagem para definir novo pixel.
     * @param p - novo pixel a ser definido.
     */
    void setPixel( int i, int j, Pixel p );

    /**
     * Define um pixel em uma posicao especifica.
     * @param position - posicao a ser definido o pixel.
     * @param p - valor a ser definido.
     */
    void setPixel( int position, Pixel p );

    /**
     * Converte de rgb para XYZ.
     * @param rgb - pixel a ser convertido.
     * @return - pixel convertido.
     */
    static Pixel rgbToXYZ( Pixel rgb );

    /**
     * Converte de XYZ para Lab.
     * @param XYZ - pixel a ser convertido.
     * @return - pixel convertido.
     */
    static Pixel XYZToLab( Pixel XYZ );

    /**
     * Converte de XYZ para rgb.
     * @param XYZ - pixel a ser convertido.
     * @return - pixel convertido.
     */
    static Pixel XYZTorgb( Pixel XYZ );

    /**
     * Converte de Lab para XYZ.
     * @param Lab - pixel a ser convertido.
     * @return - pixel convertido.
     */
    static Pixel LabToXYZ( Pixel Lab );

    /**
     * Calcula a posicao no vetor dada a posicao (i,j) a imagem.
     * @param (i,j) - posicao na imagem.
     * @return - posicao no vetor.
     */
    int computePosition( int i, int j ) const;

    /**
     * Operador de atribuicao.
     * @param img - imagem a ser copiada.
     * @return - referencia para a nova imagem.
     */
    Image& operator=( const Image& img );
private:
    /**
     * Dimensoes da imagem.
     */
    int _w, _h;

    /**
     * Pixels da imagem armazenados em sequencia.
     */
    Pixel *_data;
};

#endif /* IMAGE_H */

