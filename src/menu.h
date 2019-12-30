
/*
 
IMPORTANTE: Cuando el menu ya no sea util, ultilizar siempre la funcion finalizarMenu
 
*/
extern "C"{
#ifndef MENU_H
#define MENU_H
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
/*
* Para compilar en UNIX desactivar la definicion de WINDOWS y activar la de UNIX
* Para otros SO desactivar ambas.
*/
 /*
#ifndef WINDOWS
    #define WINDOWS
#endif
*/

#ifndef UNIX
    #define UNIX 
#endif

 
#define MENU_OK                   0L
#define MENU_ERROR                1L
#define MENU_MEMORIA_INSUFICIENTE 2L
 
#define MENU_MODO_NORMAL      0L
#define MENU_MODO_ULTIMA_CERO 1L
 
struct Menu
{
    char*         titulo;
    char*         texto;
    char**        opciones;
    int           numOpciones; /* numero de opciones */
    int           maxlen;
    unsigned long lFlags;
};
typedef struct Menu Menu;
 
void inicializarMenu(Menu* menu,char* titulo, char* texto, char* opciones[]);
int mostrarMenu(Menu* menu,int modo);
void finalizarMenu(Menu* menu);
 
#endif /* MENU_H */
}