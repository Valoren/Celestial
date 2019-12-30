#include <stdio.h>
#include <stdlib.h>
 
#include "menu.h"
#include "menu.cpp" /* solo si no se va crear un proyecto*/
 
int main(int argc, char *argv[])
{
    Menu menu;
    char *opciones[] = {"Continuar","Salir",""}; /* la ultima opcion siempre es una cadena vacia*/
    int opcion;
 
    inicializarMenu(&menu,"MENU DE PRUEBA","Escoger una opcion:",opciones);
 
    do{
        opcion=mostrarMenu(&menu,MENU_MODO_ULTIMA_CERO);
 
        switch(opcion)
        {
            case 1:
                printf("Has escogido continuar.\n");
                break;
            case 0:
                printf("Has escogido salir.\n");
                break;
        }
    }while(opcion!=0);
 
    /* ¡¡¡¡¡IMPORTANTE!!!!! */
    finalizarMenu(&menu);
 
    return 0;
}