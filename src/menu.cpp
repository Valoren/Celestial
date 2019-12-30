#include "menu.h"
extern "C"{
void inicializarMenu(Menu* menu, char* titulo, char* texto, char* opciones[])
{
   menu->lFlags = MENU_OK;
 
   if(titulo == NULL)
      menu->titulo = NULL;
   else
   {
       if((menu->titulo = (char*) malloc(strlen(titulo) * sizeof(char) + 1)))
           strcpy(menu->titulo,titulo);
      else
      {
         menu->lFlags |= MENU_ERROR | MENU_MEMORIA_INSUFICIENTE;
         menu->texto = NULL;
      }
   }
 
   if(texto == NULL)
      menu->texto = NULL;
   else
   {
      if((menu->texto = (char*) malloc(strlen(texto) * sizeof(char) + 1)))
         strcpy(menu->texto,texto);
      else
      {
         menu->lFlags |= MENU_ERROR | MENU_MEMORIA_INSUFICIENTE;
         menu->texto = NULL;
 
         if(menu->titulo)
         {
            free(menu->titulo);
            menu->titulo = NULL;
         }
      }
   }
 
   if(opciones == NULL)
      menu->opciones = NULL;
   else
   {
      menu->numOpciones = -1;
 
      while(opciones[ ++(menu->numOpciones) ][0]);
 
      if((menu->opciones = (char**) malloc(menu->numOpciones * sizeof(char*))))
      {
         int i=0;
 
         menu->maxlen = 0;
 
         for(i=0 ; i < menu->numOpciones ; i++)
         {
            if
            (
               (menu->opciones[i] =
               (char*) malloc(strlen(opciones[i]) * sizeof(char) + 1))
            )
            {
               strcpy(menu->opciones[i],opciones[i]);
 
               if(strlen(menu->opciones[i]) > menu->maxlen)
                  menu->maxlen = strlen(menu->opciones[i]);
            }
            else
            {
               int j;
 
               for(j=0 ; j<i ; j++)
                  free(menu->opciones[j]);
 
               free(menu->opciones);
 
               menu->opciones = NULL;
               menu->numOpciones = 0;
 
               if(menu->texto)
               {
                  free(menu->texto);
                  menu->texto = NULL;
               }
 
               if(menu->titulo)
               {
                  free(menu->titulo);
                  menu->titulo = NULL;
               }
 
               menu->lFlags |= MENU_ERROR | MENU_MEMORIA_INSUFICIENTE;
            }
         }
      }
      else
      {
         menu->lFlags |= MENU_ERROR | MENU_MEMORIA_INSUFICIENTE;
         menu->opciones = NULL;
 
         if(menu->texto)
         {
            free(menu->texto);
            menu->texto = NULL;
         }
 
         if(menu->titulo)
         {
            free(menu->titulo);
            menu->titulo = NULL;
         }
      }
   }
}
 
int mostrarMenu(Menu* menu,int modo)
{
   int i=0, opcion=0;
 
   if(menu->numOpciones)
   {
       do{
 
          #ifdef WINDOWS
             system("CLS");
          #elif defined UNIX
             system("clear");
          #endif
 
          if(menu->titulo)
          {
             for(i=0 ; i<80 ; i++)
                 printf("=");
 
             printf("%*s\n", 40 + strlen(menu->titulo) / 2, menu->titulo);
 
             for(i=0 ; i<80 ; i++)
                 printf("=");
          }
 
          if(menu->texto)
              printf("%-*s%s\n\n",40 - menu->maxlen / 2 - 8, " ", menu->texto);
 
          for(i=0 ; i < menu->numOpciones ; i++)
          {
             if(modo == MENU_MODO_ULTIMA_CERO && i == menu->numOpciones - 1)
             {
                if(menu->numOpciones > 10)
                   printf("%*s 0. %-*s\n", 40 - menu->maxlen / 2 - 4, " ", 40 + menu->maxlen / 2,
                                            menu->opciones[i]);
                else
                   printf("%*s0. %-*s\n", 40 - menu->maxlen / 2 - 4, " ", 40 + menu->maxlen / 2,
                                           menu->opciones[i]);
             }
             else
             {
                if(modo != MENU_MODO_ULTIMA_CERO)
                {
                   if(menu->numOpciones >= 10)
                   {
                      printf("%*s%2d. %-*s\n", 40 - menu->maxlen / 2 - 4, " ", i+1,
                                                40 + menu->maxlen / 2, menu->opciones[i]);
                   }
                   else
                   {
                      printf("%*s%d. %-*s\n", 40 - menu->maxlen / 2 - 4, " ", i+1,
                                               40 + menu->maxlen / 2, menu->opciones[i]);
                   }
                }
                else
                {
                   if(menu->numOpciones > 10)
                   {
                      printf("%*s%2d. %-*s\n", 40 - menu->maxlen / 2 - 4, " ", i+1,
                                                40 + menu->maxlen / 2, menu->opciones[i]);
                   }
                   else
                   {
                      printf("%*s%d. %-*s\n", 40 - menu->maxlen / 2 - 4, " ", i+1,
                                               40 + menu->maxlen / 2, menu->opciones[i]);
                   }
                }
             }
          }
          printf("%-*s>", 40 - menu->maxlen / 2 - 8, " ");
 
         if(!scanf("%d", &opcion))
         {
            while(getchar() != '\n');
            scanf("%d", &opcion);
         }         
 
       }while
       (
          (
             (modo != MENU_MODO_ULTIMA_CERO) ? (opcion < 1) : (opcion < 0)
          )
          ||
          (
             (modo != MENU_MODO_ULTIMA_CERO) ?
                (opcion >  menu->numOpciones)
                :
                (opcion >= menu->numOpciones)
          )
       );
 
       return opcion;
   }
 
   return -1;
}
 
void finalizarMenu(Menu* menu)
{
   if(menu->texto)
   {
      free(menu->texto);
      menu->texto = NULL;
   }
 
   if(menu->opciones)
   {
       int i=0;
 
       for(i=0 ; i < menu->numOpciones ; i++)
       {
           if(menu->opciones[i])
           {
               free(menu->opciones[i]);
               menu->opciones[i] = NULL;
           }
       }
       free(menu->opciones);
       menu->opciones = NULL;
   }
 
   menu->numOpciones = 0;
}
}