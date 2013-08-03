/*******************************************

  matrix.h

********************************************/

/*  'i'   integer
    length
    ' '   format 'c' character , 'f' double , 'i' integer
    array
    'e'   end of input
    example:
    print_vec  ( 'i' , nb , 'i' , b , 'e' ) ;
*/
void print_vec ( char format , ... ) ;

char** c_matrix2 ( int x , int y ) ;
void init_c_matrix2 ( char **a , int x , int y , char value ) ;
void free_c_matrix2 ( char **a , int x  ) ;
void print_c_matrix2 ( char **a , int x , int y ) ;
int Bin2Dec ( char *a , int len ) ;

