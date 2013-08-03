    /* message.h */
    /*************/


#ifndef __MESSAGE_H__
#define __MESSAGE_H__


typedef void (*FailFun)(char *, ...);
 

void fatal (char format_str[], ... );


/*
void warning (char format_str[], ... );
*/

void message (char format_str[], ... );


void print (char format_str[], ... );


#endif /* __MESSAGE_H__ */
