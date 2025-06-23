#ifndef _BUFFER_H_
#define _BUFFER_H_

#define SERIAL_MAX_NUM_BYTES 131072
#define MAX_STRING_BUFF SERIAL_MAX_NUM_BYTES

extern char string_buff[MAX_STRING_BUFF];

size_t getline_string(char **lineptr, int *n, int delimiter, FILE *fp);


#endif