#include "common.hh"

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char *filename,const char* mode) {
    FILE *fp=fopen(filename,mode);
    if(fp==NULL) {
        fprintf(stderr,"Can't open file '%s'\n",filename);
        exit(1);
    }
    return fp;
}

/** \brief Reads data from a file and checks that the operation was successful.
 *
 * Reads data from a file and checks the return value to ensure that the
 * operation was successful. If not successful, it prints an error message and
 * exits.
 * \param[in] ptr the memory to write to.
 * \param[in] size the size of each element to be read.
 * \param[in] count the number of elements to be read.
 * \param[in] fp the file handle to read from.
 * \param[in] p a description of what is being read, used in the error message
 *              if the operation isn't successful. */
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p) {
    if(fread(ptr,size,count,fp)!=count) {
        fprintf(stderr,"twin: can't read %s from file\n",p);
        exit(1);
    }
}

/** \brief Shifts the position indicator for an open file.
 *
 * Shifts the position indicator for an open file and checks the return value
 * to ensure that the operation was successful. If not successful, it prints an
 * error message and exits.
 * \param[in] fp the file handle to use.
 * \param[in] offset the number of bytes to offset by.
 * \param[in] origin the position used as a reference for the offset. */
void safe_fseek(FILE *fp,long int offset,int origin) {
    if(fseek(fp,offset,origin)!=0) {
        fprintf(stderr,"twin: error shifting file position by %ld bytes\n",offset);
        exit(1);
    }
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void fatal_error(const char *p,int status) {
    fprintf(stderr,"twin: %s\n",p);
    exit(status);
}
