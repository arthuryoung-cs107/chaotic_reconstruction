#include "MH_tools.hh"

void fseek_SAFE(FILE *fp,long int offset,int origin)
{
  if(fseek(fp,offset,origin)!=0)
  {
    printf("fseek_SAFE: error shifting file position by %ld bytes\n",offset);
    exit(1);
  }
}
void fread_SAFE(void *ptr,size_t size,size_t count,FILE *fp)
{
  if(fread(ptr,size,count,fp)!=count)
  {
    printf("fread_SAFE: can't read file\n");
    exit(1);
  }
}
