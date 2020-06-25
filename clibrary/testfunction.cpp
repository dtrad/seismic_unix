#include "su.h"
int testfunction()
{
  int testvalue;

  if (!getparint("testvalue", &testvalue)) testvalue = 0;

  fprintf(stderr,"\n testvalue=%d\n", testvalue);

  return(testvalue);
}


