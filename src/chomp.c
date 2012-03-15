#include "libsbmlsim/libsbmlsim.h"

void chomp(char *str)
{
  size_t pos = strlen(str);
  while (pos--) {
    if (str[pos] == '\n' || str[pos] == '\r')
      str[pos] = '\0';
    else
      break;
  }
}

