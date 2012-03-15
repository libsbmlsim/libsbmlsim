#include "libsbmlsim/libsbmlsim.h"

boolean str_is_number(const char *str)
{
  boolean isfloat = false;
  if (*str == '\0')
    return false;
  while (*str != '\0') {
    if (isdigit(*str)) {
      ;
    } else if (*str == '.') {
      if (isfloat)
        return false;
      isfloat = true;
    } else {
      return false;
    }
    str++;
  }
  return true;
}

