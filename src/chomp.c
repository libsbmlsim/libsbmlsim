#include "libsbmlsim/libsbmlsim.h"

void chomp(char *str){
  int i;
  i = 0;
  while(1){
    if(str[i] == '\n' || str[i] == '\0'){
      str[i] = '\0';
      break;
    }
    i++;
  }
}
