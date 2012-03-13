#include<stdio.h>
#include "sbml/SBMLTypes.h"
#include "header.h"

boolean str_is_number(const char *str){
  int i = 0;
  int dot_count = 0;
  while(str[i] != '\0'){
    if((str[i] < 48 || str[i] > 57) && str[i] != 46){
      return 0;
    }else if(str[i] == 46){
      dot_count++;
      if(dot_count > 1){
        return 0;
      }
    }
    i++;
  }
  return 1;
}
