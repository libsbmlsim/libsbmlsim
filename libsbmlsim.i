/* libsbmlsim.i */
%module libsbmlsim

%include "enumsimple.swg"

%{
#include "src/libsbmlsim/myResult.h"
extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
%}

//////////////////////////////////////////////////////////////////////
// (Java) double[] ->  (C) double*                                                                  
//////////////////////////////////////////////////////////////////////
%typemap(jni)   double*  "jdoubleArray" 
%typemap(jtype) double*  "double[]" 
%typemap(jstype) double*  "double[]" 

%typemap(javain) double* "$javainput"
%typemap(javaout) double* {
  return $jnicall;
}

%typemap(in) double* (jint size) {
  int i=0;
  jboolean isCopy;
  size = (*jenv)->GetArrayLength(jenv, $input);
  $1 = (double*)malloc(sizeof(double) * size);
  jdouble* jd_array = (jdouble*)(*jenv)->GetDoubleArrayElements(jenv, $input, &isCopy);

  for (i = 0; i < size; i++) {
    $1[i] = jd_array[i];
  }

  if (isCopy == JNI_TRUE) {
    (*jenv)->ReleaseDoubleArrayElements(jenv, $input, jd_array, 0);
  }
}

%typemap(out) double* {
  *(double **)&$result = $1;
}

%typemap(freearg) double* {
  free($1);
}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// (Java) String[] ->  (C) char**
//////////////////////////////////////////////////////////////////////
/* This tells SWIG to treat char ** as a special case when used as a parameter
   in a function call */
%typemap(in) char ** (jint size) {
  int i = 0;
  size = (*jenv)->GetArrayLength(jenv, $input);
  $1 = (char **) malloc((size+1)*sizeof(char *));
  /* make a copy of each string */
  for (i = 0; i<size; i++) {
    jstring j_string = (jstring)(*jenv)->GetObjectArrayElement(jenv, $input, i);
    const char * c_string = (*jenv)->GetStringUTFChars(jenv, j_string, 0);
    $1[i] = malloc((strlen(c_string)+1)*sizeof(char));
    strcpy($1[i], c_string);
    (*jenv)->ReleaseStringUTFChars(jenv, j_string, c_string);
    (*jenv)->DeleteLocalRef(jenv, j_string);
  }
  $1[i] = 0;
}

/* This cleans up the memory we malloc'd before the function call */
%typemap(freearg) char ** {
  int i;
  for (i=0; i<size$argnum-1; i++)
    free($1[i]);
  free($1);
}

/* This allows a C function to return a char ** as a Java String array */
%typemap(out) char ** {
  int i;
  int len=0;
  jstring temp_string;
  const jclass clazz = (*jenv)->FindClass(jenv, "java/lang/String");

  while ($1[len]) len++;    
  jresult = (*jenv)->NewObjectArray(jenv, len, clazz, NULL);
  /* exception checking omitted */

  for (i=0; i<len; i++) {
    temp_string = (*jenv)->NewStringUTF(jenv, *result++);
    (*jenv)->SetObjectArrayElement(jenv, jresult, i, temp_string);
    (*jenv)->DeleteLocalRef(jenv, temp_string);
  }
}

/* These 3 typemaps tell SWIG what JNI and Java types to use */
%typemap(jni) char ** "jobjectArray"
%typemap(jtype) char ** "String[]"
%typemap(jstype) char ** "String[]"

/* These 2 typemaps handle the conversion of the jtype to jstype typemap type
   and vice versa */
%typemap(javain) char ** "$javainput"
%typemap(javaout) char ** {
  return $jnicall;
}
//////////////////////////////////////////////////////////////////////

/* %include "src/libsbmlsim/myResult.h" */
typedef _myResult myResult;
struct _myResult{
  int num_of_rows;
  int num_of_columns_sp;
  int num_of_columns_param;
  int num_of_columns_comp;
%immutable;
  char *column_name_time;
  char **column_name_sp;
  char **column_name_param;
  char **column_name_comp;
%mutable;
  double *values_time;
  double *values_sp;
  double *values_param;
  double *values_comp;
};

extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
