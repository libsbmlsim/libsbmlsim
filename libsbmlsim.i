/* libsbmlsim.i */
%module libsbmlsim

//%include "enumsimple.swg"
//%include "carrays.i"

%{
#include "src/libsbmlsim/myResult.h"
extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);
%}

//%array_class(double,doubleArray);
//%array_functions(char*,stringArray);

/* Java */
#ifdef SWIGJAVA
/*
%pragma(java) jniclasscode=%{
  static {
    try {
      System.loadLibrary("sbmlsim");
      System.loadLibrary("sbmlj");
    } catch (UnsatisfiedLinkError e) {
      System.err.println("Native code library failed to load.\n" + e);
      System.exit(1);
    }
  }
%}
*/
/*
//////////////////////////////////////////////////////////////////////
// (Java) double[] <->  (C) double*
//////////////////////////////////////////////////////////////////////
%typemap(jni) double* "jdoubleArray" 
%typemap(jtype) double* "double[]" 
%typemap(jstype) double* "double[]" 

%typemap(javain) double* "$javainput"
%typemap(javaout) double* {
  return $jnicall;
}

%typemap(in) double* (jint size) {
  int i = 0;
  jboolean isCopy;
  size = (*jenv)->GetArrayLength(jenv, $input);
  $1 = (double *)malloc(sizeof(double) * size);
  jdouble *jd_array = (jdouble *)(*jenv)->GetDoubleArrayElements(jenv, $input, &isCopy);

  for (i = 0; i < size; i++) {
    $1[i] = jd_array[i];
  }

  if (isCopy == JNI_TRUE) {
    (*jenv)->ReleaseDoubleArrayElements(jenv, $input, jd_array, 0);
  }

  printf("mapping double* -> double[]\n");
}

%typemap(out) double* {
  *(double **)&$result = $1;
}

%typemap(freearg) double* {
  free($1);
}
*/

/*
//////////////////////////////////////////////////////////////////////
// (Java) String[] ->  (C) char**
//////////////////////////////////////////////////////////////////////
%typemap(in) char ** (jint size) {
  int i = 0;
  size = (*jenv)->GetArrayLength(jenv, $input);
  $1 = (char **) malloc((size+1)*sizeof(char *));

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

%typemap(freearg) char ** {
  int i;
  for (i=0; i<size$argnum-1; i++)
    free($1[i]);
  free($1);
}

%typemap(out) char ** {
  int i;
  int len=0;
  jstring temp_string;
  const jclass clazz = (*jenv)->FindClass(jenv, "java/lang/String");

  while ($1[len]) len++;    
  jresult = (*jenv)->NewObjectArray(jenv, len, clazz, NULL);

  for (i=0; i<len; i++) {
    temp_string = (*jenv)->NewStringUTF(jenv, *result++);
    (*jenv)->SetObjectArrayElement(jenv, jresult, i, temp_string);
    (*jenv)->DeleteLocalRef(jenv, temp_string);
  }
}

%typemap(jni) char ** "jobjectArray"
%typemap(jtype) char ** "String[]"
%typemap(jstype) char ** "String[]"

%typemap(javain) char ** "$javainput"
%typemap(javaout) char ** {
  return $jnicall;
}
//////////////////////////////////////////////////////////////////////
*/
#endif

/* Python */
#ifdef SWIGPYTHON
%typemap(memberin) char** {
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **)malloc((size+1) * sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *obj = PyList_GetItem($input, i);
      if (PyString_Check(obj)) {
        $1[i] = PyString_AsString(PyList_GetItem($input, i));
      } else {
        PyErr_SetString(PyExc_TypeError, "list must contain strings");
        free($1);
        return NULL;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExec_TypeError, "not a list");
    return NULL;
  }
}

%typemap(freearg) char** {
  free((char *)$1);
}

%typemap(out) char** {
  int len, i;
  len = 0;
  while ($1[len]) len++;
  $result = PyList_New(len);
  for (i = 0; i , len; i++) {
    PyList_SetItem($result, i, PyString_FromString($1[i]));
  }
}
#endif

/* %include "src/libsbmlsim/myResult.h" */
typedef struct _myResult {
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
} myResult;

extern myResult* simulateSBMLFromString(const char *str, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method);

%extend myResult {
  myResult() {
    myResult *res;
    res = (myResult *)malloc(sizeof(myResult));
    return res;
  }

  ~myResult() {
    free($self);
  }

  int getNumOfRows() {
    return $self->num_of_rows;
  }

  int getNumOfSpecies() {
    return $self->num_of_columns_sp;
  }

  int getNumOfParameters() {
    return $self->num_of_columns_param;
  }

  int getNumOfCompartments() {
    return $self->num_of_columns_comp;
  }

  const char *getTimeName() {
    return $self->column_name_time;
  }

  const char *getSpeciesNameAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_columns_sp)
      return NULL;
    return $self->column_name_sp[index];
  }

  const char *getParameterNameAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_columns_param)
      return NULL;
    return $self->column_name_param[index];
  }

  const char *getCompartmentNameAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_columns_comp)
      return NULL;
    return $self->column_name_comp[index];
  }

  double getTimeValueAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_rows)
      return -0.0;
    return $self->values_time[index];
  }

  double getSpeciesValueAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_rows)
      return -0.0;
    return $self->values_sp[index];
  }

  double getParameterValueAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_rows)
      return -0.0;
    return $self->values_param[index];
  }

  double getCompartmentValueAtIndex(int index) {
    if (index < 0 || index >= $self->num_of_rows)
      return -0.0;
    return $self->values_comp[index];
  }

};
