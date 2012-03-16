#include <stdarg.h>
#include <stdio.h>

void prg_printf(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vfprintf(stdout, fmt, args);
  va_end(args);
}
