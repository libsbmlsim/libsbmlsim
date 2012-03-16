#include <stdarg.h>
#include <stdio.h>

void dbg_printf(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}
