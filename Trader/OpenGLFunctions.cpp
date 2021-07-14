#include "OpenGLFunctions.h"

#if defined(_WIN32) || defined(__WIN32__)

#undef APIENTRY
#include <windows.h>  //For wglGetProcAddress

#define GetGLProcAddress(name) wglGetProcAddress((LPCSTR)name)
#else

#define GetGLProcAddress(name) (*glXGetProcAddress)((const GLubyte*)name)
#endif

auto GetGLFuncAddress(const char* fname) -> void* {
  void* pret = (void*)GetGLProcAddress(fname);

  if (pret == (void*)-1 || pret == (void*)1 || pret == (void*)2 || pret == (void*)3) pret = (void*)0;

  return pret;
}

GL_CREATE_SHADER_PROC glCreateShader = nullptr;

void InitGLPointers() {
  if (!(glCreateShader = reinterpret_cast<GL_CREATE_SHADER_PROC>(GetGLProcAddress("glCreateShader")))) throw;
}