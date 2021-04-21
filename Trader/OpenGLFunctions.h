#pragma once

#include <stddef.h>

// APIENTRY and GLAPI definitions
#if defined(_WIN32) || defined(__WIN32__)
// MS Windows needs APIENTRY.
#ifndef APIENTRY
#define APIENTRY __stdcall
#endif

#ifndef GLAPI
#if defined(_MSC_VER)  // VC++
#define GLAPI __declspec(dllimport)
#else
#define GLAPI extern
#endif
#endif
#endif  // _WIN32 / __WIN32__

// Other platforms are simpler
#ifndef APIENTRY
#define APIENTRY
#endif
#ifndef GLAPI
#define GLAPI extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

using GLvoid = void;
using GLubyte = unsigned char;
using GLchar = char;
using GLboolean = unsigned char;
using GLushort = unsigned short;
using GLint = int;
using GLsizei = int;
using GLenum = unsigned int;
using GLuint = unsigned int;
using GLbitfield = unsigned int;
using GLfloat = float;
using GLsizeiptr = ptrdiff_t;
using GLintptr = ptrdiff_t;


#define GL_ARRAY_BUFFER 0x8892
#define GL_BLEND 0x0BE2
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_COLOR_BUFFER_BIT 0x00004000
#define GL_COMPILE_STATUS 0x8B81
#define GL_DEPTH_BUFFER_BIT 0x00000100
#define GL_DEPTH_TEST 0x0B71
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_FALSE 0
#define GL_FIRST_VERTEX_CONVENTION 0x8E4D
#define GL_FLOAT 0x1406
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_INFO_LOG_LENGTH 0x8B84
#define GL_LEQUAL 0x0203
#define GL_LINK_STATUS 0x8B82
#define GL_LINEAR 0x2601
#define GL_NO_ERROR 0
#define GL_ONE_MINUS_SRC_ALPHA 0x0303
#define GL_RENDERER 0x1F01
#define GL_RGBA 0x1908
#define GL_SRC_ALPHA 0x0302
#define GL_STATIC_DRAW 0x88E4
#define GL_TEXTURE0 0x84C0
#define GL_TEXTURE_2D 0x0DE1
#define GL_TEXTURE_BASE_LEVEL 0x813C
#define GL_TEXTURE_MAX_LEVEL 0x813D
#define GL_TEXTURE_MAG_FILTER 0x2800
#define GL_TEXTURE_MIN_FILTER 0x2801
#define GL_TEXTURE_WRAP_S 0x2802
#define GL_TEXTURE_WRAP_T 0x2803
#define GL_TRIANGLE_STRIP 0x0005
#define GL_UNSIGNED_BYTE 0x1401
#define GL_UNSIGNED_SHORT 0x1403
#define GL_VENDOR 0x1F00
#define GL_VERSION 0x1F02
#define GL_VERTEX_SHADER 0x8B31

GLAPI void APIENTRY glBindTexture(GLenum target, GLuint texture);
GLAPI void APIENTRY glBlendFunc(GLenum sfactor, GLenum dfactor);
GLAPI void APIENTRY glClear(GLbitfield mask);
GLAPI void APIENTRY glClearColor(GLfloat red, GLfloat green, GLfloat blue,
                                 GLfloat alpha);
GLAPI void APIENTRY glDeleteTextures(GLsizei n, const GLuint *textures);
GLAPI void APIENTRY glDepthFunc(GLenum func);
GLAPI void APIENTRY glDisable(GLenum cap);
GLAPI void APIENTRY glDrawArrays(GLenum mode, GLint first, GLsizei count);
GLAPI void APIENTRY glDrawElements(GLenum mode, GLsizei count, GLenum type,
                                   const void *indices);
GLAPI void APIENTRY glEnable(GLenum cap);
GLAPI void APIENTRY glFinish(void);
GLAPI void APIENTRY glFlush(void);
GLAPI void APIENTRY glGenTextures(GLsizei n, GLuint *textures);
GLAPI GLenum APIENTRY glGetError(void);
GLAPI const GLubyte *APIENTRY glGetString(GLenum name);
GLAPI void APIENTRY glTexImage2D(GLenum target, GLint level,
                                 GLint internalformat, GLsizei width,
                                 GLsizei height, GLint border, GLenum format,
                                 GLenum type, const void *pixels);
GLAPI void APIENTRY glTexParameteri(GLenum target, GLenum pname, GLint param);
GLAPI void APIENTRY glViewport(GLint x, GLint y, GLsizei width, GLsizei height);

using GL_CREATE_SHADER_PROC = GLuint(APIENTRY*)(GLenum type);

extern GL_CREATE_SHADER_PROC glCreateShader;

#ifdef __cplusplus
}
#endif

void InitGLPointers();