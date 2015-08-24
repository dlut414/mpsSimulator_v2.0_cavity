/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#pragma once
#include "Header.h"
#include "Draw_.h"
#include "Shader.h"
#include "Controller.h"
#include "../common/typedef/Vec3.h"
#include "../common/typedef/BBox.h"

namespace REN {

	class DrawSim : public Draw_ {
		typedef Vec3<double> vec;
	public:
		DrawSim(Controller* state) : Draw_(state) {}
		~DrawSim() {}

		void init() {
			Draw_::initVAO(Draw_::vao);
			Draw_::initVBO(Draw_::vbo[0]);
			Draw_::initVBO(Draw_::vbo[1]);
			Draw_::initVBO(Draw_::vbo[2]);
			Draw_::initVBO(Draw_::vbo[3]);
			Draw_::initVBO(Draw_::vbo[4]);
			Draw_::initVBO(Draw_::vbo[5]);
			initShader();
		}

		void draw(unsigned numVert, int* type, vec* vert, vec* vel, double* scal) {
			///clear framebuffer
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) printf("fbo_def not ready\n");
			glViewport(0, 0, Draw_::stateObj->u_width, Draw_::stateObj->u_height);
			glFrontFace(GL_CCW);
			glClearDepth(1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			
			///use program0
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) printf("0 not ready\n");
			glViewport(0, 0, Draw_::stateObj->u_width, Draw_::stateObj->u_height);
			glFrontFace(GL_CCW);
			glDepthFunc(GL_LESS);

			glClearDepth(1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glUseProgram(shaderObj.programID[0]);
			
			glUniform1i(shaderObj.intID[0], Draw_::stateObj->i_visFlag);
			glUniform1f(shaderObj.floatID[0], Draw_::stateObj->f_visRange);
			glUniformMatrix4fv(shaderObj.matrixID[1], 1, GL_FALSE, &(Draw_::stateObj->m_mvpInv[0][0]));
			glUniformMatrix4fv(shaderObj.matrixID[0], 1, GL_FALSE, &(Draw_::stateObj->m_mvp[0][0]));
			glUniformMatrix4fv(shaderObj.matrixID[1], 1, GL_FALSE, &(Draw_::stateObj->m_mvpInv[0][0]));

			glBindBuffer(GL_ARRAY_BUFFER, Draw_::vbo[0]);
			glBufferData(GL_ARRAY_BUFFER, numVert*sizeof(int), type, GL_STATIC_DRAW);
			glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);

			glBindBuffer(GL_ARRAY_BUFFER, Draw_::vbo[1]);
			glBufferData(GL_ARRAY_BUFFER, numVert*sizeof(vec), vert, GL_STATIC_DRAW);
			glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, (void*)0);

			glBindBuffer(GL_ARRAY_BUFFER, Draw_::vbo[2]);
			glBufferData(GL_ARRAY_BUFFER, numVert*sizeof(vec), vel, GL_STATIC_DRAW);
			glVertexAttribPointer(2, 3, GL_DOUBLE, GL_FALSE, 0, (void*)0);

			glBindBuffer(GL_ARRAY_BUFFER, Draw_::vbo[3]);
			glBufferData(GL_ARRAY_BUFFER, numVert*sizeof(double), scal, GL_STATIC_DRAW);
			glVertexAttribPointer(3, 1, GL_DOUBLE, GL_FALSE, 0, (void*)0);

			glEnableVertexAttribArray(0);
			glEnableVertexAttribArray(1);
			glEnableVertexAttribArray(2);
			glEnableVertexAttribArray(3);

			glDrawArrays(GL_POINTS, 0, numVert);

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			glDisableVertexAttribArray(2);
			glDisableVertexAttribArray(3);
		}

	public:
		GLfloat     r_bg_chessboard[36];

		GLuint bgTex;
		GLuint frontFaceTex;
		GLuint dirTex;
		GLuint dbo[2];

		GLuint volNormTex;
		GLuint volDataTex;
		GLuint cellInFluidTex;

		Shader shaderObj;

		BBox<float> box;
		GLuint ix, iy, iz;
		GLuint icx, icy, icz;

		void* volNorm;
		void* volData;
		void* cellInFluid;

	private:
		void initShader() {
			shaderObj.programID.push_back( shaderObj.LoadShader("./renderer/shader0/vertex.glsl", "./renderer/shader0/fragment.glsl") );
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "vMvp") );
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "fMvpInv") );

			shaderObj.intID.push_back(glGetUniformLocation(shaderObj.programID[0], "flag"));
			shaderObj.floatID.push_back(glGetUniformLocation(shaderObj.programID[0], "range"));
		}
	};

}
