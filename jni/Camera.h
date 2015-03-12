#ifndef _CAMERA
#define _CAMERA

#include <Structures.h>
//#include <Vec3.h>

class Camera {
private:
	float proj_w, proj_h;
public:
	Vertex pos, dir, up;

	float near, fov_hor;
	float res_x, res_y;

	float x_step_norm, y_step_norm;
	Vertex x_step, y_step, proj_origin;


	Camera(): near(10.0f), fov_hor(1.85f) {
	}

	void updateProjVectors() {
		cross(up, dir, x_step);
		cross(x_step, dir, y_step);

		normalize(x_step);
		normalize(y_step);

		x_step *= x_step_norm;
		y_step *= y_step_norm;

		proj_origin.x = pos.x + near * dir.x - x_step.x * (res_x * 0.5f) - y_step.x * (res_y * 0.5f);
		proj_origin.y = pos.y + near * dir.y - x_step.y * (res_x * 0.5f) - y_step.y * (res_y * 0.5f);
		proj_origin.z = pos.z + near * dir.z - x_step.z * (res_x * 0.5f) - y_step.z * (res_y * 0.5f);

		//proj_origin = pos + near * dir - x_step * (res_x * 0.5f) - y_step * (res_y * 0.5f);
	}

	void setup(const int &_res_x, const int &_res_y, const float &_fov_hor,
			const Vertex &_pos, const Vertex &_dir, const Vertex &_up) {
		res_x = (float)_res_x;
		res_y = (float)_res_y;
		fov_hor = _fov_hor;

		memcpy(&pos, &_pos, sizeof(Vertex));
		memcpy(&dir, &_dir, sizeof(Vertex));
		memcpy(&up, &_up, sizeof(Vertex));

		//normalize(pos);
		normalize(dir);
		normalize(up);

		// Calculate the pixel steps.
		proj_w = 2.0f * near * tan(fov_hor * 0.5f);
		proj_h = proj_w * res_y / res_x;

		x_step_norm = proj_w / res_x;
		y_step_norm = proj_h / res_y;

		updateProjVectors();
		/*
		LOGI("---------------------------------------");
		LOGI("Camera set up");
		LOGI("pos: %f, %f, %f", pos.x, pos.y, pos.z);
		LOGI("dir: %f, %f, %f", dir.x, dir.y, dir.z);
		LOGI("up: %f, %f, %f", up.x, up.y, up.z);
		LOGI("proj_w, proj_h: %f, %f", proj_w, proj_h);
		LOGI("x_step_norm, y_step_norm: %f, %f", x_step_norm, y_step_norm);
		LOGI("proj_origin: %f, %f, %f", proj_origin.x, proj_origin.y, proj_origin.z);
		LOGI("x_step: %f, %f, %f", x_step.x, x_step.y, x_step.z);
		LOGI("y_step: %f, %f, %f", y_step.x, y_step.y, y_step.z);
		*/
	}

	void changePos(const Vertex &_pos) {
		memcpy(&pos, &_pos, sizeof(Vertex));
		//normalize(pos);
		updateProjVectors();
	}

	void changeDir(const Vertex &_dir) {
		memcpy(&dir, &_dir, sizeof(Vertex));
		normalize(dir);
		updateProjVectors();
	}

	void changeResolution(const int &_res_x, const int &_res_y) {
		res_x = (float)_res_x;
		res_y = (float)_res_y;
		setup(res_x, res_y, fov_hor, pos, dir, up);

		//updateProjVectors();
	}
};

#endif
