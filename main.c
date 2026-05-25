#include <stdio.h>
#include <linalg.h>

int main() {
	Mat A = mat_from(4, 4,
		0.0f, 0.5f, 0.0f, 0.2f,
		0.3f, 0.0f, 0.7f, 0.0f,
		0.0f, 0.5f, 0.0f, 0.0f,
		0.4f, 0.0f, 0.0f, 0.0f
	);
	Mat DM = mat_from(4, 4,
		0.7f, 0.5f, 0.0f, 0.2f,
		0.3f, 1.0f, 0.7f, 0.0f,
		0.0f, 0.5f, 5.0f, 0.0f,
		0.4f, 0.0f, 0.0f, 4.0f
	);
	Mat B = mat_diag(&DM);
	Mat C = mat_scale(&B, 1.0f/3);
	Mat D = mat_add(&A, &C);

	mat_print(&D);

	mat_free(&A);	
	mat_free(&B);	
	mat_free(&C);	
	mat_free(&D);	
}
