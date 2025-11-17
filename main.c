
#include <stdio.h>
#include <linalg.h>

int main() {
	Mat A = mat_from(3, 3,
		2.0, 1.0, 0.0,
		1.0, 1.0, 2.0,
		2.0, 0.0, 1.0
	);
	Mat B = mat_from(3, 1, 
		4.0,
		9.0,
		5.0
	);
	Mat result = mat_cramer(&A, &B);
	mat_print(&result);
}
